"""
MIMO RTU - SCADA HMI Dashboard (with NESO replay capture mode)
Reads serial telemetry from Feather 1 M6 sketch.
Displays live measurements, alarms, three-phase waveform reconstruction,
PDS compliance summary and a raw-hex DNP3 frame inspector.

NESO replay capture mode:
- Optional second serial connection to Feather 2 (DAC generator).
- "Start Replay" button issues REPLAY_START to Feather 2.
- Logs every Feather 1 DNP3 frame to a timestamped CSV, tagged with the
  current replay index and target frequency parsed from Feather 2's
  "REPLAY <i> <f>" tokens.
- On REPLAY_END, closes the CSV and renders a Matplotlib overlay
  (NESO reference vs measured f_Vb vs |error|) for sanity-checking the
  run before declaring it good.

Usage:  python SCADA_HMI.py                  (port selector dialog opens)
        python SCADA_HMI.py COM5             (Feather 1 only on COM5)
        python SCADA_HMI.py COM5 COM6        (Feather 1 = COM5, Feather 2 = COM6)

Architecture: a startup dialog lets
the operator pick the Feather's COM port from a dropdown (two dropdowns:
Feather 1 required, Feather 2 optional). Tkinter UI on the main thread,
two daemon reader threads polling 115 200 baud with automatic reconnection,
and a shared RTUData store. UI refreshes at 2 Hz via Tk.after().
"""

import tkinter as tk
from tkinter import ttk, messagebox
import serial
import serial.tools.list_ports
import threading
import time
import re
import math
import sys
import csv
import os
from datetime import datetime, timezone


# ─────────────────────────────── Data Store ───────────────────────────────
class RTUData:
    """Thread-shared snapshot of the most recent telemetry frame."""

    def __init__(self):
        # RMS values (volts / amps)
        self.Va_rms = 0.0
        self.Vb_rms = 0.0
        self.Vc_rms = 0.0
        self.Ib_rms = 0.0
        self.Ic_rms = 0.0
        # Per-phase frequency (Hz)
        self.f_Va = 0.0
        self.f_Vb = 0.0
        self.f_Vc = 0.0
        # Per-phase THD (%)
        self.THD_Va = 0.0
        self.THD_Vb = 0.0
        self.THD_Vc = 0.0
        # Inter-channel phase angles (degrees)
        self.Ph_VbIb = 0.0
        self.Ph_VcIc = 0.0
        # Phase-b power triangle
        self.Pb = 0.0
        self.Qb = 0.0
        self.Sb = 0.0
        self.PFb = 0.0
        # Phase-c power triangle
        self.Pc = 0.0
        self.Qc = 0.0
        self.Sc = 0.0
        self.PFc = 0.0
        # Alarm flags (set by firmware GGIO logic)
        self.freq_alarm = False
        self.uv_alarm = False
        self.thd_alarm = False
        # Telemetry metadata
        self.dnp3_bytes = 0
        self.frame_count = 0
        self.last_update = "Never"
        self.connected = False
        self.raw_frame = ""
        self.raw_lines = []

        # ── Replay state (driven by the Feather 2 reader) ──
        # connected            : Feather 2 USB link up
        # replay_active        : between REPLAY_START and REPLAY_END
        # replay_complete      : latched True when REPLAY_END seen (UI polls)
        # replay_index         : last index seen on "REPLAY <i> <f>" line
        # replay_target_f      : last target frequency from Feather 2 (Hz)
        # replay_total         : declared trajectory length (default 360)
        self.f2_connected = False
        self.replay_active = False
        self.replay_complete = False
        self.replay_index = -1
        self.replay_target_f = float('nan')
        self.replay_total = 360

        # Set by SCADAApp.start_replay(); cleared on close.
        self.replay_logger = None


# ───────────────────────────── Serial Parser ──────────────────────────────
_FLOAT_RE = r'[+-]?\d+\.?\d*'


def parse_float(text, key):
    """Pull a float out of an ASCII telemetry line of the form 'key: 1.234 V'."""
    pattern = re.escape(key) + r'[:\s]*(' + _FLOAT_RE + r')'
    m = re.search(pattern, text)
    return float(m.group(1)) if m else None


def parse_line(line, data):
    """Update `data` in-place from one ASCII line emitted by the M6 sketch."""
    line = line.strip()
    if not line:
        return

    # Debug ring buffer (last 50 lines)
    data.raw_lines.append(line)
    if len(data.raw_lines) > 50:
        data.raw_lines.pop(0)

    # DNP3 frame header: "[DNP3 FRAME 28 bytes]"
    if 'DNP3 FRAME' in line:
        m = re.search(r'(\d+)\s*bytes', line)
        if m:
            data.dnp3_bytes = int(m.group(1))
            data.frame_count += 1
            data.last_update = time.strftime("%H:%M:%S")
            # If we're in a replay run, log this frame to CSV. The parser
            # runs in the Feather 1 reader thread; ReplayLogger is internally
            # locked so writes from any thread are safe.
            logger = getattr(data, 'replay_logger', None)
            if logger is not None and data.replay_active:
                try:
                    logger.write_row(data)
                except Exception as e:
                    print(f"[parse_line] CSV write failed: {e}")
        data.raw_frame = ""

    # Hex frame bytes: lines starting "XX " (two hex digits then space)
    if re.match(r'^([0-9A-Fa-f]{2}\s)+[0-9A-Fa-f]{2}\s*$', line):
        data.raw_frame += line + " "

    # RMS values
    for key in ('Va_rms', 'Vb_rms', 'Vc_rms', 'Ib_rms', 'Ic_rms'):
        v = parse_float(line, key)
        if v is not None:
            setattr(data, key, v)

    # Frequencies
    for key in ('f_Va', 'f_Vb', 'f_Vc'):
        v = parse_float(line, key)
        if v is not None:
            setattr(data, key, v)

    # THD
    for key in ('THD_Va', 'THD_Vb', 'THD_Vc'):
        v = parse_float(line, key)
        if v is not None:
            setattr(data, key, v)

    # Phase angles
    v = parse_float(line, 'Ph_Vb/Ib')
    if v is not None:
        data.Ph_VbIb = v
    v = parse_float(line, 'Ph_Vc/Ic')
    if v is not None:
        data.Ph_VcIc = v

    # Phase-b power line: "Pb: 0.012 W Qb: 0.011 VAR Sb: 0.016 VA PFb: 0.770"
    if 'Pb:' in line and 'Qb:' in line:
        for key in ('Pb', 'Qb', 'Sb', 'PFb'):
            v = parse_float(line, key)
            if v is not None:
                setattr(data, key, v)

    # Phase-c power line
    if 'Pc:' in line and 'Qc:' in line:
        for key in ('Pc', 'Qc', 'Sc', 'PFc'):
            v = parse_float(line, key)
            if v is not None:
                setattr(data, key, v)

    # Alarm flags
    m = re.search(r'FREQ_alarm:\s*(\d)', line)
    if m:
        data.freq_alarm = (m.group(1) == '1')
    m = re.search(r'UV_alarm:\s*(\d)', line)
    if m:
        data.uv_alarm = (m.group(1) == '1')
    m = re.search(r'THD_alarm:\s*(\d)', line)
    if m:
        data.thd_alarm = (m.group(1) == '1')


# ───────────────────────── Serial Reader Thread ───────────────────────────
class SerialReader(threading.Thread):
    """Background thread: opens the port, reads ASCII lines, parses into `data`,
    auto-reconnects every 2 s on failure."""

    def __init__(self, port, baud, data):
        super().__init__(daemon=True)
        self.port = port
        self.baud = baud
        self.data = data
        self.running = True

    def run(self):
        while self.running:
            try:
                ser = serial.Serial(self.port, self.baud, timeout=1)
                self.data.connected = True
                print(f"[SerialReader] connected to {self.port} @ {self.baud}")
                while self.running:
                    try:
                        raw = ser.readline()
                        if not raw:
                            continue
                        line = raw.decode('utf-8', errors='ignore')
                        if line:
                            parse_line(line, self.data)
                    except serial.SerialException:
                        raise  # drop out to reconnect
                    except Exception as e:
                        # Don't kill the loop on a single bad line
                        print(f"[SerialReader] parse error: {e}")
            except Exception as e:
                self.data.connected = False
                print(f"[SerialReader] error: {e} - retrying in 2 s...")
                time.sleep(2)


# ─────────────────── Feather 2 Reader / Commander ────────────────────────
class Feather2Reader(threading.Thread):
    """Background thread for Feather 2 (DAC generator board).

    Inbound: parses "REPLAY <i> <f>" lines (one per second during a replay)
    and the terminating "REPLAY_END" token. Updates RTUData replay fields.

    Outbound: send_command("REPLAY_START") / send_command("REPLAY_ABORT")
    from the UI thread. A small lock guards the serial write so the
    background reader doesn't get a partial-write collision.
    """

    REPLAY_RE = re.compile(r'^REPLAY\s+(\d+)\s+(' + _FLOAT_RE + r')\s*$')

    def __init__(self, port, baud, data):
        super().__init__(daemon=True)
        self.port = port
        self.baud = baud
        self.data = data
        self.running = True
        self._ser = None
        self._ser_lock = threading.Lock()

    def send_command(self, cmd):
        """Send a single ASCII command line to Feather 2. Returns True on
        successful write, False if the port isn't open or the write failed."""
        with self._ser_lock:
            if self._ser is None or not self._ser.is_open:
                return False
            try:
                payload = (cmd.strip() + '\n').encode('ascii')
                self._ser.write(payload)
                self._ser.flush()
                print(f"[Feather2] sent: {cmd.strip()}")
                return True
            except Exception as e:
                print(f"[Feather2] send_command failed: {e}")
                return False

    def _handle_line(self, line):
        line = line.strip()
        if not line:
            return
        # "REPLAY <i> <f>" — one per second during replay
        m = self.REPLAY_RE.match(line)
        if m:
            try:
                self.data.replay_index = int(m.group(1))
                self.data.replay_target_f = float(m.group(2))
                if not self.data.replay_active:
                    # First REPLAY token after REPLAY_START: arm the active flag
                    self.data.replay_active = True
                    self.data.replay_complete = False
            except ValueError:
                pass
            return
        # Terminator
        if line == 'REPLAY_END':
            self.data.replay_active = False
            self.data.replay_complete = True
            print("[Feather2] REPLAY_END received")
            return
        # Optional: echo of acknowledgements, boot banners, errors -
        # surface in console for troubleshooting without polluting RTUData.
        if line.startswith(('REPLAY_', '#', '[')):
            print(f"[Feather2] {line}")

    def run(self):
        while self.running:
            try:
                with self._ser_lock:
                    self._ser = serial.Serial(self.port, self.baud, timeout=1)
                self.data.f2_connected = True
                print(f"[Feather2Reader] connected to {self.port} @ {self.baud}")
                while self.running:
                    try:
                        raw = self._ser.readline()
                        if not raw:
                            continue
                        line = raw.decode('utf-8', errors='ignore')
                        if line:
                            self._handle_line(line)
                    except serial.SerialException:
                        raise
                    except Exception as e:
                        print(f"[Feather2Reader] parse error: {e}")
            except Exception as e:
                self.data.f2_connected = False
                with self._ser_lock:
                    try:
                        if self._ser is not None:
                            self._ser.close()
                    except Exception:
                        pass
                    self._ser = None
                print(f"[Feather2Reader] error: {e} - retrying in 2 s...")
                time.sleep(2)


# ─────────────────────────── Replay CSV Logger ────────────────────────────
# Column order is the same one used by replay_metrics.py - keep them in sync.
CSV_COLUMNS = [
    "iso_time", "monotonic_s", "frame_count",
    "replay_index", "f_target_hz",
    "Va_rms", "Vb_rms", "Vc_rms", "Ib_rms", "Ic_rms",
    "f_Va", "f_Vb", "f_Vc",
    "THD_Va", "THD_Vb", "THD_Vc",
    "Ph_VbIb", "Ph_VcIc",
    "Pb", "Qb", "Sb", "PFb",
    "Pc", "Qc", "Sc", "PFc",
    "FREQ_alarm", "UV_alarm", "THD_alarm",
    "dnp3_bytes",
]


class ReplayLogger:
    """Streams one CSV row per Feather 1 DNP3 frame during a replay run.

    Lifecycle:
        log = ReplayLogger.start(out_dir)        # opens file, writes header
        log.write_row(data)                       # called per F1 frame
        log.close()                               # closes file
        log.path                                  # final CSV path

    Thread-safety: writes are protected by an internal lock because rows
    are appended from the SCADAApp.update_ui callback (UI thread), but the
    file may also be closed by the Feather 2 reader on REPLAY_END.
    """

    def __init__(self, path):
        self.path = path
        self._fh = None
        self._writer = None
        self._lock = threading.Lock()
        self._rows_written = 0
        self._t0 = time.monotonic()

    @classmethod
    def start(cls, out_dir="."):
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        path = os.path.join(out_dir, f"replay_{stamp}.csv")
        logger = cls(path)
        logger._fh = open(path, "w", newline="", encoding="utf-8")
        logger._writer = csv.DictWriter(logger._fh, fieldnames=CSV_COLUMNS)
        logger._writer.writeheader()
        logger._fh.flush()
        print(f"[ReplayLogger] opened {path}")
        return logger

    def write_row(self, data):
        with self._lock:
            if self._fh is None or self._fh.closed:
                return
            row = {
                "iso_time": datetime.now(timezone.utc).isoformat(timespec='milliseconds'),
                "monotonic_s": f"{time.monotonic() - self._t0:.3f}",
                "frame_count": data.frame_count,
                "replay_index": data.replay_index,
                "f_target_hz": (f"{data.replay_target_f:.4f}"
                                if not math.isnan(data.replay_target_f) else ""),
                "Va_rms": f"{data.Va_rms:.4f}",
                "Vb_rms": f"{data.Vb_rms:.4f}",
                "Vc_rms": f"{data.Vc_rms:.4f}",
                "Ib_rms": f"{data.Ib_rms:.5f}",
                "Ic_rms": f"{data.Ic_rms:.5f}",
                "f_Va": f"{data.f_Va:.4f}",
                "f_Vb": f"{data.f_Vb:.4f}",
                "f_Vc": f"{data.f_Vc:.4f}",
                "THD_Va": f"{data.THD_Va:.3f}",
                "THD_Vb": f"{data.THD_Vb:.3f}",
                "THD_Vc": f"{data.THD_Vc:.3f}",
                "Ph_VbIb": f"{data.Ph_VbIb:.3f}",
                "Ph_VcIc": f"{data.Ph_VcIc:.3f}",
                "Pb": f"{data.Pb:.5f}",
                "Qb": f"{data.Qb:.5f}",
                "Sb": f"{data.Sb:.5f}",
                "PFb": f"{data.PFb:.4f}",
                "Pc": f"{data.Pc:.5f}",
                "Qc": f"{data.Qc:.5f}",
                "Sc": f"{data.Sc:.5f}",
                "PFc": f"{data.PFc:.4f}",
                "FREQ_alarm": int(bool(data.freq_alarm)),
                "UV_alarm": int(bool(data.uv_alarm)),
                "THD_alarm": int(bool(data.thd_alarm)),
                "dnp3_bytes": data.dnp3_bytes,
            }
            self._writer.writerow(row)
            self._fh.flush()
            self._rows_written += 1

    def close(self):
        with self._lock:
            if self._fh is not None and not self._fh.closed:
                self._fh.close()
                print(f"[ReplayLogger] closed {self.path} "
                      f"({self._rows_written} rows)")


# ─────────────────────────── Overlay Plot Helper ──────────────────────────
def load_neso_reference(csv_path):
    """Load the NESO reference trajectory from a CSV with columns
    'index,frequency_hz'. Returns (indices, frequencies) or (None, None)
    if the file doesn't exist or is malformed."""
    if not os.path.isfile(csv_path):
        return None, None
    try:
        idx, freq = [], []
        with open(csv_path, newline='', encoding='utf-8') as fh:
            r = csv.DictReader(fh)
            for row in r:
                idx.append(int(row['index']))
                freq.append(float(row['frequency_hz']))
        return idx, freq
    except Exception as e:
        print(f"[load_neso_reference] failed: {e}")
        return None, None


def generate_overlay_plot(csv_path, neso_csv_path="neso_trajectory.csv"):
    """Render the sanity-check overlay for one captured replay CSV.

    Three traces: NESO reference (truth), measured f_Vb (clean DAC2 channel),
    and Feather 2's emitted target frequency. Annotates max |error| in mHz.
    Saves PNG next to the CSV. Returns the PNG path, or None on failure.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")  # safe from worker contexts and headless boxes
        import matplotlib.pyplot as plt
    except ImportError:
        print("[overlay] matplotlib not installed - skipping plot. "
              "pip install matplotlib to enable.")
        return None

    # Load captured CSV
    rows = []
    with open(csv_path, newline='', encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            try:
                ridx = int(row['replay_index'])
                if ridx < 0:
                    continue
                rows.append({
                    'i':       ridx,
                    'f_Vb':    float(row['f_Vb']),
                    'f_target': float(row['f_target_hz']) if row['f_target_hz'] else float('nan'),
                })
            except (ValueError, KeyError):
                continue
    if not rows:
        print("[overlay] no replay-tagged rows in CSV - nothing to plot.")
        return None

    rows.sort(key=lambda r: r['i'])
    indices = [r['i'] for r in rows]
    f_meas  = [r['f_Vb'] for r in rows]
    f_targ  = [r['f_target'] for r in rows]

    # Load NESO ground truth
    neso_i, neso_f = load_neso_reference(neso_csv_path)
    using_neso = neso_i is not None

    # Build per-index NESO reference aligned to captured indices
    if using_neso:
        neso_map = dict(zip(neso_i, neso_f))
        ref_f = [neso_map.get(i, float('nan')) for i in indices]
        ref_label = "NESO reference (truth)"
    else:
        # Fall back to Feather 2 emitted targets - weaker sanity check
        ref_f = f_targ
        ref_label = "Feather 2 emitted target (fallback - neso_trajectory.csv not found)"

    # Compute |error|
    err_mhz = [(abs(m - r) * 1000.0) if not math.isnan(r) else float('nan')
               for m, r in zip(f_meas, ref_f)]
    finite_err = [e for e in err_mhz if not math.isnan(e)]
    max_err = max(finite_err) if finite_err else float('nan')
    mean_err = sum(finite_err) / len(finite_err) if finite_err else float('nan')

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6.5), sharex=True,
                                    gridspec_kw={'height_ratios': [3, 1]})

    ax1.plot(indices, ref_f, color='#1f77b4', linewidth=1.8,
             label=ref_label)
    ax1.plot(indices, f_meas, color='#d62728', linewidth=1.2,
             linestyle='--', label='Measured f_Vb (clean DAC2 channel)')
    if using_neso:
        ax1.plot(indices, f_targ, color='#2ca02c', linewidth=0.8,
                 linestyle=':', label='Feather 2 emitted target')
    ax1.axhline(49.5, color='#888888', linestyle=':', linewidth=0.7)
    ax1.axhline(50.5, color='#888888', linestyle=':', linewidth=0.7)
    ax1.text(indices[-1], 49.5, ' 49.5 Hz Grid Code',
             color='#666666', fontsize=8, va='bottom', ha='right')
    ax1.set_ylabel("Frequency (Hz)")
    ax1.set_title(f"NESO replay sanity check — {os.path.basename(csv_path)}\n"
                  f"max |error| = {max_err:.1f} mHz, "
                  f"mean |error| = {mean_err:.1f} mHz")
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='lower right', fontsize=9)

    # ── Error panel: |error| as percentage of nominal 50 Hz ────────────────
    # 500 mHz = 1 % of 50 Hz; written as /50*100 to make the "% of nominal"
    # meaning obvious to a reader following the code.
    err_pct = [
        (abs(m - r) / 50.0 * 100.0) if not math.isnan(r) else float('nan')
        for m, r in zip(f_meas, ref_f)
    ]
    finite_err_pct = [e for e in err_pct if not math.isnan(e)]
    max_err_pct = max(finite_err_pct) if finite_err_pct else float('nan')

    # Five performance regions, for display only. This is an earlier partition:
    # the paper's Section 3.1 partition ends the nadir region at index 200, not
    # 180 (nadir and early recovery 120-200, late recovery 201-300). Shading is
    # drawn first with zorder=0 so the error trace paints on top.
    REGIONS = [
        ('Pre-event nominal',   0,   45,  '#e8f0fe'),
        ('RoCoF descent',       46,  119, '#fde8e8'),
        ('Nadir',               120, 180, '#fdead7'),
        ('Late recovery',       181, 300, '#e8f4ea'),
        ('Settled',             301, 359, '#f0eef7'),
    ]
    for _label, i_start, i_end, colour in REGIONS:
        ax2.axvspan(i_start, i_end, alpha=0.35, color=colour, zorder=0)

    # Error trace on top of the shading
    ax2.plot(indices, err_pct, color='#9467bd', linewidth=1.0, zorder=3)

    # IEC 61400-21 compliance reference lines (as % of nominal 50 Hz)
    ax2.axhline(0.02, color='#888888', linestyle=':', linewidth=0.7, zorder=2)
    ax2.text(indices[-1], 0.02, ' IEC 61400-21 Class I (±0.02 %)',
             color='#666666', fontsize=8, va='bottom', ha='right')
    ax2.axhline(0.20, color='#aaaaaa', linestyle=':', linewidth=0.7, zorder=2)
    ax2.text(indices[-1], 0.20, ' Class II (±0.2 %)',
             color='#666666', fontsize=8, va='bottom', ha='right')

    # Lock the y-axis before labelling so region labels land at a fixed spot.
    y_top = max(max_err_pct * 1.1, 0.25) if not math.isnan(max_err_pct) else 0.25
    ax2.set_ylim(0.0, y_top)

    # Region labels near the top of the panel
    for label, i_start, i_end, _colour in REGIONS:
        mid = (i_start + i_end) / 2.0
        ax2.text(mid, y_top * 0.92, label,
                 fontsize=7, ha='center', va='top', color='#333333',
                 style='italic', alpha=0.85, zorder=4)

    ax2.set_xlabel("Replay index (1 s per step)")
    ax2.set_ylabel("|error|  (% of nominal 50 Hz)")
    ax2.grid(True, alpha=0.3, zorder=1)

    fig.tight_layout()
    png_path = os.path.splitext(csv_path)[0] + ".png"
    fig.savefig(png_path, dpi=140)
    plt.close(fig)
    print(f"[overlay] saved {png_path} "
          f"(max |err| {max_err:.1f} mHz, mean {mean_err:.1f} mHz)")
    return png_path


# ──────────────────────────────── GUI ─────────────────────────────────────
class SCADAApp:
    """Tkinter dashboard. Compliance limits: GB Grid Code and EREC G5/5 (below)."""

    # GB Grid Code (manuscript reference [26]; limits checked against Rev. 37, 13 April 2026)
    F_NOM_MIN = 49.5
    F_NOM_MAX = 50.5
    F_STAT_MIN = 47.0
    F_STAT_MAX = 52.0
    # EREC G5/5 planning limit
    THD_LIMIT = 5.0
    # GB Grid Code PF floor (Change 16 applied below)
    PF_FLOOR = 0.95
    # V_rms nominal at ADC node (value set from an earlier draft's design section)
    V_NOM = 0.928
    # 90 % undervoltage threshold
    V_UV = 0.835

    # Colour palette (GitHub-dark style)
    C_BG = '#0d1117'
    C_PANEL = '#161b22'
    C_BORDER = '#30363d'
    C_ACCENT = '#00ff88'
    C_HEADER = '#58a6ff'
    C_LABEL = '#8b949e'
    C_OK_BG = '#0d3520'
    C_OK_FG = '#3fb950'
    C_BAD_BG = '#5c0d0d'
    C_BAD_FG = '#f85149'
    C_COMP = '#d2a8ff'

    def __init__(self, root, data, f2_reader=None, out_dir="."):
        self.root = root
        self.data = data
        self.f2_reader = f2_reader  # None if Feather 2 not in use
        self.out_dir = out_dir
        self.wave_phase = 0.0
        self._last_logged_count = 0

        root.title("MIMO RTU - SCADA HMI Dashboard")
        root.geometry("1100x870")
        root.configure(bg=self.C_BG)
        root.resizable(True, True)

        self._build_styles()
        self._build_title_bar()
        self._build_main_layout()

        # 2 Hz refresh cadence
        root.after(500, self.update_ui)

    # ----- styling -----
    def _build_styles(self):
        st = ttk.Style()
        st.theme_use('clam')
        st.configure('Title.TLabel',
                     background=self.C_PANEL, foreground=self.C_ACCENT,
                     font=('Consolas', 16, 'bold'))
        st.configure('Header.TLabel',
                     background=self.C_BG, foreground=self.C_HEADER,
                     font=('Consolas', 10, 'bold'))
        st.configure('Value.TLabel',
                     background=self.C_BG, foreground=self.C_ACCENT,
                     font=('Consolas', 10))
        st.configure('Label.TLabel',
                     background=self.C_BG, foreground=self.C_LABEL,
                     font=('Consolas', 9))
        st.configure('AlarmOK.TLabel',
                     background=self.C_OK_BG, foreground=self.C_OK_FG,
                     font=('Consolas', 11, 'bold'), padding=6)
        st.configure('AlarmBAD.TLabel',
                     background=self.C_BAD_BG, foreground=self.C_BAD_FG,
                     font=('Consolas', 11, 'bold'), padding=6)
        st.configure('Status.TLabel',
                     background=self.C_PANEL, foreground=self.C_LABEL,
                     font=('Consolas', 9))
        st.configure('Frame.TFrame', background=self.C_BG)

    # ----- title bar -----
    def _build_title_bar(self):
        bar = tk.Frame(self.root, bg=self.C_PANEL, height=46)
        bar.pack(fill='x')
        bar.pack_propagate(False)
        ttk.Label(bar, text="  MIMO RTU - SCADA HMI Dashboard",
                  style='Title.TLabel').pack(side='left', padx=10, pady=8)
        self.conn_label = ttk.Label(bar, text="DISCONNECTED",
                                    style='Title.TLabel',
                                    foreground=self.C_BAD_FG)
        self.conn_label.pack(side='right', padx=20, pady=8)

    # ----- main layout -----
    def _build_main_layout(self):
        main = ttk.Frame(self.root, style='Frame.TFrame')
        main.pack(fill='both', expand=True, padx=8, pady=4)

        left = ttk.Frame(main, style='Frame.TFrame')
        left.pack(side='left', fill='y', padx=4)
        right = ttk.Frame(main, style='Frame.TFrame')
        right.pack(side='right', fill='both', expand=True, padx=4)

        # ----- left column: numeric readouts -----
        self._section(left, "VOLTAGES (RMS)")
        self.lbl_Va = self._row(left, "Va_rms", "V")
        self.lbl_Vb = self._row(left, "Vb_rms", "V")
        self.lbl_Vc = self._row(left, "Vc_rms", "V")

        self._section(left, "CURRENTS (RMS)")
        self.lbl_Ia = self._row(left, "Ia_rms", "A")
        self.lbl_Ib = self._row(left, "Ib_rms", "A")
        self.lbl_Ic = self._row(left, "Ic_rms", "A")
        # Ia is not sampled by Feather 1 (M1 covers 5 channels: Va, Vb, Vc,
        # Ib, Ic - see the channel map, Figure 8 of the paper). Show a persistent
        # "not sampled" placeholder in the muted label colour so the
        # three-phase panel is visually complete without misrepresenting
        # what is measured.
        self.lbl_Ia.configure(text="\u2014  (not sampled)",
                              foreground=self.C_LABEL)

        self._section(left, "FREQUENCY")
        self.lbl_fa = self._row(left, "f_Va", "Hz")
        self.lbl_fb = self._row(left, "f_Vb", "Hz")
        self.lbl_fc = self._row(left, "f_Vc", "Hz")

        self._section(left, "THD")
        self.lbl_THDa = self._row(left, "THD_Va", "%")
        self.lbl_THDb = self._row(left, "THD_Vb", "%")
        self.lbl_THDc = self._row(left, "THD_Vc", "%")

        self._section(left, "PHASE & POWER")
        self.lbl_Phb = self._row(left, "Ph_Vb/Ib", "deg")
        self.lbl_Pb = self._row(left, "Pb", "W")
        self.lbl_PFb = self._row(left, "PFb", "")
        self.lbl_Phc = self._row(left, "Ph_Vc/Ic", "deg")
        self.lbl_Pc = self._row(left, "Pc", "W")
        self.lbl_PFc = self._row(left, "PFc", "")

        # ----- right column: alarms, waveform, compliance, hex -----
        self._section(right, "ALARMS")
        tile_row = ttk.Frame(right, style='Frame.TFrame')
        tile_row.pack(fill='x', padx=4)
        self.tile_freq = self._tile(tile_row, "FREQ")
        self.tile_uv = self._tile(tile_row, "U/V ")
        self.tile_thd = self._tile(tile_row, "THD ")
        self.tile_pf = self._tile(tile_row, "PF  ")

        self._section(right, "NESO REPLAY (capture mode)")
        rep_row = tk.Frame(right, bg=self.C_BG)
        rep_row.pack(fill='x', padx=4, pady=(2, 0))
        self.btn_replay_start = tk.Button(
            rep_row, text="Start Replay", command=self._on_replay_start,
            bg='#238636', fg='#ffffff', font=('Consolas', 10, 'bold'),
            borderwidth=0, padx=18, pady=4, activebackground='#2ea043',
            activeforeground='#ffffff', state='disabled')
        self.btn_replay_start.pack(side='left', padx=(0, 6))
        self.btn_replay_abort = tk.Button(
            rep_row, text="Abort", command=self._on_replay_abort,
            bg='#5c0d0d', fg='#ffffff', font=('Consolas', 10),
            borderwidth=0, padx=14, pady=4, activebackground='#7c1818',
            activeforeground='#ffffff', state='disabled')
        self.btn_replay_abort.pack(side='left', padx=(0, 12))
        self.replay_link_label = ttk.Label(
            rep_row, text="Feather 2: not connected", style='Status.TLabel')
        self.replay_link_label.pack(side='left', padx=4)
        self.replay_status_label = tk.Label(
            right, text="Status: idle  (connect Feather 2 to enable replay)",
            bg=self.C_BG, fg=self.C_LABEL, font=('Consolas', 9),
            anchor='w', justify='left')
        self.replay_status_label.pack(fill='x', padx=4, pady=(2, 4))

        self._section(right, "THREE-PHASE WAVEFORM")
        self.canvas = tk.Canvas(right, width=520, height=170, bg=self.C_BG,
                                highlightthickness=1,
                                highlightbackground=self.C_BORDER)
        self.canvas.pack(padx=4, pady=4, fill='x')

        self._section(right, "PDS COMPLIANCE")
        self.compliance_text = tk.Text(right, height=4, width=64,
                                       bg=self.C_BG, fg=self.C_COMP,
                                       font=('Consolas', 9),
                                       borderwidth=1, relief='solid')
        self.compliance_text.pack(padx=4, pady=4, fill='x')

        self._section(right, "DNP3 FRAME (RAW HEX)")
        self.dnp3_text = tk.Text(right, height=6, width=64,
                                 bg=self.C_BG, fg=self.C_ACCENT,
                                 font=('Consolas', 9),
                                 borderwidth=1, relief='solid')
        self.dnp3_text.pack(padx=4, pady=4, fill='both', expand=True)

        # ----- bottom status bar -----
        bar = tk.Frame(self.root, bg=self.C_PANEL, height=24)
        bar.pack(fill='x', side='bottom')
        self.status_label = ttk.Label(bar, text="Frames: 0  |  Last: Never",
                                      style='Status.TLabel')
        self.status_label.pack(side='left', padx=10, pady=3)

    # ----- helpers -----
    def _section(self, parent, title):
        ttk.Label(parent, text=title, style='Header.TLabel').pack(
            anchor='w', padx=4, pady=(10, 2))

    def _row(self, parent, key, unit):
        f = ttk.Frame(parent, style='Frame.TFrame')
        f.pack(fill='x', padx=4, pady=1)
        ttk.Label(f, text=f"{key}:", style='Label.TLabel',
                  width=12).pack(side='left')
        val = ttk.Label(f, text=f"0.00 {unit}".strip(), style='Value.TLabel')
        val.pack(side='left')
        return val

    def _tile(self, parent, name):
        t = ttk.Label(parent, text=f"{name}  OK", style='AlarmOK.TLabel',
                      width=12)
        t.pack(side='left', padx=4, pady=4)
        return t

    # ----- 2 Hz refresh -----
    def update_ui(self):
        d = self.data

        # connection indicator
        if d.connected:
            self.conn_label.configure(text="CONNECTED",
                                      foreground=self.C_OK_FG)
        else:
            self.conn_label.configure(text="DISCONNECTED",
                                      foreground=self.C_BAD_FG)

        # numeric readouts
        self.lbl_Va.configure(text=f"{d.Va_rms:.3f} V")
        self.lbl_Vb.configure(text=f"{d.Vb_rms:.3f} V")
        self.lbl_Vc.configure(text=f"{d.Vc_rms:.3f} V")
        self.lbl_Ib.configure(text=f"{d.Ib_rms:.4f} A")
        self.lbl_Ic.configure(text=f"{d.Ic_rms:.4f} A")
        self.lbl_fa.configure(text=f"{d.f_Va:.2f} Hz")
        self.lbl_fb.configure(text=f"{d.f_Vb:.2f} Hz")
        self.lbl_fc.configure(text=f"{d.f_Vc:.2f} Hz")
        self.lbl_THDa.configure(text=f"{d.THD_Va:.2f} %")
        self.lbl_THDb.configure(text=f"{d.THD_Vb:.2f} %")
        self.lbl_THDc.configure(text=f"{d.THD_Vc:.2f} %")
        self.lbl_Phb.configure(text=f"{d.Ph_VbIb:.2f} deg")
        self.lbl_Phc.configure(text=f"{d.Ph_VcIc:.2f} deg")
        self.lbl_Pb.configure(text=f"{d.Pb:.4f} W")
        self.lbl_Pc.configure(text=f"{d.Pc:.4f} W")
        self.lbl_PFb.configure(text=f"{d.PFb:.3f}")
        self.lbl_PFc.configure(text=f"{d.PFc:.3f}")

        # firmware-flagged alarms
        self._set_tile(self.tile_freq, "FREQ", d.freq_alarm)
        self._set_tile(self.tile_uv,   "U/V ", d.uv_alarm)
        self._set_tile(self.tile_thd,  "THD ", d.thd_alarm)

        # PF check - Change 16 correction: floor 0.95, label "FAIL <0.95"
        active_pf = [p for p in (d.PFb, d.PFc) if p > 0]
        pf_ok = all(p >= self.PF_FLOOR for p in active_pf) if active_pf else True
        self._set_pf_tile(self.tile_pf, pf_ok)

        # status bar
        self.status_label.configure(
            text=f"Frames: {d.frame_count}  |  Bytes: {d.dnp3_bytes}  "
                 f"|  Last: {d.last_update}")

        # waveform + compliance + DNP3 hex view
        self._redraw_waveform()
        self._update_compliance(pf_ok)
        self._update_dnp3_view()
        self._update_replay_panel()

        self.root.after(500, self.update_ui)

    # ----- replay control -----
    def _on_replay_start(self):
        """Open a CSV logger and tell Feather 2 to begin the trajectory."""
        if self.f2_reader is None or not self.data.f2_connected:
            messagebox.showwarning(
                "Feather 2 not connected",
                "Cannot start replay - Feather 2 serial link is down.")
            return
        if self.data.replay_active:
            return  # already running
        try:
            logger = ReplayLogger.start(self.out_dir)
        except Exception as e:
            messagebox.showerror("CSV error",
                                 f"Could not open replay CSV:\n{e}")
            return
        # Reset latch and arm logger BEFORE issuing the command so we don't
        # miss the very first REPLAY token.
        self.data.replay_complete = False
        self.data.replay_index = -1
        self.data.replay_target_f = float('nan')
        self.data.replay_logger = logger
        self._last_logged_count = self.data.frame_count

        ok = self.f2_reader.send_command("REPLAY_START")
        if not ok:
            self.data.replay_logger = None
            try:
                logger.close()
                os.remove(logger.path)
            except Exception:
                pass
            messagebox.showerror(
                "Send failed",
                "Could not write REPLAY_START to Feather 2.")
            return
        print("[SCADAApp] replay armed")

    def _on_replay_abort(self):
        """Send REPLAY_ABORT to Feather 2 and close the logger."""
        if self.f2_reader is not None:
            self.f2_reader.send_command("REPLAY_ABORT")
        # Force-close the logger; the Feather 2 thread will set
        # replay_complete=True when it sees REPLAY_END, but on a manual
        # abort we don't wait for that.
        self.data.replay_active = False
        if self.data.replay_logger is not None:
            try:
                self.data.replay_logger.close()
            except Exception:
                pass
            self.data.replay_logger = None
        self.replay_status_label.configure(
            text="Status: aborted by operator  (CSV closed)",
            fg=self.C_BAD_FG)
        print("[SCADAApp] replay aborted")

    def _finalise_replay_run(self):
        """Called once per run, when REPLAY_END is detected. Closes the CSV
        and generates the sanity-check overlay PNG."""
        logger = self.data.replay_logger
        if logger is None:
            return
        path = logger.path
        try:
            logger.close()
        except Exception as e:
            print(f"[_finalise] close failed: {e}")
        self.data.replay_logger = None
        self.replay_status_label.configure(
            text=f"Status: complete  ({os.path.basename(path)})  "
                 f"- rendering overlay...",
            fg=self.C_OK_FG)
        self.root.update_idletasks()

        # Generate the overlay - runs synchronously in the UI thread because
        # matplotlib is not safe to call from worker threads on some platforms.
        # The replay is finished, so blocking the UI for 1-2 seconds is fine.
        neso_ref = os.path.join(self.out_dir, "neso_trajectory.csv")
        png = generate_overlay_plot(path, neso_csv_path=neso_ref)
        if png:
            self.replay_status_label.configure(
                text=f"Status: complete  CSV={os.path.basename(path)}  "
                     f"PNG={os.path.basename(png)}",
                fg=self.C_OK_FG)
        else:
            self.replay_status_label.configure(
                text=f"Status: complete  CSV={os.path.basename(path)}  "
                     f"(overlay skipped - see console)",
                fg=self.C_OK_FG)

    def _update_replay_panel(self):
        """Refresh button-enable states and the status line each tick."""
        d = self.data
        f2_up = (self.f2_reader is not None) and d.f2_connected

        # Link label
        if self.f2_reader is None:
            self.replay_link_label.configure(text="Feather 2: not configured")
        elif f2_up:
            self.replay_link_label.configure(text="Feather 2: connected",
                                             foreground=self.C_OK_FG)
        else:
            self.replay_link_label.configure(text="Feather 2: link down",
                                             foreground=self.C_BAD_FG)

        # Button enable/disable
        if d.replay_active:
            self.btn_replay_start.configure(state='disabled')
            self.btn_replay_abort.configure(state='normal')
        elif f2_up:
            self.btn_replay_start.configure(state='normal')
            self.btn_replay_abort.configure(state='disabled')
        else:
            self.btn_replay_start.configure(state='disabled')
            self.btn_replay_abort.configure(state='disabled')

        # Status line during a run
        if d.replay_active and d.replay_index >= 0:
            err_mhz = ((d.f_Vb - d.replay_target_f) * 1000.0
                       if (d.f_Vb > 0 and not math.isnan(d.replay_target_f))
                       else float('nan'))
            err_str = (f"{err_mhz:+.0f} mHz" if not math.isnan(err_mhz)
                       else "n/a")
            self.replay_status_label.configure(
                text=(f"Status: RUNNING  "
                      f"{d.replay_index + 1}/{d.replay_total}   "
                      f"target {d.replay_target_f:.4f} Hz   "
                      f"f_Vb {d.f_Vb:.4f} Hz   "
                      f"Δ {err_str}"),
                fg=self.C_ACCENT)

        # Detect REPLAY_END (latched by Feather 2 reader)
        if d.replay_complete and d.replay_logger is not None:
            d.replay_complete = False  # consume the latch
            self._finalise_replay_run()
        elif (not d.replay_active and d.replay_logger is None
              and not f2_up and self.f2_reader is None):
            # Initial "no feather 2" hint - leave the existing text alone
            pass

    def _set_tile(self, tile, name, is_alarm):
        if is_alarm:
            tile.configure(text=f"{name}  ALARM", style='AlarmBAD.TLabel')
        else:
            tile.configure(text=f"{name}  OK   ", style='AlarmOK.TLabel')

    def _set_pf_tile(self, tile, pf_ok):
        if pf_ok:
            tile.configure(text="PF    OK   ", style='AlarmOK.TLabel')
        else:
            tile.configure(text="PF    FAIL <0.95", style='AlarmBAD.TLabel')

    def _redraw_waveform(self):
        d = self.data
        c = self.canvas
        c.delete('all')
        w = max(int(c.winfo_width()), 100)
        h = max(int(c.winfo_height()), 50)
        mid = h // 2
        # zero line
        c.create_line(0, mid, w, mid, fill=self.C_BORDER, dash=(2, 4))

        # animation phase advance (purely visual)
        self.wave_phase = (self.wave_phase + 0.10) % (2 * math.pi)

        Vmax = max(d.Va_rms, d.Vb_rms, d.Vc_rms, 1e-6)
        amp = (h / 2) - 8

        def draw_phase(rms, offset_deg, colour):
            if rms <= 0:
                return
            scale = rms / Vmax
            pts = []
            for x in range(0, w, 2):
                t = (x / w) * 4 * math.pi  # 2 cycles across canvas
                y = mid - amp * scale * math.sin(
                    t + math.radians(offset_deg) + self.wave_phase)
                pts.extend([x, y])
            if len(pts) >= 4:
                c.create_line(*pts, fill=colour, width=2, smooth=True)

        draw_phase(d.Va_rms,    0, self.C_BAD_FG)   # red
        draw_phase(d.Vb_rms, -120, self.C_OK_FG)    # green
        draw_phase(d.Vc_rms, -240, self.C_HEADER)   # blue

    def _update_compliance(self, pf_ok):
        d = self.data
        active_f = [f for f in (d.f_Va, d.f_Vb, d.f_Vc) if f > 0]
        f_avg = sum(active_f) / len(active_f) if active_f else 0.0
        thd_max = max(d.THD_Va, d.THD_Vb, d.THD_Vc)

        # checks
        f_grid_ok = (self.F_NOM_MIN <= f_avg <= self.F_NOM_MAX) if f_avg else False
        f_stat_ok = (self.F_STAT_MIN <= f_avg <= self.F_STAT_MAX) if f_avg else False
        thd_ok = (thd_max <= self.THD_LIMIT)

        t = self.compliance_text
        t.delete('1.0', tk.END)
        t.insert(tk.END,
                 f"Grid Code  (49.5-50.5 Hz):  "
                 f"{'PASS' if f_grid_ok else 'FAIL'}   f_avg = {f_avg:.3f} Hz\n")
        t.insert(tk.END,
                 f"Statutory  (47.0-52.0 Hz):  "
                 f"{'WITHIN' if f_stat_ok else 'OUT  '}   (Grid Code [26])\n")
        t.insert(tk.END,
                 f"EREC G5/5  (THD <= 5 %):    "
                 f"{'PASS' if thd_ok else 'FAIL'}   THD_max = {thd_max:.2f} %\n")
        t.insert(tk.END,
                 f"GBGC PF    (>= 0.95):       "
                 f"{'PASS' if pf_ok else 'FAIL <0.95'}   "
                 f"PFb = {d.PFb:.3f}, PFc = {d.PFc:.3f}")

    def _update_dnp3_view(self):
        d = self.data
        self.dnp3_text.delete('1.0', tk.END)
        hex_bytes = d.raw_frame.strip().split()
        if not hex_bytes:
            self.dnp3_text.insert(tk.END,
                                  "(no DNP3 frame received yet)")
            return
        # 16 bytes per line
        for i in range(0, len(hex_bytes), 16):
            self.dnp3_text.insert(tk.END,
                                  ' '.join(hex_bytes[i:i + 16]) + '\n')


# ─────────────────────────── Port Discovery ───────────────────────────────
def _guess_feather_index(ports):
    """Return the index of the most likely Feather/USB-CDC port in `ports`,
    or 0 if no clear match. Used to pre-select the dropdown."""
    keywords = ('cp210', 'ch340', 'ch910', 'silicon labs', 'wch',
                'usb-serial', 'usb serial', 'feather', 'esp32')
    for i, p in enumerate(ports):
        desc = (p.description or '').lower()
        manu = (getattr(p, 'manufacturer', '') or '').lower()
        if any(k in desc or k in manu for k in keywords):
            return i
    return 0


# ─────────────────────────── Port Selector Dialog ─────────────────────────
class PortSelector:
    """Modal startup dialog. Lets the user pick a COM port from a dropdown
    before the main HMI window opens. Pre-highlights the most likely Feather
    port; a Refresh button re-scans after plugging in the board."""

    BG = '#0d1117'
    PANEL = '#161b22'
    ACCENT = '#00ff88'
    HEADER = '#58a6ff'
    LABEL = '#8b949e'
    OK_FG = '#3fb950'
    BAD_FG = '#f85149'
    BTN_BG = '#21262d'
    BTN_FG = '#c9d1d9'
    CONNECT_BG = '#238636'

    def __init__(self):
        self.selected_port = None       # Feather 1 (required)
        self.selected_port_f2 = None    # Feather 2 (optional, may be None)
        self._port_map = {}

        self.root = tk.Tk()
        self.root.title("MIMO RTU - Select Serial Port(s)")
        self.root.geometry("560x340")
        self.root.configure(bg=self.BG)
        self.root.resizable(False, False)

        # Re-use the same dark style language as the main HMI
        st = ttk.Style(self.root)
        st.theme_use('clam')
        st.configure('SelTitle.TLabel', background=self.BG,
                     foreground=self.ACCENT, font=('Consolas', 13, 'bold'))
        st.configure('SelLabel.TLabel', background=self.BG,
                     foreground=self.HEADER, font=('Consolas', 10))
        st.configure('SelStatus.TLabel', background=self.BG,
                     foreground=self.LABEL, font=('Consolas', 9))
        st.configure('Sel.TCombobox', fieldbackground=self.PANEL,
                     foreground=self.ACCENT)

        ttk.Label(self.root, text="MIMO RTU - SCADA HMI",
                  style='SelTitle.TLabel').pack(pady=(14, 4))

        # ----- Feather 1 (required) -----
        ttk.Label(self.root,
                  text="Feather 1 (measurement board) - required:",
                  style='SelLabel.TLabel').pack(pady=(8, 2))
        self.combo = ttk.Combobox(self.root, font=('Consolas', 10),
                                   width=60, state='readonly',
                                   style='Sel.TCombobox')
        self.combo.pack(pady=2, padx=20)

        # ----- Feather 2 (optional) -----
        ttk.Label(self.root,
                  text="Feather 2 (DAC generator) - optional, "
                       "needed for NESO replay:",
                  style='SelLabel.TLabel').pack(pady=(10, 2))
        self.combo_f2 = ttk.Combobox(self.root, font=('Consolas', 10),
                                      width=60, state='readonly',
                                      style='Sel.TCombobox')
        self.combo_f2.pack(pady=2, padx=20)

        self.status = ttk.Label(self.root, text="", style='SelStatus.TLabel')
        self.status.pack(pady=(8, 4))

        btns = tk.Frame(self.root, bg=self.BG)
        btns.pack(pady=12)
        tk.Button(btns, text="Refresh", command=self._refresh,
                  bg=self.BTN_BG, fg=self.BTN_FG, font=('Consolas', 10),
                  borderwidth=0, padx=14, pady=4,
                  activebackground='#30363d',
                  activeforeground=self.BTN_FG).pack(side='left', padx=6)
        tk.Button(btns, text="Connect", command=self._connect,
                  bg=self.CONNECT_BG, fg='#ffffff',
                  font=('Consolas', 10, 'bold'),
                  borderwidth=0, padx=22, pady=4,
                  activebackground='#2ea043',
                  activeforeground='#ffffff').pack(side='left', padx=6)
        tk.Button(btns, text="Cancel", command=self._cancel,
                  bg=self.BTN_BG, fg=self.BTN_FG, font=('Consolas', 10),
                  borderwidth=0, padx=14, pady=4,
                  activebackground='#30363d',
                  activeforeground=self.BTN_FG).pack(side='left', padx=6)

        self.root.protocol("WM_DELETE_WINDOW", self._cancel)
        # Enter key triggers Connect
        self.root.bind('<Return>', lambda e: self._connect())

        self._refresh()

    def _refresh(self):
        ports = list(serial.tools.list_ports.comports())
        if not ports:
            self.combo['values'] = []
            self.combo.set('')
            self.combo_f2['values'] = ['(none - skip Feather 2)']
            self.combo_f2.set('(none - skip Feather 2)')
            self._port_map = {}
            self.status.configure(
                text="No serial ports found. Plug in the Feather(s) and click Refresh.",
                foreground=self.BAD_FG)
            return
        labels = []
        self._port_map = {}
        for p in ports:
            label = f"{p.device}   -   {(p.description or 'Unknown device').strip()}"
            labels.append(label)
            self._port_map[label] = p.device
        self.combo['values'] = labels
        # Feather 2 dropdown includes "(none)" so the user can run without it
        labels_f2 = ['(none - skip Feather 2)'] + labels
        self.combo_f2['values'] = labels_f2

        primary = _guess_feather_index(ports)
        self.combo.current(primary)
        # Pre-select the next port (different from primary) for Feather 2 if
        # there is one - the two Feathers usually enumerate adjacently.
        if len(ports) >= 2:
            secondary = (primary + 1) % len(ports)
            self.combo_f2.current(secondary + 1)  # +1 for the "(none)" entry
        else:
            self.combo_f2.current(0)  # (none)

        self.status.configure(
            text=f"Found {len(ports)} port(s). Pre-selected the most likely "
                 f"pair - change if needed, then Connect.",
            foreground=self.OK_FG)

    def _connect(self):
        sel = self.combo.get()
        if not sel:
            self.status.configure(text="No Feather 1 port selected.",
                                  foreground=self.BAD_FG)
            return
        self.selected_port = self._port_map.get(sel, sel.split()[0])
        # Feather 2 is optional
        sel_f2 = self.combo_f2.get()
        if sel_f2 and sel_f2 != '(none - skip Feather 2)':
            f2_dev = self._port_map.get(sel_f2, sel_f2.split()[0])
            if f2_dev == self.selected_port:
                self.status.configure(
                    text="Feather 1 and Feather 2 cannot share the same port.",
                    foreground=self.BAD_FG)
                return
            self.selected_port_f2 = f2_dev
        else:
            self.selected_port_f2 = None
        self.root.destroy()

    def _cancel(self):
        self.selected_port = None
        self.selected_port_f2 = None
        self.root.destroy()

    def run(self):
        self.root.mainloop()
        return self.selected_port, self.selected_port_f2


# ──────────────────────────────── Main ────────────────────────────────────
def main():
    # Command-line override still works for scripted use:
    #   python SCADA_HMI.py COM5            (Feather 1 only)
    #   python SCADA_HMI.py COM5 COM6       (Feather 1 + Feather 2)
    f2_port = None
    if len(sys.argv) > 1:
        port = sys.argv[1]
        if len(sys.argv) > 2:
            f2_port = sys.argv[2]
        print(f"Using ports from command line: F1={port}, F2={f2_port}")
    else:
        port, f2_port = PortSelector().run()
        if port is None:
            print("Cancelled - no port selected.")
            sys.exit(0)
        print(f"Selected ports: F1={port}, F2={f2_port}")

    data = RTUData()
    reader = SerialReader(port, 115200, data)
    reader.start()

    f2_reader = None
    if f2_port:
        f2_reader = Feather2Reader(f2_port, 115200, data)
        f2_reader.start()

    root = tk.Tk()
    out_dir = os.path.dirname(os.path.abspath(__file__))
    app = SCADAApp(root, data, f2_reader=f2_reader, out_dir=out_dir)

    def on_close():
        reader.running = False
        if f2_reader is not None:
            f2_reader.running = False
        # Make sure any open CSV is flushed if the user closes mid-run
        if data.replay_logger is not None:
            try:
                data.replay_logger.close()
            except Exception:
                pass
        root.destroy()

    root.protocol("WM_DELETE_WINDOW", on_close)
    root.mainloop()


if __name__ == '__main__':
    main()
