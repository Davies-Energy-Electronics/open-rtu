#!/usr/bin/env python3
"""
capture_v2.py - capture the unmodified V2 build's serial output.

WHY THIS EXISTS AND log_capture.py DOES NOT DO IT
log_capture.py is built for the characterisation sketches: it echoes lines
beginning with '#', counts bare-integer or "t_us,code" sample rows, and stops
on "# end". The released V2 firmware emits none of those. Its output is

    [DNP3 FRAME | 64 bytes | Src:0x0001 Dst:0x0003]
    05 64 3A 44 ...
    [G30V5 ANALOGUE]
      f_Va:49.9987Hz  f_Vb:50.0013Hz  f_Vc:50.0002Hz

so log_capture.py collects it correctly into the file but prints NOTHING to
the screen and never terminates, which looks exactly like a dead port. This
script shows progress, counts what matters, and stops on its own.

It reuses log_capture.py's port-opening logic, which matters: pyserial
asserts DTR and RTS on open, and on the ESP32 auto-reset circuit that holds
the board in reset or drops it into the bootloader. Both lines must be
deasserted before the port is opened.

Usage:
    python capture_v2.py COM6 runE_v2.txt
    python capture_v2.py COM6 runE_v2.txt --seconds 300
    python capture_v2.py COM6 runE_v2.txt --min-windows 500

Needs pyserial:  pip install pyserial
"""
from __future__ import annotations
import argparse, re, sys, time

try:
    import serial
except ImportError:
    sys.exit("pyserial is not installed.  Run:   pip install pyserial")

# Reuse the tested open_port() so the DTR/RTS handling is identical.
try:
    from log_capture import open_port
except Exception:
    open_port = None

BAUD = 115200
FVB = re.compile(r"f_Vb\s*[:=]\s*([0-9]+\.?[0-9]*)", re.I)


def _open(port, baud):
    if open_port is not None:
        return open_port(port, baud)
    ser = serial.Serial()
    ser.port = port; ser.baudrate = baud; ser.timeout = 2
    ser.dtr = False; ser.rts = False
    ser.open()
    try:
        ser.dtr = False; ser.rts = False
        ser.setDTR(False); ser.setRTS(True); time.sleep(0.12); ser.setRTS(False)
    except Exception:
        pass
    return ser


def capture(port, outfile, seconds, min_windows, baud):
    print(f"\n  Opening {port} at {baud} baud...")
    try:
        ser = _open(port, baud)
    except serial.SerialException as e:
        sys.exit(f"\n  Could not open {port}: {e}\n"
                 f"  Close the Arduino IDE serial monitor and try again.\n")
    time.sleep(0.3)
    ser.reset_input_buffer()
    print("  Board reset. Expect the V2 banner within about 5 seconds.")
    print(f"  Running for up to {seconds:.0f} s, or {min_windows} f_Vb windows.\n")

    lines, windows, banner = [], 0, False
    t0 = time.time()
    t_last = t0
    try:
        while True:
            raw = ser.readline()
            now = time.time()
            if raw:
                t_last = now
                line = raw.decode("utf-8", errors="replace").rstrip("\r\n")
                if line.strip():
                    lines.append(line)
                    if not banner and "DNP3" in line.upper():
                        banner = True
                        print("  V2 banner seen - the right firmware is running.")
                    if FVB.search(line):
                        windows += 1
                        if windows % 50 == 0:
                            print(f"  ... {windows} windows  ({now - t0:.0f} s)")
            if now - t0 >= seconds:
                print(f"\n  Reached {seconds:.0f} s.")
                break
            if min_windows and windows >= min_windows:
                print(f"\n  Reached {min_windows} windows.")
                break
            if not raw and now - t_last > 20:
                print("\n  No data for 20 s - stopping.")
                break
    except KeyboardInterrupt:
        print("\n  Interrupted by user.")
    finally:
        ser.close()

    with open(outfile, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")

    print("\n  " + "=" * 58)
    print(f"  Wrote {outfile}   ({len(lines)} lines, {time.time() - t0:.0f} s)")
    print(f"  f_Vb windows captured: {windows}")
    ok = True
    if not banner:
        print("\n  ** NO V2 BANNER SEEN.")
        print("     Either the wrong sketch is on the board, or the baud is wrong.")
        print("     V2 prints at 115200. The characterisation sketch prints at")
        print("     921600 - listening at the wrong rate gives silence or junk.")
        ok = False
    if windows < 100:
        print(f"\n  ** ONLY {windows} WINDOWS. Want at least a few hundred for a")
        print("     dispersion figure. Re-run for longer.")
        ok = False
    print("  " + "=" * 58)
    print("  Capture looks good.\n" if ok else "  RE-RUN THIS CAPTURE.\n")
    return 0 if ok else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("port")
    ap.add_argument("outfile")
    ap.add_argument("--seconds", type=float, default=300.0)
    ap.add_argument("--min-windows", type=int, default=600)
    ap.add_argument("--baud", type=int, default=BAUD)
    a = ap.parse_args(argv)
    return capture(a.port, a.outfile, a.seconds, a.min_windows, a.baud)


if __name__ == "__main__":
    raise SystemExit(main())
