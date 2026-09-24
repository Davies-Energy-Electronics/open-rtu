#!/usr/bin/env python3
"""
log_capture.py - reliable serial logger for characterise.ino captures.

The Arduino IDE's serial monitor drops lines at 921600 baud, which silently
corrupts a capture: you get a file that looks fine and a standard deviation
that means nothing. This reads the port directly and checks what it got.

Usage:
    python log_capture.py COM6 A1_quiet_single.txt
    python log_capture.py COM6 runB_1.txt --skip-warmup
    python log_capture.py --list

    # with the mode guard added 24 Sep 2026 - always use it:
    python log_capture.py COM6 A2_quiet_roundrobin_run1.txt ^
        --expect-mode roundrobin --expect-label A2_quiet_roundrobin_run1

Opening the port resets the ESP32, so the logger always sees the run from
boot. Warm-up is 600 s; --skip-warmup sends the keypress that shortcuts it,
which is what you want for every capture after the first.

Needs pyserial:  pip install pyserial
"""
from __future__ import annotations
import argparse, os, sys, time

try:
    import serial
    from serial.tools import list_ports
except ImportError:
    sys.exit("pyserial is not installed.  Run:   pip install pyserial")

BAUD     = 921600
EXPECT_N = 16384


def _is_sample_row(line):
    """True for a data row, in either format the sketches emit."""
    t = line.strip()
    if not t or t.startswith("#"):
        return False
    if t.lstrip("-").isdigit():
        return True
    parts = t.split(",")
    if len(parts) == 2 and all(p.strip().lstrip("-").isdigit() for p in parts):
        return True
    return False


def list_serial_ports():
    ports = list(list_ports.comports())
    if not ports:
        print("No serial ports found. Is the board plugged in?")
        return
    print("\n  Available serial ports:")
    for p in ports:
        print(f"    {p.device:<10} {p.description}")
    print()


def open_port(port, baud, reset=True):
    """Open without letting pyserial assert DTR/RTS.

    pyserial raises both lines on open. On the ESP32 auto-reset circuit RTS
    drives EN and DTR drives GPIO0, so the default open holds the board in
    reset or drops it into the download bootloader - and nothing is ever
    transmitted. Both lines must be deasserted BEFORE the port is opened.

    With that done, the esptool run-mode reset is: GPIO0 high (DTR False)
    throughout, EN pulsed low (RTS True) and released (RTS False).
    """
    ser = serial.Serial()
    ser.port = port
    ser.baudrate = baud
    ser.timeout = 2
    ser.dtr = False
    ser.rts = False
    ser.open()

    # Some drivers ignore the pre-open settings; assert them again.
    try:
        ser.dtr = False
        ser.rts = False
    except Exception:
        pass

    if reset:
        try:
            ser.setDTR(False)      # GPIO0 high -> normal boot, not download
            ser.setRTS(True)       # EN low  -> held in reset
            time.sleep(0.12)
            ser.setRTS(False)      # EN high -> run
        except Exception:
            pass
    return ser


def probe(port, baud, seconds):
    """Dump whatever raw bytes arrive. Diagnostic only."""
    ser = open_port(port, baud)
    print(f"\n  Probing {port} at {baud} baud for {seconds:.0f} s - raw bytes.\n")
    t0 = time.time()
    total = 0
    try:
        while time.time() - t0 < seconds:
            n = ser.in_waiting
            if n:
                data = ser.read(n)
                total += len(data)
                sys.stdout.write(data.decode("utf-8", errors="replace"))
                sys.stdout.flush()
            else:
                time.sleep(0.05)
    except KeyboardInterrupt:
        pass
    finally:
        ser.close()
    print(f"\n\n  {total} bytes received.")
    if total == 0:
        print("  Nothing at all. Try the other baud rate (--baud 115200), check")
        print("  the IDE serial monitor is closed, and confirm the port.\n")
    else:
        print("  Data is flowing. If it looks like garbage the baud is wrong.\n")
    return 0 if total else 1


def capture(port, outfile, skip_warmup, baud, timeout_s, expect=None,
            expect_mode=None, expect_label=None, overwrite=False):
    expect = expect or EXPECT_N

    # Added 24 Sep 2026. A2_quiet_roundrobin_run1.txt was written twice in one
    # session because the same command was re-issued; the mode and label both
    # matched, so every guard passed and the first capture was silently lost.
    # An existing file is now refused before the board is even reset - there is
    # no point spending 600 s on a capture that will overwrite one you kept.
    if os.path.exists(outfile) and not overwrite:
        sys.exit(f"\n  {outfile} already exists.\n"
                 f"  Refusing to overwrite it. Either you have already taken\n"
                 f"  this run, or the filename has not been advanced. Rename\n"
                 f"  the old file, or pass --overwrite if you mean it.\n")

    print(f"\n  Opening {port} at {baud} baud...")
    try:
        ser = open_port(port, baud)
    except serial.SerialException as e:
        sys.exit(f"\n  Could not open {port}: {e}\n"
                 f"  Close the Arduino IDE serial monitor and try again.\n")

    time.sleep(0.3)
    ser.reset_input_buffer()
    print("  Board reset (DTR/RTS deasserted). Waiting for output.")
    print("  You should see '# waiting 600 s...' within about 4 seconds.\n")

    lines, numeric, skipped = [], 0, False
    header = {}
    t_last = time.time()
    t0 = time.time()

    try:
        while True:
            raw = ser.readline()
            if not raw:
                if time.time() - t_last > timeout_s:
                    print(f"\n  No data for {timeout_s:.0f} s - stopping.")
                    break
                continue
            t_last = time.time()

            try:
                line = raw.decode("utf-8", errors="replace").strip()
            except Exception:
                continue
            if not line:
                continue

            lines.append(line)

            if line.strip() == "t_us,code":
                print(f"  {line}   <- data follows")
                continue

            if line.startswith("#"):
                body = line.lstrip("#").strip()
                if "," in body:
                    k, _, v = body.partition(",")
                    header[k.strip()] = v.strip()
                print(f"  {line}")
                if line == "# end":
                    break
                if skip_warmup and not skipped and "warm" in line.lower():
                    ser.write(b"\n")
                    ser.flush()
                    skipped = True
                    print("  -> sent keypress to skip the warm-up")
                continue

            # A sample row is either a bare integer (characterise.ino,
            # adc_latency.ino) or "t_us,code" (characterise_v2loop.ino).
            # Added 23 Sep 2026: without the second form the counter stayed at
            # zero for every v2loop capture and this tool reported SAMPLE COUNT
            # WRONG on a perfectly good file - which would have retired the one
            # guard that catches dropped lines at 921600 baud.
            if _is_sample_row(line):
                numeric += 1
                if numeric % 2000 == 0:
                    print(f"  ... {numeric} samples")
    except KeyboardInterrupt:
        print("\n  Interrupted.")
    finally:
        ser.close()

    with open(outfile, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")

    elapsed = time.time() - t0
    print(f"\n  {'='*58}")
    print(f"  Wrote {outfile}   ({len(lines)} lines, {elapsed:.0f} s)")
    print(f"  Numeric samples: {numeric}   (expected {expect})")

    ok = True
    if numeric != expect:
        print(f"\n  ** SAMPLE COUNT WRONG - the log dropped data or the run")
        print(f"     was cut short. Do NOT analyse this file. Re-run it.")
        ok = False

    mean = header.get("mean_code")
    if mean is not None:
        try:
            m = float(mean)
            print(f"  Mean code: {m:.1f}", end="")
            if m < 200 or m > 3000:
                print("   ** RAILED - the bias is not reaching the ADC.")
                ok = False
            else:
                print("   (sensible - bias present)")
        except ValueError:
            pass

    if header.get("sd_millivolts"):
        print(f"  On-board sd estimate: {header['sd_millivolts']} mV"
              f"   (the analyser recomputes this properly)")

    # Added 24 Sep 2026. On 13 September three captures were
    # taken without re-editing the sketch between uploads, so A2 and C1 both
    # repeated A1's mode and the two measurements were silently lost. The
    # sketch prints what it actually did; the operator is the only thing that
    # was checking, and the operator forgot. So the tool checks now. Refusing
    # a capture is cheap; discovering a month later that three runs were the
    # same run is not.
    got_mode = header.get("mode")
    if expect_mode is not None:
        if got_mode is None:
            print(f"\n  ** NO MODE IN HEADER - cannot confirm this is a"
                  f" '{expect_mode}' capture.")
            ok = False
        elif got_mode.strip().lower() != expect_mode.strip().lower():
            print(f"\n  ** WRONG MODE. The sketch reports '{got_mode}',"
                  f" you asked for '{expect_mode}'.")
            print(f"     The flash did not take, or the #defines were not"
                  f" edited. This is the")
            print(f"     13 September failure. Re-edit, reflash, re-run.")
            ok = False
        else:
            print(f"  Mode: {got_mode}   (matches --expect-mode)")
    elif got_mode:
        print(f"  Mode: {got_mode}   (unchecked - pass --expect-mode to verify)")

    got_label = header.get("label")
    if expect_label is not None:
        if got_label is None or got_label.strip() != expect_label.strip():
            print(f"\n  ** WRONG LABEL. The sketch reports '{got_label}',"
                  f" you asked for '{expect_label}'.")
            print(f"     RUN_LABEL was not updated before this flash, so this"
                  f" capture would be")
            print(f"     filed under the wrong name. Re-edit, reflash, re-run.")
            ok = False
        else:
            print(f"  Label: {got_label}   (matches --expect-label)")
    elif got_label:
        print(f"  Label: {got_label}   (unchecked - pass --expect-label to verify)")

    print(f"  {'='*58}")
    print("  Capture looks good.\n" if ok else "  RE-RUN THIS CAPTURE.\n")
    return 0 if ok else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("port", nargs="?", help="e.g. COM6")
    ap.add_argument("outfile", nargs="?", help="e.g. A1_quiet_single.txt")
    ap.add_argument("--list", action="store_true", help="list serial ports and exit")
    ap.add_argument("--probe", type=float, metavar="SECONDS", default=None,
                    help="dump raw bytes for N seconds instead of capturing "
                         "(diagnostic, when nothing arrives)")
    ap.add_argument("--skip-warmup", action="store_true",
                    help="shortcut the 600 s thermal settle (use after the first run)")
    ap.add_argument("--baud", type=int, default=BAUD)
    ap.add_argument("--expect", type=int, default=None,
                    help="expected sample count (default 16384; latency runs use 20000)")
    ap.add_argument("--timeout", type=float, default=90.0,
                    help="give up after this many seconds of silence (default 90)")
    ap.add_argument("--overwrite", action="store_true",
                    help="allow an existing output file to be replaced")
    ap.add_argument("--expect-mode", default=None, metavar="MODE",
                    help="fail the capture unless the sketch's '# mode,' header "
                         "matches, e.g. single, roundrobin, sweep_single")
    ap.add_argument("--expect-label", default=None, metavar="LABEL",
                    help="fail the capture unless the sketch's '# label,' header "
                         "matches RUN_LABEL exactly")
    a = ap.parse_args(argv)

    if a.list or not a.port:
        list_serial_ports()
        if not a.port:
            print("  Then:  python log_capture.py COM6 A1_quiet_single.txt\n")
        return 0
    if a.probe:
        return probe(a.port, a.baud, a.probe)

    if not a.outfile:
        sys.exit("  Give an output filename, e.g. A1_quiet_single.txt")

    return capture(a.port, a.outfile, a.skip_warmup, a.baud, a.timeout, a.expect,
                   a.expect_mode, a.expect_label, a.overwrite)


if __name__ == "__main__":
    raise SystemExit(main())
