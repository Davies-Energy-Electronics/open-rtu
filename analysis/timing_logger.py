import serial, time
from datetime import datetime

f1 = serial.Serial()
f1.port, f1.baudrate, f1.timeout = "COM6", 115200, 2
f1.dtr = False; f1.rts = False        # release control lines - do not reset/hold the board
f1.open()
fname = f"runtime_capture_{datetime.now().strftime('%Y%m%d')}.csv"
print("COM6 open (DTR/RTS released).")
print(">>> NOW PRESS THE RESET BUTTON ON FEATHER 1 <<<")
ntim = 0; t0 = time.time()
with open(fname, "w", encoding="utf-8", newline="") as fh:
    while (time.time() - t0) < 1900:
        line = f1.readline().decode(errors="replace")
        if not line: continue
        fh.write(line); fh.flush()
        if ntim == 0 and not line.startswith("TIMING"):
            print(f"  [boot] {line.strip()[:90]}")
        if line.startswith("TIMING"):
            ntim += 1
            if ntim in (1, 100) or ntim % 100 == 0:
                print(f"  {ntim} TIMING iterations, t={(time.time()-t0)/60:4.1f} min")
        if "TIMING_CAPTURE_COMPLETE" in line:
            print("  Sentinel seen - capture complete."); break
f1.close()
print(f"Done: {ntim} TIMING lines -> {fname}")
