import serial, time
from datetime import datetime

F1, DURATION = "COM6", 400          # seconds, matches the canonical 400 s capture shape
f1 = serial.Serial(F1, 115200, timeout=1)
print("COM6 open - 5 s settle after DTR reset...")
time.sleep(5)

for run in range(1, 4):
    input(f"\n=== RUN {run}/3 - press Enter to START LOGGING, then send REPLAY_START in the Arduino monitor ===")
    fname = f"replay_{datetime.now().strftime('%Y%m%d_%H%M%S')}_v6_run{run}.csv"
    f1.reset_input_buffer()
    nlines = nrep = 0; t0 = time.time()
    with open(fname, "w", encoding="utf-8", newline="") as fh:
        print(f"  LOGGING {F1} -> {fname} for {DURATION}s ... send REPLAY_START now")
        while (time.time() - t0) < DURATION:
            line = f1.readline().decode(errors="replace")
            if not line: continue
            fh.write(line); nlines += 1
            if line.startswith("REPORT"): nrep += 1
            if nlines <= 5: print(f"  [F1 live] {line.strip()[:100]}")
            elif nlines % 2500 == 0:
                print(f"    {nlines} lines ({nrep} REPORT), t={time.time()-t0:6.1f}s")
    print(f"  RUN {run} done: {time.time()-t0:.1f}s | {nlines} lines | {nrep} REPORT | -> {fname}")
    if run < 3: print("  Wait >2 min bench re-settle before starting the next run.")
f1.close()
print("\nAll three captures complete.")
