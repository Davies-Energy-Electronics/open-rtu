import serial, time, threading
from datetime import datetime

F1, F2 = "COM6", "COM5"
RUNS, GAP, TIMEOUT = 3, 150, 480

f1 = serial.Serial(F1, 115200, timeout=1)
f2 = serial.Serial(F2, 115200, timeout=1)
print("Ports open - waiting 5 s for DTR-reset settle...")
time.sleep(5); f1.reset_input_buffer(); f2.reset_input_buffer()

end_seen = threading.Event()
def watch_f2():
    while not end_seen.is_set():
        try: line = f2.readline().decode(errors="replace").strip()
        except Exception: break
        if line: print(f"  [F2] {line[:100]}")
        if "REPLAY_END" in line: end_seen.set()

for run in range(1, RUNS + 1):
    input(f"\n=== RUN {run}/{RUNS} - press Enter to start capture ===")
    fname = f"replay_{datetime.now().strftime('%Y%m%d_%H%M%S')}_v6_run{run}.csv"
    end_seen.clear(); f1.reset_input_buffer(); f2.reset_input_buffer()
    t = threading.Thread(target=watch_f2, daemon=True); t.start()
    nlines = nrep = 0; t0 = time.time()
    with open(fname, "w", encoding="utf-8", newline="") as fh:
        f2.write(b"REPLAY_START\r\n"); f2.flush()
        print(f"  REPLAY_START sent -> logging {F1} to {fname}")
        while not end_seen.is_set() and (time.time() - t0) < TIMEOUT:
            line = f1.readline().decode(errors="replace")
            if line:
                fh.write(line); nlines += 1
                if line.startswith("REPORT"): nrep += 1
                if nlines <= 5: print(f"  [F1 live] {line.strip()[:100]}")
                elif nlines % 2500 == 0:
                    print(f"    {nlines} lines ({nrep} REPORT), t={time.time()-t0:6.1f}s")
            elif nlines == 0 and (time.time() - t0) > 10:
                print("  !! 10 s and NOTHING from F1 - Ctrl+C, ports are wrong or F1 is dead"); t0 -= 5
            if "REPLAY_END" in line: end_seen.set()
    dur = time.time() - t0
    status = "REPLAY_END seen" if end_seen.is_set() else "TIMEOUT (no REPLAY_END!)"
    print(f"  RUN {run} done: {status} | {dur:.1f}s | {nlines} lines | {nrep} REPORT | -> {fname}")
    if run < RUNS:
        print(f"  Re-settle {GAP}s before next run...")
        for r in range(GAP, 0, -30): print(f"    {r}s remaining"); time.sleep(min(30, r))
f1.close(); f2.close()
print("\nAll runs complete.")
