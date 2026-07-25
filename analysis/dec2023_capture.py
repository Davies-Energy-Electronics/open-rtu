import serial, time, threading
from datetime import datetime

f1 = serial.Serial("COM6", 115200, timeout=1)

f2 = serial.Serial()
f2.port, f2.baudrate, f2.timeout = "COM5", 115200, 1
f2.dtr = False; f2.rts = False          # do NOT enter bootloader
f2.open()
f2.rts = True; time.sleep(0.2); f2.rts = False   # clean run-mode reset (EN pulse, IO0 high)
print("F2 reset to RUN mode - waiting 10 s for boot + DC midpoint + AC start...")
t0 = time.time()
while time.time() - t0 < 12:
    ln = f2.readline().decode(errors="replace").strip()
    if ln: print(f"  [F2 boot] {ln[:80]}")
    if "AC running" in ln: break
f1.reset_input_buffer(); f2.reset_input_buffer()

end = threading.Event()
def watch():
    while not end.is_set():
        try: ln = f2.readline().decode(errors="replace").strip()
        except Exception: break
        if ln and (ln.startswith("REPLAY_") or int(ln.split()[1])%60==0 if ln.startswith("REPLAY ") else True):
            print(f"  [F2] {ln[:60]}")
        if "REPLAY_END" in ln: end.set()
threading.Thread(target=watch, daemon=True).start()

fname = f"replay_{datetime.now().strftime('%Y%m%d_%H%M%S')}_dec2023peak_v6.csv"
import re
pat = re.compile(r"f_Vb=([0-9.]+)")
recent = []; nr = 0; t0 = time.time()
with open(fname, "w", encoding="utf-8", newline="") as fh:
    f2.write(b"REPLAY_START\r\n"); f2.flush()
    print(f"REPLAY_START sent -> logging COM6 to {fname}")
    while not end.is_set() and (time.time()-t0) < 420:
        line = f1.readline().decode(errors="replace")
        if not line: continue
        fh.write(line)
        if line.startswith("REPORT"):
            nr += 1
            m = pat.search(line)
            if m: recent.append(float(m.group(1))); recent = recent[-100:]
            if nr == 400:                      # ~8 s in: sanity gate
                import statistics
                med = statistics.median(recent)
                if 45 <= med <= 55:
                    print(f"  SANITY OK: median f_Vb {med:.3f} Hz - healthy sine, continuing")
                else:
                    print(f"  !! SANITY FAIL: median f_Vb {med:.1f} Hz - F2 not generating. Ctrl+C and shout.")
            if nr % 2500 == 0: print(f"  [F1] {nr} REPORT, t={time.time()-t0:5.1f}s")
        if "REPLAY_END" in line: end.set()
f1.close(); f2.close()
print(f"Done: {'REPLAY_END seen' if end.is_set() else 'TIMEOUT'} | {nr} REPORT -> {fname}")
