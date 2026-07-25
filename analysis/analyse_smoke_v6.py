import re
lines = open("smoke_v6.txt").readlines()
report_lines = [l for l in lines if "REPORT," in l]
diag_lines = [l for l in lines if "DIAG," in l]
print(f"REPORT lines: {len(report_lines)}")
print(f"DIAG lines:   {len(diag_lines)}")
if report_lines:
    times = [int(re.search(r"REPORT,(\d+),", l).group(1)) for l in report_lines]
    deltas = [times[i+1] - times[i] for i in range(len(times)-1)]
    mean_d = sum(deltas)/len(deltas) if deltas else 0
    print(f"Cadence: n={len(deltas)}, mean={mean_d:.1f} ms, min={min(deltas) if deltas else 0}, max={max(deltas) if deltas else 0}")
    f_Vb_vals = []
    for l in report_lines:
        m = re.search(r"f_Vb=([-\d.]+)", l)
        if m: f_Vb_vals.append(float(m.group(1)))
    if f_Vb_vals:
        nonzero = [f for f in f_Vb_vals if f > 0.1]
        print(f"f_Vb: n={len(f_Vb_vals)}, nonzero={len(nonzero)}, mean_nonzero={sum(nonzero)/len(nonzero) if nonzero else 0:.4f} Hz")
if diag_lines:
    last_diag = diag_lines[-1]
    print(f"Last DIAG: {last_diag.strip()}")
# Overall verdict
if report_lines and mean_d and 15 <= mean_d <= 30 and any(f > 0.1 for f in f_Vb_vals):
    print("V6 SMOKE PASS: 20 ms cadence + real f_Vb. Proceed to B.4.")
elif not report_lines:
    print("V6 SMOKE FAIL: no REPORT lines. Check banner in serial output.")
elif all(f < 0.1 for f in f_Vb_vals):
    print("V6 SMOKE PARTIAL: cadence OK, but f_Vb still zero. M2 port may need work.")
elif mean_d > 30:
    print(f"V6 SMOKE PARTIAL: cadence {mean_d:.1f} ms slower than 20 ms.")
else:
    print("V6 SMOKE PARTIAL: examine smoke_v6.txt manually.")
