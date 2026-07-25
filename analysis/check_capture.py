import csv
rows = list(csv.DictReader(open("replay_20260714_174352_3phase_thd_varied.csv")))
mono_prev = float(rows[0]["monotonic_s"])
idx_prev = int(rows[0]["replay_index"])
gaps = []
for r in rows[1:]:
    m = float(r["monotonic_s"]); i = int(r["replay_index"])
    if (m - mono_prev) > 3.0 or (i - idx_prev) > 3:
        gaps.append((mono_prev, m, m-mono_prev, idx_prev, i, i-idx_prev))
    mono_prev, idx_prev = m, i
total_s = float(rows[-1]["monotonic_s"]) - float(rows[0]["monotonic_s"])
print("Total rows:", len(rows))
print("Duration: {:.1f}s".format(total_s))
print("Sleep/skip gaps > 3s or 3 indices:", len(gaps))
for g in gaps:
    print("  {:.2f}s -> {:.2f}s (dt={:.2f}s), idx {} -> {} (di={})".format(*g))
print("PASS" if len(gaps)==0 else "FAIL - has gaps")
