
import serial, time, datetime

port = "COM6"  # CHANGE to Feather 1 COM port

ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

outfile = f"replay_{ts}_parallel_raw.txt"

ser = serial.Serial(port, 115200, timeout=1)

print(f"Capturing to {outfile} — send REPLAY_START on Feather 2 now.")

print("Capture runs until you press Ctrl+C or 400s elapses.")

start = time.time()

report_count = 0

with open(outfile, "w", encoding="utf-8") as f:

    try:

        while time.time() - start < 400:  # 6.7 min max

            line = ser.readline().decode("utf-8", errors="replace")

            if line:

                f.write(line)

                f.flush()

                if "REPORT," in line: report_count += 1

                if report_count % 500 == 0 and report_count > 0:

                    elapsed = time.time() - start

                    print(f"  {elapsed:.0f}s: {report_count} REPORT lines captured")

    except KeyboardInterrupt:

        print(f"\nStopped by user. {report_count} REPORT lines captured.")

ser.close()

print(f"Done. File: {outfile}")

