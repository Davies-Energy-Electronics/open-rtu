
import serial, time

port = "COM6"  # Replace with Feather 1's COM port

ser = serial.Serial(port, 115200, timeout=1)

with open("smoke_v6.txt", "w", encoding="utf-8") as f:

    print(f"Logging from {port} to smoke_v6.txt for 30 seconds - press RST on Feather 1 now...")

    start = time.time()

    while time.time() - start < 30:

        line = ser.readline().decode("utf-8", errors="replace")

        if line:

            f.write(line)

            print(line, end="")

ser.close()

print("Done. smoke_v6.txt saved.")

