/*
 * Feather 2 — Va/Vb DAC generator with NESO trajectory replay mode
 *
 * Normal mode: 50 Hz fixed (Va on DAC_VA at 0°, Vb on DAC_VB at -120°)
 * Replay mode: steps through 360 NESO 9 Aug 2019 frequency samples at 1 Hz,
 *              reprogramming the hardware timer alarm period at each step.
 *
 * Serial commands:
 *   REPLAY_START   — begin replay from index 0
 *   REPLAY_ABORT   — return to 50 Hz nominal immediately
 *
 * Serial output during replay:
 *   REPLAY_BEGIN
 *   REPLAY <i> <f_target>   (one line per second, 360 lines total)
 *   REPLAY_END
 */

#include <math.h>
#include "neso_trajectory.h"

#define DAC_VA     25
#define DAC_VB     26
#define LUT_SIZE   40
#define AMPLITUDE  100

static uint8_t lutVa[LUT_SIZE], lutVb[LUT_SIZE];
volatile int   dacIdx = 0;
hw_timer_t    *timer  = NULL;

// === Replay-mode state ===
enum ReplayMode { MODE_NORMAL = 0, MODE_REPLAY = 1 };
volatile ReplayMode currentMode = MODE_NORMAL;
volatile uint16_t   replayIndex = 0;
uint32_t            lastReplayStepMs = 0;
const float         REPLAY_NOMINAL_FREQ = 50.0f;
const uint32_t      REPLAY_STEP_MS      = 1000;  // 1 second per NESO sample

void IRAM_ATTR onTimer() {
  dacWrite(DAC_VA, lutVa[dacIdx]);
  dacWrite(DAC_VB, lutVb[dacIdx]);
  dacIdx = (dacIdx + 1) % LUT_SIZE;
}

// Reprogram the hardware-timer alarm period for the requested fundamental.
// Does NOT reset dacIdx — the sine waveform continues from its current phase.
void setHardwareTimerForFrequency(float freq) {
  // period_us = 1e6 / (freq * LUT_SIZE), rounded to nearest integer
  uint64_t period_us = (uint64_t)(1000000.0 / (freq * (float)LUT_SIZE) + 0.5);
  timerAlarm(timer, period_us, true, 0);
}

void handleReplayCommand() {
  if (!Serial.available()) return;
  String cmd = Serial.readStringUntil('\n');
  cmd.trim();
  if (cmd == "REPLAY_START" && currentMode == MODE_NORMAL) {
    currentMode = MODE_REPLAY;
    replayIndex = 0;
    lastReplayStepMs = millis() - REPLAY_STEP_MS;  // trigger first step now
    Serial.println("REPLAY_BEGIN");
  } else if (cmd == "REPLAY_ABORT" && currentMode == MODE_REPLAY) {
    currentMode = MODE_NORMAL;
    setHardwareTimerForFrequency(REPLAY_NOMINAL_FREQ);
    Serial.println("REPLAY_ABORTED");
  }
}

void setup() {
  Serial.begin(115200);
  delay(1000);

  for (int i = 0; i < LUT_SIZE; i++) {
    float a = 2.0f * PI * i / LUT_SIZE;
    lutVa[i] = (uint8_t)(127 + AMPLITUDE * sinf(a));
    lutVb[i] = (uint8_t)(127 + AMPLITUDE * sinf(a - 2.0f * PI / 3.0f));
  }

  dacWrite(DAC_VA, 127);
  dacWrite(DAC_VB, 127);
  Serial.println("Feather 2 — DC midpoint 5s then AC");
  delay(5000);

  timer = timerBegin(1000000);
  timerAttachInterrupt(timer, &onTimer);
  timerAlarm(timer, 500, true, 0);  // 500 us = 2 kHz ISR = 50 Hz sine

  Serial.println("AC running.");
  Serial.println("Send REPLAY_START to begin NESO 9 Aug 2019 replay.");
}

void loop() {
  handleReplayCommand();

  if (currentMode == MODE_REPLAY) {
    uint32_t now = millis();
    if (now - lastReplayStepMs >= REPLAY_STEP_MS) {
      lastReplayStepMs = now;
      if (replayIndex < NESO_TRAJECTORY_LENGTH) {
        float f_target = NESO_TRAJECTORY[replayIndex];
        setHardwareTimerForFrequency(f_target);
        Serial.print("REPLAY ");
        Serial.print(replayIndex);
        Serial.print(" ");
        Serial.println(f_target, 4);
        replayIndex++;
      } else {
        setHardwareTimerForFrequency(REPLAY_NOMINAL_FREQ);
        currentMode = MODE_NORMAL;
        Serial.println("REPLAY_END");
      }
    }
  }

  delay(10);  // yield to background tasks
}
