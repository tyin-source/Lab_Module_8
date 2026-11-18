// Phys 39 Module 8 - Angstrom-method rod controller (Arduino side)
// (extends the Module 7 sketch, arduino/module_07_rod_controller in the Module 7
// repository, which stays unchanged as the working checkpoint)
//
// Changes from Module 7:
//   - The Arduino now generates the drive itself. The sinusoidal setpoint, the
//     raw sinusoidal PWM, the feedforward command, and the P/PI law are all
//     computed here, once per DATA record (about 1 s), from the same Arduino
//     clock that timestamps the record. The Python GUI only sends settings and a
//     heartbeat, and displays and saves the records. This keeps the requested
//     waveform, the applied PWM, and the six temperatures on one time base.
//   - The phase of the sinusoid is accumulated, phase += 2 pi dt / period, so a
//     period change during a run continues the waveform without a jump.
//   - Five control modes (control_mode field):
//       MANUAL    signed_pwm = manual PWM (constant; 0 for a baseline)
//       PWM_SINE  signed_pwm = u_a + u_b sin(phase)   raw open-loop PWM (S18.3-S18.4):
//                                                     u_a = 0 zero-centered, u_a > u_b positive-only
//                                                     (no cooling compensation in this mode)
//       FF        u = (Tset - T_room) / chi_H          static feedforward (S18.5-S18.6);
//                                                     = u_a + u_b sin(phase) with
//                                                     u_a = (Ta - T_room)/chi_H, u_b = Tb/chi_H
//       P         u = Kp e                             e = Tset - T_TEC (S18.7)
//       PI        u = Kp e + Ki q,  q = sum of e dt    conditional-integration anti-windup
//     with Tset = Ta + Tb sin(phase). In FF, P, and PI every negative (cooling)
//     command is multiplied by the measured heating-to-cooling ratio r (static
//     actuator compensation, S18.5 and S18.7). With FFADD 1 the feedforward term
//     is added to the P/PI command before that scaling (optional).
//     The result is clamped to +/- u_max and rounded.
//   - One DATA record per average (header printed once at start-up):
//       DATA,time_ms,run_state,control_mode,step,Ta_C,Tb_C,period_s,setpoint_C,
//            feedforward_pwm,TEC_C,T1_C,T2_C,T3_C,T4_C,T5_C,controller_pwm,
//            heating_to_cooling_ratio,signed_pwm
//     (one line; wrapped here). Fields that do not apply in the current mode
//     print nan; the number and order of fields never change. See
//     data/module_08/README.md for every field and its placeholder.
//   - Every transition (start, stop, mode or setting change, safety, timeout)
//     increments `step` and prints one line
//       EVENT,<time_ms>,<step>,<description>
//     The GUI copies EVENT lines into the run's metadata file.
//   - Safety: the Module 7 limit on the TEC thermistor stays, and a second limit
//     on the rod-base thermistor T1 is added. Either one (or a missing TEC or T1
//     reading) during a run zeros both H-bridge outputs and sets run_state to
//     SAFETY. A run also stops (run_state TIMEOUT) if no command, including the
//     GUI's once-per-record heartbeat HB, arrives for COMMAND_TIMEOUT_MS. After
//     SAFETY or TIMEOUT the TEC stays off until a new START, which is refused
//     while a limit is still exceeded.
//   - Serial speed raised from 9600 to 115200 baud: a record is now about 130
//     characters, and at 9600 baud printing it would hold up the ADC sampling.
//
// Commands from the Python GUI (one per line, keywords case-insensitive):
//   HB                        heartbeat (no other effect)
//   START | STOP              start the selected mode / zero the output
//   MODE MANUAL|PWM_SINE|FF|P|PI
//   SET MEAN <C>              Ta, mean of the requested TEC temperature
//   SET AMP <C>               Tb, amplitude of the requested TEC temperature
//   SET PERIOD <s>            forcing period
//   SET KP <PWM/C> | SET KI <PWM/(C s)> | SET UMAX <0-255>
//   SET RATIO <r>             heating-to-cooling ratio chi_H / chi_C
//   SET CHIH <C/PWM>          heating susceptibility used by the feedforward
//   SET TROOM <C>             room temperature used by the feedforward
//   SET PWMAMP <counts> | SET PWMOFF <counts>   u_b (amplitude) and u_a (offset) of PWM_SINE
//   SET MANUAL <counts>       signed PWM of MANUAL (cooling negative)
//   SET FFADD 0|1             add the feedforward term to P/PI
//   SET T1LIMIT <C>           rod-base safety limit in force, at most 56 C
//   SET TECLIMIT <C>          TEC safety limit in force, at most 60 C
//                             (lowered only to test the shutdown; both return to
//                             56 and 60 C at every power-up or reset)
//   ZERO                      zero the PI integral q
//   GET                       print the current settings (CONFIG line)
// A malformed command or an out-of-range value is rejected with an ERROR line
// and changes nothing (it does not count as a heartbeat either).
//
// Pin 9 / Pin 10: H-bridge control signals (unchanged since Module 4).

#include <math.h>

// ------------------------------------------------------------------
// CHANNELS (unchanged from Module 7)
// ------------------------------------------------------------------
const int NUM_CHANNELS = 6;
// Index 0 is the TEC-base thermistor; indices 1-5 are rod thermistors T1-T5.
const int THERMISTOR_PINS[NUM_CHANNELS] = {A0, A1, A2, A3, A4, A5};
const int TEC_CHANNEL = 0;
const int ROD_BASE_CHANNEL = 1;   // T1, first rod hole, x = 0

const int HBRIDGE_PIN_1 = 9;
const int HBRIDGE_PIN_2 = 10;

// 500 readings of all six channels, one every 2000 us: one average and one
// DATA record about every second, as in Modules 6 and 7.
const int ADC_SAMPLES = 500;
const unsigned long SAMPLE_INTERVAL_US = 2000;

// Thermistor constants (Module 2) and divider resistors (Module 7: all six within
// 0.3 % of 100 kOhm, so the nominal value is kept for every channel).
const float SERIES_RESISTOR_OHM[NUM_CHANNELS] = {
  100000.0,  // A0  TEC
  100000.0,  // A1  T1 (rod, x = 0, first hole)
  100000.0,  // A2  T2
  100000.0,  // A3  T3
  100000.0,  // A4  T4
  100000.0   // A5  T5
};
const float NOMINAL_RESISTANCE = 100000.0;
const float NOMINAL_TEMPERATURE_C = 25.0;
const float BETA_COEFFICIENT = 4540.0;
const float ADC_MAX = 1023.0;

// ------------------------------------------------------------------
// SAFETY LIMITS
// ------------------------------------------------------------------
// TEC-base thermistor: the Module 4-7 software limit. The hardware thermal
// switch (about 70 C) in series with the TEC remains the final protection.
const float TEC_LIMIT_C = 60.0;
// Rod-base thermistor T1: S18.1 requires the measured rod-base temperature to stay
// at or below 56 C. Keep margin for overshoot in the requested waveform.
const float ROD_BASE_LIMIT_C = 56.0;
// The limits in force. SET T1LIMIT / SET TECLIMIT can lower them (never raise them
// above the constants above) so the shutdown can be tested without uploading a
// changed sketch: the schedule's TEST_T1_LIMIT and TEST_TEC_LIMIT rows set a limit
// 2 C above the present reading, heat gently, and check that the run stops.
float rodBaseLimitC = ROD_BASE_LIMIT_C;
float tecLimitC = TEC_LIMIT_C;

// Limits on the requested waveform, checked when a setting arrives.
const float SETPOINT_MIN_C = 10.0;    // Ta - Tb may not go below this
const float SETPOINT_MAX_C = 56.0;    // Ta + Tb may not go above this (40 + 16 C, S18.5)
const float MIN_PERIOD_S = 100.0;     // no forcing period shorter than 100 s (S18.1)
const float MAX_PERIOD_S = 3600.0;

// Stop a run if no command (the GUI sends HB after every record) arrives for
// this long: Python crashed, GUI frozen, or USB unplugged.
const unsigned long COMMAND_TIMEOUT_MS = 3000;

// ------------------------------------------------------------------
// EXPERIMENTALLY VERIFIED HEAT/COOL MAPPING (unchanged since Module 4)
// ------------------------------------------------------------------
const int HEAT_ACTIVE_PIN = 9;
const int COOL_ACTIVE_PIN = 10;

// ------------------------------------------------------------------
// RUN STATES AND CONTROL MODES
// ------------------------------------------------------------------
// Plain int constants instead of enums: the Arduino IDE writes function
// prototypes above the sketch, where an enum type would not yet be known.
const int STATE_OFF = 0;
const int STATE_RUN = 1;
const int STATE_SAFETY = 2;
const int STATE_TIMEOUT = 3;
const char* const STATE_NAMES[] = {"OFF", "RUN", "SAFETY", "TIMEOUT"};

const int MODE_MANUAL = 0;
const int MODE_PWM_SINE = 1;
const int MODE_FF = 2;
const int MODE_P = 3;
const int MODE_PI = 4;
const int NUM_MODES = 5;
const char* const MODE_NAMES[] = {"MANUAL", "PWM_SINE", "FF", "P", "PI"};

// ------------------------------------------------------------------
// SETTINGS (defaults; the GUI sends its own values when it connects)
// ------------------------------------------------------------------
int controlMode = MODE_MANUAL;
float meanC = 40.0;          // Ta (C)   S18.5 suggested value
float ampC = 16.0;           // Tb (C)   S18.5 suggested value
float periodS = 100.0;       // forcing period (s)
float kp = 6.0;              // PWM per C       (Module 7 Run 2)
float ki = 0.04;             // PWM per (C s)   (Module 7 Run 2)
float ratioHC = 2.82;        // chi_H / chi_C   (Module 4)
float chiH = 0.485;          // C per PWM count (Module 4)
float roomC = 22.0;          // T_room for the feedforward (C)
float pwmAmp = 50.0;         // u_b of PWM_SINE (counts), S18.3 initial value
float pwmOffset = 0.0;       // u_a of PWM_SINE (counts)
float manualPwm = 0.0;       // signed PWM of MANUAL (counts)
float uMax = 255.0;          // PWM limit (counts)
bool ffAdd = false;          // add the feedforward term to P/PI

// ------------------------------------------------------------------
// STATE
// ------------------------------------------------------------------
int runState = STATE_OFF;
int step = 0;                // run-segment counter, +1 at every transition
float phase = 0.0;           // rad, accumulated
float integralQ = 0.0;       // C s
unsigned long lastControlMs = 0;
unsigned long lastCommandMs = 0;

float temperaturesC[NUM_CHANNELS];

// Values of the most recent control update, printed in the DATA record.
float recSetpoint = NAN;
float recFeedforward = NAN;
float recController = NAN;
float recRatio = NAN;
int recSignedPwm = 0;

unsigned long adcSum[NUM_CHANNELS];
int adcCount = 0;
unsigned long lastSampleUs = 0;

String commandBuffer = "";

// ------------------------------------------------------------------
// THERMISTOR MEASUREMENT (unchanged from Module 7)
// ------------------------------------------------------------------
float adcToTemperatureC(float adc, float seriesResistor) {
  if (adc <= 0 || adc >= ADC_MAX) {
    return NAN;
  }

  float resistance = seriesResistor * adc / (ADC_MAX - adc);

  float steinhart = resistance / NOMINAL_RESISTANCE;
  steinhart = log(steinhart);
  steinhart /= BETA_COEFFICIENT;
  steinhart += 1.0 / (NOMINAL_TEMPERATURE_C + 273.15);
  steinhart = 1.0 / steinhart;
  steinhart -= 273.15;

  return steinhart;
}

// Read all six channels once every SAMPLE_INTERVAL_US (each twice, keeping the
// second reading: the 100 kOhm dividers are high-impedance sources). After
// ADC_SAMPLES readings, update temperaturesC and return true.
bool takeSample() {
  unsigned long nowUs = micros();
  if (nowUs - lastSampleUs < SAMPLE_INTERVAL_US) {
    return false;
  }
  lastSampleUs = nowUs;

  for (int ch = 0; ch < NUM_CHANNELS; ch++) {
    analogRead(THERMISTOR_PINS[ch]);               // settling reading, discarded
    adcSum[ch] += analogRead(THERMISTOR_PINS[ch]);  // kept
  }
  adcCount++;

  if (adcCount < ADC_SAMPLES) {
    return false;
  }

  for (int ch = 0; ch < NUM_CHANNELS; ch++) {
    temperaturesC[ch] = adcToTemperatureC(adcSum[ch] / (float)ADC_SAMPLES,
                                          SERIES_RESISTOR_OHM[ch]);
    adcSum[ch] = 0;
  }
  adcCount = 0;
  return true;
}

// ------------------------------------------------------------------
// H-BRIDGE OUTPUT
// ------------------------------------------------------------------
void zeroOutput() {
  digitalWrite(HBRIDGE_PIN_1, LOW);
  digitalWrite(HBRIDGE_PIN_2, LOW);
  analogWrite(HBRIDGE_PIN_1, 0);
  analogWrite(HBRIDGE_PIN_2, 0);
}

// Apply a signed PWM: positive heats, negative cools. Only one H-bridge input
// gets PWM at a time; the other is held LOW first, so both are never active.
void applySignedPwm(int u) {
  if (runState != STATE_RUN || u == 0) {
    zeroOutput();
    return;
  }
  int pwmPin = (u > 0) ? HEAT_ACTIVE_PIN : COOL_ACTIVE_PIN;
  int lowPin = (u > 0) ? COOL_ACTIVE_PIN : HEAT_ACTIVE_PIN;
  int magnitude = (u > 0) ? u : -u;
  if (magnitude > 255) magnitude = 255;

  analogWrite(lowPin, 0);
  digitalWrite(lowPin, LOW);
  analogWrite(pwmPin, magnitude);
}

// ------------------------------------------------------------------
// OUTPUT LINES
// ------------------------------------------------------------------
void printValue(float value, int decimals) {
  if (isnan(value)) {
    Serial.print(F("nan"));
  } else {
    Serial.print(value, decimals);
  }
}

void printEventStart() {
  step++;
  Serial.print(F("EVENT,"));
  Serial.print(millis());
  Serial.print(F(","));
  Serial.print(step);
  Serial.print(F(","));
}

// All settings on one line (used by GET and in the START event).
void printSettings() {
  Serial.print(F("mode="));     Serial.print(MODE_NAMES[controlMode]);
  Serial.print(F(" Ta_C="));    Serial.print(meanC, 2);
  Serial.print(F(" Tb_C="));    Serial.print(ampC, 2);
  Serial.print(F(" period_s="));Serial.print(periodS, 1);
  Serial.print(F(" Kp="));      Serial.print(kp, 3);
  Serial.print(F(" Ki="));      Serial.print(ki, 4);
  Serial.print(F(" u_max="));   Serial.print(uMax, 0);
  Serial.print(F(" r="));       Serial.print(ratioHC, 3);
  Serial.print(F(" chi_H="));   Serial.print(chiH, 4);
  Serial.print(F(" T_room_C="));Serial.print(roomC, 2);
  Serial.print(F(" u_a="));     Serial.print(pwmOffset, 1);
  Serial.print(F(" u_b="));     Serial.print(pwmAmp, 1);
  Serial.print(F(" manual_pwm="));Serial.print(manualPwm, 0);
  Serial.print(F(" ff_add="));  Serial.print(ffAdd ? 1 : 0);
  Serial.print(F(" T1_limit_C="));  Serial.print(rodBaseLimitC, 2);
  Serial.print(F(" TEC_limit_C=")); Serial.print(tecLimitC, 2);
}

void printDataRecord() {
  bool hasSetpoint = (controlMode == MODE_FF || controlMode == MODE_P ||
                      controlMode == MODE_PI);
  bool hasPeriod = (controlMode != MODE_MANUAL);

  Serial.print(F("DATA,"));
  Serial.print(millis());                                  // time_ms
  Serial.print(F(","));
  Serial.print(STATE_NAMES[runState]);                     // run_state
  Serial.print(F(","));
  Serial.print(MODE_NAMES[controlMode]);                   // control_mode
  Serial.print(F(","));
  Serial.print(step);                                      // step
  Serial.print(F(","));
  printValue(hasSetpoint ? meanC : NAN, 2);                // Ta_C
  Serial.print(F(","));
  printValue(hasSetpoint ? ampC : NAN, 2);                 // Tb_C
  Serial.print(F(","));
  printValue(hasPeriod ? periodS : NAN, 1);                // period_s
  Serial.print(F(","));
  printValue(recSetpoint, 2);                              // setpoint_C
  Serial.print(F(","));
  printValue(recFeedforward, 1);                           // feedforward_pwm
  for (int ch = 0; ch < NUM_CHANNELS; ch++) {
    Serial.print(F(","));
    printValue(temperaturesC[ch], 2);                      // TEC_C, T1_C ... T5_C
  }
  Serial.print(F(","));
  printValue(recController, 1);                            // controller_pwm
  Serial.print(F(","));
  printValue(recRatio, 3);                                 // heating_to_cooling_ratio
  Serial.print(F(","));
  Serial.println(recSignedPwm);                            // signed_pwm
}

// ------------------------------------------------------------------
// SAFETY
// ------------------------------------------------------------------
// Returns 0 if safe, otherwise a code for the first problem found.
int unsafeReason() {
  float tec = temperaturesC[TEC_CHANNEL];
  float t1 = temperaturesC[ROD_BASE_CHANNEL];
  if (isnan(tec)) return 1;
  if (tec > tecLimitC) return 2;
  if (isnan(t1)) return 3;
  if (t1 > rodBaseLimitC) return 4;
  return 0;
}

void printUnsafeReason(int reason) {
  switch (reason) {
    case 1: Serial.print(F("no valid TEC reading (check A0)")); break;
    case 2: Serial.print(F("TEC "));
            Serial.print(temperaturesC[TEC_CHANNEL], 2);
            Serial.print(F(" C above the "));
            Serial.print(tecLimitC, 2);
            Serial.print(F(" C limit")); break;
    case 3: Serial.print(F("no valid T1 reading (check A1)")); break;
    case 4: Serial.print(F("T1 "));
            Serial.print(temperaturesC[ROD_BASE_CHANNEL], 2);
            Serial.print(F(" C above the "));
            Serial.print(rodBaseLimitC, 2);
            Serial.print(F(" C rod-base limit")); break;
    default: Serial.print(F("unknown")); break;
  }
}

// Called with every new set of averages, before the control update.
void checkSafety() {
  if (runState != STATE_RUN) {
    return;                    // output is already zero
  }
  int reason = unsafeReason();
  if (reason == 0) {
    return;
  }
  runState = STATE_SAFETY;
  zeroOutput();
  printEventStart();
  Serial.print(F("SAFETY SHUTDOWN: "));
  printUnsafeReason(reason);
  Serial.println(F(". Both H-bridge outputs set to 0; send START to resume."));
}

// ------------------------------------------------------------------
// CONTROL UPDATE (once per DATA record)
// ------------------------------------------------------------------
float scaleCooling(float u) {
  // Static actuator compensation: cooling commands are multiplied by r.
  return (u < 0) ? u * ratioHC : u;
}

float clampU(float u) {
  if (u > uMax) return uMax;
  if (u < -uMax) return -uMax;
  return u;
}

void updateControl() {
  recSetpoint = NAN;
  recFeedforward = NAN;
  recController = NAN;
  recRatio = NAN;
  recSignedPwm = 0;

  if (runState != STATE_RUN) {
    zeroOutput();
    return;
  }

  unsigned long nowMs = millis();
  float dt = (nowMs - lastControlMs) / 1000.0;
  lastControlMs = nowMs;
  if (controlMode != MODE_MANUAL) {
    phase += TWO_PI * dt / periodS;
    if (phase >= TWO_PI) {
      phase = fmod(phase, TWO_PI);
    }
  }
  float s = sin(phase);
  float u = 0.0;

  if (controlMode == MODE_MANUAL) {
    u = manualPwm;
  } else if (controlMode == MODE_PWM_SINE) {
    u = pwmOffset + pwmAmp * s;
  } else {
    float tset = meanC + ampC * s;
    float uff = (tset - roomC) / chiH;          // before the cooling scaling
    recSetpoint = tset;
    recRatio = ratioHC;

    if (controlMode == MODE_FF) {
      recFeedforward = uff;
      u = scaleCooling(uff);
    } else {
      float e = tset - temperaturesC[TEC_CHANNEL];
      float ffPart = ffAdd ? uff : 0.0;
      recFeedforward = ffAdd ? uff : NAN;

      if (controlMode == MODE_PI) {
        // Conditional integration (Module 6): keep q if the updated command
        // would be beyond the limit and this error pushes it farther.
        float qNext = integralQ + e * dt;
        float uNext = scaleCooling(ffPart + kp * e + ki * qNext);
        bool hold = (uNext > uMax && e > 0) || (uNext < -uMax && e < 0);
        if (!hold) {
          integralQ = qNext;
        }
      } else {
        integralQ = 0.0;
      }
      float uc = kp * e + ki * integralQ;
      recController = uc;
      u = scaleCooling(ffPart + uc);
    }
  }

  u = clampU(u);
  recSignedPwm = (int)(u >= 0 ? u + 0.5 : u - 0.5);
  applySignedPwm(recSignedPwm);
}

// ------------------------------------------------------------------
// COMMAND PARSER
// ------------------------------------------------------------------
// A number: optional leading '-', digits, at most one '.'.
bool parseNumber(String text, float &value) {
  text.trim();
  if (text.length() == 0) return false;
  int dots = 0;
  int digits = 0;
  for (unsigned int i = 0; i < text.length(); i++) {
    char c = text.charAt(i);
    if (c == '.') {
      dots++;
    } else if (isDigit(c)) {
      digits++;
    } else if (!(c == '-' && i == 0)) {
      return false;
    }
  }
  if (dots > 1 || digits == 0) return false;
  value = text.toFloat();
  return true;
}

void reject(const __FlashStringHelper *why, String line) {
  Serial.print(F("ERROR: "));
  Serial.print(why);
  Serial.print(F(": "));
  Serial.println(line);
}

bool waveformAllowed(float mean, float amp) {
  return amp >= 0 && mean - amp >= SETPOINT_MIN_C && mean + amp <= SETPOINT_MAX_C;
}

// Returns true if the command was valid (and therefore counts as a heartbeat).
bool handleSet(String name, String valueText, String line) {
  float v;
  if (!parseNumber(valueText, v)) {
    reject(F("not a number"), line);
    return false;
  }
  float old;
  if (name == "MEAN") {
    if (!waveformAllowed(v, ampC)) { reject(F("Ta - Tb or Ta + Tb outside the setpoint limits"), line); return false; }
    old = meanC; meanC = v;
  } else if (name == "AMP") {
    if (!waveformAllowed(meanC, v)) { reject(F("Ta - Tb or Ta + Tb outside the setpoint limits"), line); return false; }
    old = ampC; ampC = v;
  } else if (name == "PERIOD") {
    if (v < MIN_PERIOD_S || v > MAX_PERIOD_S) { reject(F("period outside the allowed range"), line); return false; }
    old = periodS; periodS = v;          // phase is continuous: no jump
  } else if (name == "KP") {
    if (v < 0 || v > 100) { reject(F("Kp outside 0-100"), line); return false; }
    old = kp; kp = v;
  } else if (name == "KI") {
    if (v < 0 || v > 5) { reject(F("Ki outside 0-5"), line); return false; }
    old = ki; ki = v;
  } else if (name == "UMAX") {
    if (v < 0 || v > 255) { reject(F("u_max outside 0-255"), line); return false; }
    old = uMax; uMax = v;
  } else if (name == "RATIO") {
    if (v < 1 || v > 10) { reject(F("r outside 1-10"), line); return false; }
    old = ratioHC; ratioHC = v;
  } else if (name == "CHIH") {
    if (v < 0.01 || v > 5) { reject(F("chi_H outside 0.01-5"), line); return false; }
    old = chiH; chiH = v;
  } else if (name == "TROOM") {
    if (v < 0 || v > 40) { reject(F("T_room outside 0-40"), line); return false; }
    old = roomC; roomC = v;
  } else if (name == "PWMAMP") {
    if (v < 0 || v > 255) { reject(F("u_b outside 0-255"), line); return false; }
    old = pwmAmp; pwmAmp = v;
  } else if (name == "PWMOFF") {
    if (v < -255 || v > 255) { reject(F("u_a outside -255-255"), line); return false; }
    old = pwmOffset; pwmOffset = v;
  } else if (name == "MANUAL") {
    if (v < -255 || v > 255) { reject(F("manual PWM outside -255-255"), line); return false; }
    old = manualPwm; manualPwm = v;
  } else if (name == "T1LIMIT") {
    if (v < 0 || v > ROD_BASE_LIMIT_C) { reject(F("T1 limit outside 0-56"), line); return false; }
    old = rodBaseLimitC; rodBaseLimitC = v;
  } else if (name == "TECLIMIT") {
    if (v < 0 || v > TEC_LIMIT_C) { reject(F("TEC limit outside 0-60"), line); return false; }
    old = tecLimitC; tecLimitC = v;
  } else if (name == "FFADD") {
    if (v != 0 && v != 1) { reject(F("FFADD must be 0 or 1"), line); return false; }
    old = ffAdd ? 1 : 0; ffAdd = (v == 1);
  } else {
    reject(F("unknown setting"), line);
    return false;
  }
  printEventStart();
  Serial.print(F("SET "));
  Serial.print(name);
  Serial.print(F(" "));
  Serial.print(v, 4);
  Serial.print(F(" (was "));
  Serial.print(old, 4);
  Serial.println(F(")"));
  return true;
}

bool handleCommand(String line) {
  line.trim();
  line.toUpperCase();
  if (line.length() == 0) return false;

  if (line == "HB") {
    return true;
  }
  if (line == "GET") {
    Serial.print(F("CONFIG,"));
    Serial.print(millis());
    Serial.print(F(","));
    Serial.print(step);
    Serial.print(F(","));
    printSettings();
    Serial.println();
    return true;
  }
  if (line == "START") {
    if (runState == STATE_RUN) {
      reject(F("already running"), line);
      return false;
    }
    int reason = unsafeReason();
    if (reason != 0) {
      Serial.print(F("ERROR: START refused: "));
      printUnsafeReason(reason);
      Serial.println();
      return false;
    }
    runState = STATE_RUN;
    phase = 0.0;
    integralQ = 0.0;
    lastControlMs = millis();
    printEventStart();
    Serial.print(F("START "));
    printSettings();
    Serial.println();
    return true;
  }
  if (line == "STOP") {
    runState = STATE_OFF;
    zeroOutput();
    printEventStart();
    Serial.print(F("STOP (output 0); q = "));
    Serial.println(integralQ, 2);
    return true;
  }
  if (line == "ZERO") {
    printEventStart();
    Serial.print(F("ZERO integral (q was "));
    Serial.print(integralQ, 2);
    Serial.println(F(")"));
    integralQ = 0.0;
    return true;
  }
  if (line.startsWith("MODE ")) {
    String name = line.substring(5);
    name.trim();
    for (int m = 0; m < NUM_MODES; m++) {
      if (name == MODE_NAMES[m]) {
        int old = controlMode;
        controlMode = m;
        integralQ = 0.0;
        printEventStart();
        Serial.print(F("MODE "));
        Serial.print(MODE_NAMES[m]);
        Serial.print(F(" (was "));
        Serial.print(MODE_NAMES[old]);
        Serial.println(F("; q = 0)"));
        return true;
      }
    }
    reject(F("unknown mode"), line);
    return false;
  }
  if (line.startsWith("SET ")) {
    String rest = line.substring(4);
    rest.trim();
    int space = rest.indexOf(' ');
    if (space < 0) {
      reject(F("SET needs a name and a value"), line);
      return false;
    }
    return handleSet(rest.substring(0, space), rest.substring(space + 1), line);
  }
  reject(F("unknown command"), line);
  return false;
}

// ------------------------------------------------------------------
// SETUP
// ------------------------------------------------------------------
void setup() {
  pinMode(HBRIDGE_PIN_1, OUTPUT);
  pinMode(HBRIDGE_PIN_2, OUTPUT);
  zeroOutput();                          // safety: PWM starts at zero

  Serial.begin(115200);
  commandBuffer.reserve(64);
  for (int ch = 0; ch < NUM_CHANNELS; ch++) {
    adcSum[ch] = 0;
    temperaturesC[ch] = NAN;
  }

  Serial.println(F("Module 8 Angstrom rod controller (drive and P/PI computed on the Arduino)"));
  Serial.print(F("SAFETY: TEC limit (C): "));
  Serial.println(TEC_LIMIT_C, 1);
  Serial.print(F("SAFETY: rod-base T1 limit (C): "));
  Serial.println(ROD_BASE_LIMIT_C, 1);
  Serial.print(F("SAFETY: command timeout (s): "));
  Serial.println(COMMAND_TIMEOUT_MS / 1000.0, 1);
  Serial.print(F("LIMITS: setpoint "));
  Serial.print(SETPOINT_MIN_C, 1);
  Serial.print(F("-"));
  Serial.print(SETPOINT_MAX_C, 1);
  Serial.print(F(" C, period >= "));
  Serial.print(MIN_PERIOD_S, 0);
  Serial.println(F(" s"));

  // Header, printed once. Every later DATA line is one record.
  Serial.println(F("DATA,time_ms,run_state,control_mode,step,Ta_C,Tb_C,period_s,setpoint_C,"
                   "feedforward_pwm,TEC_C,T1_C,T2_C,T3_C,T4_C,T5_C,controller_pwm,"
                   "heating_to_cooling_ratio,signed_pwm"));

  // Build the first average before accepting commands, so START always has
  // measured temperatures to check.
  while (!takeSample()) {
  }
  lastCommandMs = millis();
}

// ------------------------------------------------------------------
// LOOP
// ------------------------------------------------------------------
void loop() {
  // ---- incoming commands, one line at a time ----
  while (Serial.available() > 0) {
    char c = (char)Serial.read();
    if (c == '\n' || c == '\r') {
      if (commandBuffer.length() > 0) {
        if (handleCommand(commandBuffer)) {
          lastCommandMs = millis();
        }
        commandBuffer = "";
      }
    } else if (commandBuffer.length() < 64) {
      commandBuffer += c;
    } else {
      commandBuffer = "";
      Serial.println(F("ERROR: command line too long, ignored"));
    }
  }

  // ---- stop the run if the GUI stops sending commands ----
  if (runState == STATE_RUN && millis() - lastCommandMs > COMMAND_TIMEOUT_MS) {
    runState = STATE_TIMEOUT;
    zeroOutput();
    printEventStart();
    Serial.print(F("TIMEOUT: no command for "));
    Serial.print(COMMAND_TIMEOUT_MS / 1000.0, 1);
    Serial.println(F(" s. Both H-bridge outputs set to 0; send START to resume."));
  }

  // ---- one new set of averages: safety, control, record ----
  if (takeSample()) {
    checkSafety();
    updateControl();
    printDataRecord();
  }
}
