// Host test of arduino/module_08_angstrom_controller (no Arduino needed).
//
// The sketch is compiled together with tests/firmware_host/arduino_mock.h and run
// against a simulated TEC + rod:
//   TEC   one lump, dT/dt = (T_room + chi u - T) / tau, with chi = chi_H for
//         heating and chi_C for cooling (deliberately not the Module 4 values,
//         to stand in for the TEC loaded by the rod)
//   rod   five nodes, each relaxing toward the previous one and losing heat to
//         the room, so the oscillation shrinks and lags along the rod
// The harness sends commands as the GUI does: HB after every record, other commands
// one every 60 ms. The serial port is modelled at 115200 baud with the Uno's 63-byte
// receive and 64-byte transmit buffers (arduino_mock.h). It checks the record
// format, every control mode, the cooling scaling, phase continuity at a period
// change, the safety limits, START refusal, the heartbeat timeout, that a burst of
// commands overflows the receive buffer while the paced commands arrive intact, and
// that the two H-bridge inputs are never driven together.
//
// Build and run (from the repository root):
//   g++ -std=c++17 -O2 -o /tmp/fw_test tests/firmware_host/firmware_host_test.cpp && /tmp/fw_test

#include "arduino_mock.h"

#include <cassert>
#include <sstream>
#include <vector>

MockSerial Serial;
unsigned long long g_micros = 0;
int g_pwm[32] = {0};
std::deque<char> g_wire;
unsigned long long g_wire_next_us = 0;
long g_rx_dropped = 0;

// ---------------- simulated plant ----------------
const double T_ROOM = 22.0;
const double CHI_H = 0.40, CHI_C = 0.14, TAU_TEC = 90.0;
double g_temp[6] = {T_ROOM, T_ROOM, T_ROOM, T_ROOM, T_ROOM, T_ROOM};
bool g_missing[6] = {false, false, false, false, false, false};
int g_bothActive = 0;

int harnessAnalogRead(int pin) {
  int ch = pin - 14;
  if (g_missing[ch]) return 1023;                      // open thermistor: rail
  double tK = g_temp[ch] + 273.15;
  double r = 100000.0 * std::exp(4540.0 * (1.0 / tK - 1.0 / 298.15));
  double adc = 1023.0 * r / (100000.0 + r);           // R = Rs adc / (1023 - adc)
  return (int)std::lround(adc);
}

void stepPlant(double dt) {
  int heat = g_pwm[9], cool = g_pwm[10];
  if (heat > 0 && cool > 0) g_bothActive++;
  double u = heat > 0 ? heat : -cool;
  double chi = u >= 0 ? CHI_H : CHI_C;
  g_temp[0] += dt / TAU_TEC * (T_ROOM + chi * u - g_temp[0]);
  // rod: T1 couples to the TEC, each node to the previous one, all lose to room
  double prev = g_temp[0];
  for (int i = 1; i < 6; i++) {
    double tc = (i == 1) ? 25.0 : 45.0;
    double next = g_temp[i] + dt * ((prev - g_temp[i]) / tc - (g_temp[i] - T_ROOM) / 900.0);
    prev = g_temp[i];
    g_temp[i] = next;
  }
}

#include "../../arduino/module_08_angstrom_controller/module_08_angstrom_controller.ino"

// ---------------- harness ----------------
struct Record {
  std::vector<std::string> f;
  double num(int i) const { return f[i] == "nan" ? NAN : std::atof(f[i].c_str()); }
};
const int N_FIELDS = 18;
enum { TIME, STATE, MODE, STEP, TA, TB, PERIOD, SETPOINT, FF, TEC, T1, T2, T3, T4, T5,
       CTRL, RATIO, SIGNED };

std::vector<Record> g_records;
std::vector<std::string> g_events, g_errors, g_other;
std::string g_header;
std::string g_pending;
bool g_heartbeat = true;
int g_failures = 0;

void check(const char *name, bool ok, const std::string &detail = "") {
  std::printf("  %s %s%s%s\n", ok ? "ok  " : "FAIL", name, detail.empty() ? "" : ": ",
              detail.c_str());
  if (!ok) g_failures++;
}

// Put one line on the wire to the Arduino (it arrives at 115200 baud).
void wireWrite(const std::string &line) {
  if (g_wire.empty()) g_wire_next_us = g_micros + CHAR_US;
  for (char c : line) g_wire.push_back(c);
  g_wire.push_back('\n');
}

// As the GUI: HB at once, every other command queued and sent one every 60 ms.
std::deque<std::string> g_outbox;
unsigned long long g_next_send_us = 0;
const unsigned long long SPACING_US = 60000;

void send(const std::string &cmd) {
  if (cmd == "HB") {
    wireWrite(cmd);
  } else {
    g_outbox.push_back(cmd);
  }
}

void pumpOutbox() {
  if (!g_outbox.empty() && g_micros >= g_next_send_us) {
    wireWrite(g_outbox.front());
    g_outbox.pop_front();
    g_next_send_us = g_micros + SPACING_US;
  }
}

void collectLines() {
  g_pending += Serial.out;
  Serial.out.clear();
  size_t nl;
  while ((nl = g_pending.find('\n')) != std::string::npos) {
    std::string line = g_pending.substr(0, nl);
    g_pending.erase(0, nl + 1);
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.rfind("DATA,time_ms", 0) == 0) {
      g_header = line;
    } else if (line.rfind("DATA,", 0) == 0) {
      Record r;
      std::stringstream ss(line.substr(5));
      std::string item;
      while (std::getline(ss, item, ',')) r.f.push_back(item);
      g_records.push_back(r);
      if (g_heartbeat) send("HB");
    } else if (line.rfind("EVENT,", 0) == 0) {
      g_events.push_back(line);
    } else if (line.rfind("ERROR", 0) == 0) {
      g_errors.push_back(line);
    } else {
      g_other.push_back(line);
    }
  }
}

// One loop() pass; the plant advances by the time the pass took (prints can block).
void stepOnce() {
  unsigned long long before = g_micros;
  g_micros += 20;                         // a fast loop() pass on the Uno
  wireTransfer();
  loop();
  stepPlant((g_micros - before) * 1e-6);
  collectLines();
  pumpOutbox();
}

// Run the sketch until n more DATA records have been printed.
void runRecords(size_t n) {
  size_t target = g_records.size() + n;
  while (g_records.size() < target) stepOnce();
}

// Run until every queued command has reached the sketch and been handled.
void flushCommands() {
  while (!g_outbox.empty() || !g_wire.empty() || !Serial.in.empty()) stepOnce();
  for (int i = 0; i < 200; i++) stepOnce();
}

// mean, amplitude, and phase of the fundamental at the given period, by least squares
void fitSine(const std::vector<double> &t, const std::vector<double> &y, double period,
             double &mean, double &amp, double &ph) {
  double w = 2 * M_PI / period;
  double S[3][3] = {{0}}, b[3] = {0};
  for (size_t i = 0; i < t.size(); i++) {
    double v[3] = {1, std::sin(w * t[i]), std::cos(w * t[i])};
    for (int r = 0; r < 3; r++) {
      b[r] += v[r] * y[i];
      for (int c = 0; c < 3; c++) S[r][c] += v[r] * v[c];
    }
  }
  // Gaussian elimination
  for (int k = 0; k < 3; k++)
    for (int r = k + 1; r < 3; r++) {
      double m = S[r][k] / S[k][k];
      for (int c = k; c < 3; c++) S[r][c] -= m * S[k][c];
      b[r] -= m * b[k];
    }
  double x[3];
  for (int r = 2; r >= 0; r--) {
    double s = b[r];
    for (int c = r + 1; c < 3; c++) s -= S[r][c] * x[c];
    x[r] = s / S[r][r];
  }
  mean = x[0];
  amp = std::hypot(x[1], x[2]);
  ph = std::atan2(x[2], x[1]);   // y = mean + amp sin(w t + ph)
}

int main() {
  setup();
  collectLines();
  runRecords(5);

  std::puts("Record format and idle state");
  check("header has the 18 Module 8 fields",
        g_header == "DATA,time_ms,run_state,control_mode,step,Ta_C,Tb_C,period_s,setpoint_C,"
                    "feedforward_pwm,TEC_C,T1_C,T2_C,T3_C,T4_C,T5_C,controller_pwm,"
                    "heating_to_cooling_ratio,signed_pwm");
  bool allFields = true;
  for (auto &r : g_records) allFields = allFields && r.f.size() == N_FIELDS;
  check("every record has 18 fields", allFields);
  const Record &r0 = g_records.back();
  check("starts OFF, MANUAL, step 0, PWM 0",
        r0.f[STATE] == "OFF" && r0.f[MODE] == "MANUAL" && r0.f[STEP] == "0" &&
            r0.f[SIGNED] == "0" && r0.f[PERIOD] == "nan" && r0.f[SETPOINT] == "nan");
  double dtRec = (g_records[4].num(TIME) - g_records[0].num(TIME)) / 4000.0;
  check("one record per second", std::fabs(dtRec - 1.0) < 0.02, std::to_string(dtRec) + " s");
  check("room temperature read back on all channels",
        std::fabs(r0.num(TEC) - T_ROOM) < 0.05 && std::fabs(r0.num(T5) - T_ROOM) < 0.05);

  std::puts("Rejected commands");
  size_t ev0 = g_events.size(), err0 = g_errors.size();
  send("SET PERIOD 50");        // below the minimum period
  send("SET MEAN 41");          // 41 + 16 above the 56 C setpoint limit
  send("SET KP abc");
  send("MODE FAST");
  send("JUMP");
  flushCommands();
  check("five bad commands give five ERROR lines and no event",
        g_errors.size() - err0 == 5 && g_events.size() == ev0);
  check("settings unchanged", periodS == 100.0f && meanC == 40.0f);

  std::puts("Zero-centered sinusoidal PWM (S18.3)");
  send("MODE PWM_SINE");
  send("SET PWMAMP 40");
  send("SET PWMOFF 0");
  send("SET PERIOD 100");
  send("START");
  flushCommands();
  size_t a = g_records.size();
  runRecords(200);
  int umin = 999, umax = -999;
  double usum = 0;
  for (size_t i = a + 100; i < g_records.size(); i++) {
    int u = (int)g_records[i].num(SIGNED);
    umin = std::min(umin, u); umax = std::max(umax, u); usum += u;
  }
  check("signed PWM spans -40 ... +40 with mean near 0",
        umin == -40 && umax == 40 && std::fabs(usum / 100) < 1.0,
        std::to_string(umin) + " / " + std::to_string(umax));
  check("Ta, Tb, setpoint are placeholders in PWM_SINE, period is shown",
        g_records.back().f[TA] == "nan" && g_records.back().f[SETPOINT] == "nan" &&
            g_records.back().num(PERIOD) == 100.0);
  check("run_state RUN and step advanced",
        g_records.back().f[STATE] == "RUN" && g_records.back().num(STEP) == 5);

  std::puts("Positive-only sinusoidal PWM (S18.4)");
  send("SET PWMOFF 40");
  flushCommands();
  a = g_records.size();
  runRecords(120);
  umin = 999;
  for (size_t i = a; i < g_records.size(); i++) umin = std::min(umin, (int)g_records[i].num(SIGNED));
  check("never cools with B = A", umin >= 0, std::to_string(umin));

  std::puts("Heartbeat timeout");
  g_heartbeat = false;
  runRecords(6);
  g_heartbeat = true;
  const Record &rt = g_records.back();
  check("run stops with TIMEOUT and zero output",
        rt.f[STATE] == "TIMEOUT" && rt.f[SIGNED] == "0" && g_pwm[9] == 0 && g_pwm[10] == 0);
  check("TIMEOUT event printed", g_events.back().find("TIMEOUT") != std::string::npos);

  std::puts("Feedforward with cooling scaling (S18.5-S18.6)");
  send("MODE FF");
  send("SET TROOM 22");
  send("SET CHIH 0.485");
  send("SET RATIO 2.82");
  send("SET MEAN 30");
  send("SET AMP 10");
  send("START");
  flushCommands();
  a = g_records.size();
  runRecords(150);
  bool ffOk = true, sawCooling = false;
  for (size_t i = a + 1; i < g_records.size(); i++) {
    const Record &r = g_records[i];
    double uff = (r.num(SETPOINT) - 22.0) / 0.485;
    double expect = uff < 0 ? uff * 2.82 : uff;
    if (std::fabs(r.num(FF) - uff) > 0.06 || std::fabs(r.num(SIGNED) - expect) > 0.6) ffOk = false;
    if (r.num(SIGNED) < 0) sawCooling = true;
    if (r.f[CTRL] != "nan" || r.num(RATIO) != 2.82) ffOk = false;
  }
  check("feedforward_pwm = (Tset - T_room)/chi_H, cooling scaled by r", ffOk && sawCooling);

  std::puts("PI tracking and phase continuity (S18.7, S18.10)");
  send("STOP");
  send("MODE PI");
  send("SET MEAN 35");
  send("SET AMP 5");
  send("SET PERIOD 100");
  send("SET KP 6");
  send("SET KI 0.04");
  send("START");
  flushCommands();
  runRecords(1200);
  std::vector<double> ts, tec, sp, t1, t2, t3;
  for (size_t i = g_records.size() - 500; i < g_records.size(); i++) {
    ts.push_back(g_records[i].num(TIME) / 1000.0);
    tec.push_back(g_records[i].num(TEC));
    sp.push_back(g_records[i].num(SETPOINT));
    t1.push_back(g_records[i].num(T1));
    t2.push_back(g_records[i].num(T2));
    t3.push_back(g_records[i].num(T3));
  }
  double mS, aS, pS, mT, aT, pT, m1, a1, p1, m2, a2, p2, m3, a3, p3;
  fitSine(ts, sp, 100, mS, aS, pS);
  fitSine(ts, tec, 100, mT, aT, pT);
  fitSine(ts, t1, 100, m1, a1, p1);
  fitSine(ts, t2, 100, m2, a2, p2);
  fitSine(ts, t3, 100, m3, a3, p3);
  // lag between neighbours, wrapped to (-pi, pi]; summed to unwrap
  auto wrap = [](double x) { return std::remainder(x, 2 * M_PI); };
  double lag12 = wrap(p1 - p2), lag23 = wrap(p2 - p3);
  char buf[200];
  std::snprintf(buf, sizeof buf, "setpoint %.2f +/- %.2f, TEC %.2f +/- %.2f, lag %.0f deg",
                mS, aS, mT, aT, (pS - pT) * 180 / M_PI);
  check("PI keeps the TEC mean near 35 C with a sinusoidal response",
        std::fabs(mT - 35.0) < 0.5 && aT > 1.0, buf);
  std::snprintf(buf, sizeof buf, "amplitudes T1/T2/T3 %.2f/%.2f/%.2f C, lags T1->T2 %.0f, T2->T3 %.0f deg",
                a1, a2, a3, lag12 * 180 / M_PI, lag23 * 180 / M_PI);
  check("rod oscillation shrinks and lags downstream",
        a2 < a1 && a3 < a2 && lag12 > 0 && lag23 > 0, buf);

  size_t before = g_records.size();
  send("SET PERIOD 200");
  flushCommands();
  runRecords(20);
  double maxJump = 0;
  for (size_t i = before - 10; i < g_records.size(); i++)
    maxJump = std::max(maxJump, std::fabs(g_records[i].num(SETPOINT) - g_records[i - 1].num(SETPOINT)));
  check("setpoint continuous across the 100 s -> 200 s change", maxJump < 0.35,
        "largest step " + std::to_string(maxJump) + " C");
  check("new period and step in the records",
        g_records.back().num(PERIOD) == 200.0 &&
            g_records.back().num(STEP) == g_records[before - 1].num(STEP) + 1);

  std::puts("Lowered limits (automatic shutdown test)");
  {
    send("STOP");
    flushCommands();
    for (int k = 0; k < 600; k++) runRecords(1);           // cool toward room
    size_t err0 = g_errors.size();
    send("SET T1LIMIT 57");                                   // above 56: refused
    send("SET TECLIMIT 61");                                  // above 60: refused
    flushCommands();
    bool refused = g_errors.size() == err0 + 2 && rodBaseLimitC == 56.0f && tecLimitC == 60.0f;
    double t1_now = g_records.back().num(T1);
    char cmd[40];
    std::snprintf(cmd, sizeof cmd, "SET T1LIMIT %.4f", t1_now + 2.0);
    send(cmd);
    send("MODE MANUAL");
    send("SET MANUAL 40");
    send("START");
    flushCommands();
    runRecords(2);                                            // START has taken effect
    size_t r0 = g_records.size();
    bool started = g_records.back().f[STATE] == "RUN";
    for (int k = 0; k < 600 && g_records.back().f[STATE] == "RUN"; k++) runRecords(1);
    const Record &rl = g_records.back();
    bool tripped_t1 = started && rl.f[STATE] == "SAFETY" && rl.num(T1) > t1_now + 2.0
                      && rl.num(T1) < t1_now + 2.3 && rl.num(TEC) < 60.0
                      && g_events.back().find("rod-base limit") != std::string::npos;
    check("SET T1LIMIT above 56 C or SET TECLIMIT above 60 C is refused", refused);
    check("lowered T1 limit (T1 + 2 C) stops a gentle heating run",
          tripped_t1, "T1 " + rl.f[T1] + " C after " + std::to_string(g_records.size() - r0) + " s");
    send("SET T1LIMIT 56");
    send("SET MANUAL 0");
    send("START");
    flushCommands();
    runRecords(2);
    check("restored to 56 C, START is accepted again",
          rodBaseLimitC == 56.0f && g_records.back().f[STATE] == "RUN");
  }

  std::puts("Safety limits");
  send("STOP");
  send("MODE MANUAL");
  send("SET MANUAL 255");
  send("START");
  flushCommands();
  size_t s0 = g_records.size();
  for (int k = 0; k < 2000 && g_records.back().f[STATE] == "RUN"; k++) runRecords(1);
  const Record &rs = g_records.back();
  check("heating stops with SAFETY at the first limit",
        rs.f[STATE] == "SAFETY" && rs.f[SIGNED] == "0" && g_pwm[9] == 0,
        "TEC " + rs.f[TEC] + " C, T1 " + rs.f[T1] + " C after " +
            std::to_string(g_records.size() - s0) + " s");
  check("SAFETY event names the limit",
        g_events.back().find("SAFETY SHUTDOWN") != std::string::npos);
  err0 = g_errors.size();
  send("START");
  flushCommands();
  check("START refused while over the limit",
        g_errors.size() == err0 + 1 && g_errors.back().find("START refused") != std::string::npos);
  send("SET MANUAL 0");
  send("STOP");
  flushCommands();
  runRecords(5);
  g_missing[1] = true;
  runRecords(2);
  err0 = g_errors.size();
  send("START");
  flushCommands();
  check("START refused with T1 unplugged (nan)",
        g_records.back().f[T1] == "nan" && g_errors.size() == err0 + 1);
  g_missing[1] = false;

  std::puts("Serial input");
  {
    // 14 settings at once, as a GUI without pacing would send them
    const char *burst[] = {"SET KP 6.0000", "SET KI 0.0400", "SET UMAX 255.0000",
                           "SET RATIO 2.8200", "SET CHIH 0.4850", "SET TROOM 22.0000",
                           "SET PWMOFF 0.0000", "SET PWMAMP 50.0000", "SET MANUAL 0.0000",
                           "SET FFADD 0.0000", "SET PERIOD 100.0000", "SET AMP 16.0000",
                           "SET MEAN 40.0000", "SET KP 6.0000"};
    long dropped0 = g_rx_dropped;
    size_t ev0 = g_events.size();
    for (const char *c : burst) wireWrite(c);
    flushCommands();
    long burst_dropped = g_rx_dropped - dropped0;
    size_t burst_events = g_events.size() - ev0;
    check("a burst of 14 settings overflows the 63-byte receive buffer",
          burst_dropped > 0 && burst_events < 14,
          std::to_string(burst_dropped) + " characters lost, " + std::to_string(burst_events) +
              " of 14 accepted");
    dropped0 = g_rx_dropped;
    ev0 = g_events.size();
    for (const char *c : burst) send(c);
    flushCommands();
    bool echoes = g_events.size() - ev0 == 14;
    for (size_t i = 0; echoes && i < 14; i++) {
      std::string want = std::string(burst[i]).substr(4);        // "KP 6.0000"
      echoes = g_events[ev0 + i].find("SET " + want) != std::string::npos;
    }
    check("the same settings one every 60 ms: none lost, every echo matches",
          g_rx_dropped == dropped0 && echoes,
          std::to_string(g_events.size() - ev0) + " accepted");
  }

  std::puts("Events and H-bridge");
  bool stepsOk = true;
  int last = 0;
  for (auto &e : g_events) {
    std::stringstream ss(e);
    std::string tag, tms, st;
    std::getline(ss, tag, ',');
    std::getline(ss, tms, ',');
    std::getline(ss, st, ',');
    int s = std::atoi(st.c_str());
    if (s != last + 1) stepsOk = false;
    last = s;
  }
  check("every EVENT increments step by one", stepsOk, std::to_string(g_events.size()) + " events");
  check("H-bridge inputs never driven together", g_bothActive == 0);

  std::printf("%s (%zu records, %zu events)\n", g_failures ? "FAIL" : "PASS", g_records.size(),
              g_events.size());
  return g_failures ? 1 : 0;
}
