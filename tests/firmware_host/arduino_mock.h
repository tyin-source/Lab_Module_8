// Minimal stand-in for the Arduino core, so the Module 8 sketch can be compiled
// and exercised on a computer (tests/firmware_host/firmware_host_test.cpp).
// Only what the sketch uses is provided. Time is simulated: the test harness
// advances g_micros; analogRead() returns the ADC value of a temperature set by
// the harness's thermal model; analogWrite() records the PWM on each pin.
// The serial port behaves like the Uno's at 115200 baud: characters sent to the
// Arduino arrive one every 87 us into a 63-byte receive buffer (further characters
// are lost and counted), and printing blocks once the 64-byte transmit buffer is
// full, so a long print delays the sketch while input keeps arriving.
#pragma once

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <deque>
#include <math.h>
#include <algorithm>
#include <string>

using std::isnan;

typedef unsigned char byte;

// ---------------- pins and constants ----------------
const int A0 = 14, A1 = 15, A2 = 16, A3 = 17, A4 = 18, A5 = 19;
const int LOW = 0, HIGH = 1, OUTPUT = 1, INPUT = 0, DEC = 10;
#ifndef TWO_PI
#define TWO_PI 6.283185307179586476925286766559
#endif

// ---------------- simulated hardware (defined in the harness) ----------------
extern unsigned long long g_micros;
extern int g_pwm[32];
int harnessAnalogRead(int pin);

const unsigned long long CHAR_US = 87;      // one character at 115200 baud (10 bits)
const size_t RX_CAPACITY = 63;              // HardwareSerial ring buffer of 64
const size_t TX_CAPACITY = 64;
extern std::deque<char> g_wire;             // characters on their way to the Arduino
extern unsigned long long g_wire_next_us;   // arrival time of g_wire.front()
extern long g_rx_dropped;
void wireTransfer();

// Each call moves the clock by the Uno's 4 us resolution, so busy-wait loops such
// as the one at the end of setup() terminate.
inline unsigned long micros() { g_micros += 4; wireTransfer(); return (unsigned long)g_micros; }
inline unsigned long millis() { wireTransfer(); return (unsigned long)(g_micros / 1000ULL); }
inline void pinMode(int, int) {}
inline void digitalWrite(int pin, int value) { g_pwm[pin] = value ? 255 : 0; }
inline void analogWrite(int pin, int value) { g_pwm[pin] = value; }
inline int analogRead(int pin) { return harnessAnalogRead(pin); }
inline bool isDigit(char c) { return std::isdigit((unsigned char)c) != 0; }

// ---------------- F() strings ----------------
class __FlashStringHelper;
#define F(x) (reinterpret_cast<const __FlashStringHelper *>(x))

// ---------------- String ----------------
class String {
 public:
  std::string s;
  String() {}
  String(const char *c) : s(c) {}
  String(const std::string &x) : s(x) {}
  void reserve(unsigned int n) { s.reserve(n); }
  unsigned int length() const { return (unsigned int)s.size(); }
  char charAt(unsigned int i) const { return i < s.size() ? s[i] : 0; }
  void trim() {
    size_t a = s.find_first_not_of(" \t\r\n");
    size_t b = s.find_last_not_of(" \t\r\n");
    s = (a == std::string::npos) ? "" : s.substr(a, b - a + 1);
  }
  void toUpperCase() { for (auto &c : s) c = (char)std::toupper((unsigned char)c); }
  bool startsWith(const char *p) const { return s.rfind(p, 0) == 0; }
  String substring(unsigned int a) const { return a >= s.size() ? String() : String(s.substr(a)); }
  String substring(unsigned int a, unsigned int b) const { return String(s.substr(a, b - a)); }
  int indexOf(char c, unsigned int from = 0) const {
    size_t p = s.find(c, from);
    return p == std::string::npos ? -1 : (int)p;
  }
  float toFloat() const { return (float)std::atof(s.c_str()); }
  long toInt() const { return std::atol(s.c_str()); }
  String &operator+=(char c) { s += c; return *this; }
  bool operator==(const char *o) const { return s == o; }
  bool operator==(const String &o) const { return s == o.s; }
  bool operator!=(const char *o) const { return s != o; }
  const char *c_str() const { return s.c_str(); }
};

// ---------------- Serial ----------------
class MockSerial {
 public:
  std::string out;               // everything printed, read by the harness
  std::deque<char> in;           // bytes waiting to be read by the sketch
  unsigned long long tx_done_us = 0;   // when the transmit buffer will be empty

  void emit(const std::string &x) {
    for (char c : x) {
      // block while the transmit buffer is full
      if (tx_done_us > g_micros + TX_CAPACITY * CHAR_US) {
        g_micros = tx_done_us - TX_CAPACITY * CHAR_US;
        wireTransfer();
      }
      tx_done_us = std::max(tx_done_us, g_micros) + CHAR_US;
      out += c;
    }
  }

  void begin(long) {}
  int available() { return (int)in.size(); }
  int read() {
    if (in.empty()) return -1;
    char c = in.front();
    in.pop_front();
    return (unsigned char)c;
  }

  void print(const char *x) { emit(x); }
  void print(const __FlashStringHelper *x) { emit(reinterpret_cast<const char *>(x)); }
  void print(const String &x) { emit(x.s); }
  void print(char c) { emit(std::string(1, c)); }
  void print(int x, int = DEC) { emit(std::to_string(x)); }
  void print(long x, int = DEC) { emit(std::to_string(x)); }
  void print(unsigned long x, int = DEC) { emit(std::to_string(x)); }
  void print(double x, int digits = 2) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.*f", digits, x);
    emit(buf);
  }
  template <typename T> void println(T x) { print(x); emit("\r\n"); }
  template <typename T> void println(T x, int d) { print(x, d); emit("\r\n"); }
  void println() { emit("\r\n"); }
};

extern MockSerial Serial;

// Move the characters whose arrival time has come into the receive buffer.
inline void wireTransfer() {
  while (!g_wire.empty() && g_wire_next_us <= g_micros) {
    char c = g_wire.front();
    g_wire.pop_front();
    if (Serial.in.size() < RX_CAPACITY) {
      Serial.in.push_back(c);
    } else {
      g_rx_dropped++;
    }
    g_wire_next_us += CHAR_US;
  }
}
