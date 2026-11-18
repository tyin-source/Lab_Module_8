// Drive the course's six-sensor TEC-rod simulator headlessly and save CSV records.
//
// Loads https://sethfraden.github.io/Phys39F26-course/simulators/angstrom-six-sensor/
// in headless Chromium, sets the controls through the page's own inputs (material,
// diameter, length, Angstrom sinusoid setpoint, Kp, Ki, P/I/D switches), clicks
// Apply / Reset, pauses the wall-clock loop, and then advances the page's own
// advance(dt) in 1 s steps, recording one row per second:
//   time_s, setpoint_C, u_percent, pwm (= 2.55 u), Tplate_C (TEC rod-side plate),
//   T0_C (the simulator's PID sensor at x = 0, with its 2 s lag), and the rod temperature
//   at our five sensor positions (TEC face + hole 1, 3, 5, 6, 7).
// The model code is the page's; nothing in it is changed.
//
// Usage (from the repository root; needs Node and Playwright with Chromium):
//   node tools/course_simulator_driver.js data/module_08/course_simulator
const fs = require('fs');
const path = require('path');
const crypto = require('crypto');
const { chromium } = require('playwright');

const URL = 'https://sethfraden.github.io/Phys39F26-course/simulators/angstrom-six-sensor/';
const LEN_MM = 285;                       // slider step 5 mm; our rod is 286.6 mm
const X1_MM = 17.6;                       // TEC face to hole 1 (28.66 - 26.90 cm)
const POS_MM = [0, 24, 50, 76, 102].map(x => X1_MM + x);   // T1-T5 (holes 1, 3, 5, 6, 7)
const PWM_PER_PERCENT = 2.55;

// [label, period_s, Kp (PWM/C), Ki (PWM/(C s)), total simulated time s]
const RUNS = [
  ['s19_p100',  100, 6, 0, 14400],
  ['s19_pi100', 100, 6, 0.04, 14400],
  ['s19_pi100_adj', 100, 30, 0.3, 14400],
  ['pi200', 200, 30, 0.3, 14400],
  ['pi300', 300, 30, 0.3, 18000],
  ['pi400', 400, 30, 0.3, 18000],
  ['pi500', 500, 30, 0.3, 18000],
];

(async () => {
  const outDir = process.argv[2] || 'course_simulator_out';
  fs.mkdirSync(outDir, { recursive: true });
  const browser = await chromium.launch();
  const page = await browser.newPage();
  const resp = await page.goto(URL, { waitUntil: 'load' });
  const html = await resp.text();
  const sha = crypto.createHash('sha256').update(html).digest('hex');
  const meta = { url: URL, fetched: new Date().toISOString(), sha256: sha, length_mm: LEN_MM,
                 diameter_mm: 12.7, positions_from_tec_face_mm: POS_MM,
                 gain_conversion: 'Kp[%/C] = Kp[PWM/C] / 2.55, Ki likewise', runs: [] };
  for (const [label, period, kpPwm, kiPwm, total] of RUNS) {
    const kp = kpPwm / PWM_PER_PERCENT, ki = kiPwm / PWM_PER_PERCENT;
    const rows = await page.evaluate(({ period, kp, ki, total, LEN_MM, POS_MM }) => {
      const set = (id, v, ev = 'change') => { const e = document.getElementById(id);
        e.value = String(v); e.dispatchEvent(new Event(ev)); };
      const chk = (id, on) => { const e = document.getElementById(id);
        if (e.checked !== on) { e.checked = on; e.dispatchEvent(new Event('change')); } };
      set('modelMode', 'coupled'); set('material', 'aluminum');
      set('diam', 12.7, 'input'); set('len', LEN_MM, 'input');
      chk('manualChk', false);
      set('setpointMode', 'sine');
      set('driveTa', 40); set('driveTb', 16); set('drivePeriod', period);
      chk('pOn', true); chk('iOn', ki > 0); chk('dOn', false);
      set('kp', kp); set('ki', ki);
      sim.paused = true;
      resetSim();                          // Apply / Reset
      const out = [];
      const fr = POS_MM.map(x => x / (sim.L * 1000));
      for (let t = 1; t <= total; t++) {
        advance(1.0);
        out.push([sim.time, pid.sp, pid.u, sim.Tc, sim.Ts[0], ...fr.map(f => rodTempAt(f))]);
      }
      return { out, L: sim.L, d: sim.d, kp: pid.kp, ki: pid.ki, iOn: pid.iOn, dOn: pid.dOn,
               Ta: setpointDrive.Ta, Tb: setpointDrive.Tb, P: setpointDrive.period,
               mat: sim.matKey, mode: lump.mode, manual: drive.manual };
    }, { period, kp, ki, total, LEN_MM, POS_MM });
    const header = 'time_s,setpoint_C,u_percent,pwm,Tplate_C,T0_C,T1_C,T2_C,T3_C,T4_C,T5_C';
    const lines = rows.out.map(r => [r[0].toFixed(2), r[1].toFixed(4), r[2].toFixed(3),
      (r[2] * 2.55).toFixed(2), ...r.slice(3).map(v => v.toFixed(4))].join(','));
    fs.writeFileSync(path.join(outDir, `course_sim_${label}.csv`), header + '\n' + lines.join('\n') + '\n');
    meta.runs.push({ label, period_s: period, Kp_PWM_per_C: kpPwm, Ki_PWM_per_Cs: kiPwm,
      Kp_percent_per_C: rows.kp, Ki_percent_per_Cs: rows.ki, I_on: rows.iOn, D_on: rows.dOn,
      Ta_C: rows.Ta, Tb_C: rows.Tb, model: rows.mode, material: rows.mat, manual: rows.manual,
      length_m: rows.L, diameter_m: rows.d, simulated_s: total });
    console.log(`${label}: ${rows.out.length} s, Kp ${rows.kp.toFixed(3)} %/C, Ki ${rows.ki.toFixed(4)} %/(C s), P ${rows.P} s, L ${rows.L} m`);
  }
  fs.writeFileSync(path.join(outDir, 'course_sim_meta.json'), JSON.stringify(meta, null, 2) + '\n');
  await browser.close();
})();
