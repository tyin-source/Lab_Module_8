// S18.3 and S18.4 in the course simulator's one-lump mode (?s18=3), headless.
//
// For each case: Thermal model = One-lump TEC, Room / Initial T / tau / chi heating /
// chi cooling typed into the page's inputs, Sinusoidal PWM with u_a, u_b, period, then
// Apply / Reset; the page's own RK4 one-lump model is advanced and its TEC temperature
// recorded every 0.1 s, plus the page's analytical overlay where it exists.
//
// Usage (from the repository root): node tools/course_simulator_lumped.js data/module_08/course_simulator
const fs = require('fs');
const path = require('path');
const { chromium } = require('playwright');
const URL = 'https://sethfraden.github.io/Phys39F26-course/simulators/angstrom-six-sensor/?s18=3';

// [label, period, u_a, u_b, room, initial, tau, chi_h, chi_c, total s]
const CASES = [];
for (const P of [100, 200, 300]) for (const cc of [0.23, 0.46])
  CASES.push([`s18_3_P${P}_chic${cc}`, P, 0, 50, 22, 22, 80, 0.46, cc, 12 * Math.max(P, 100)]);
// S18.4 with our measured one-lump parameters (sim_parameters.csv: one_lump_*), room 21.73 C
const ONE = { chi_h: 0.4115, chi_c: 0.1732, tau: 59.1, room: 21.73 };
CASES.push(['s18_4_zero100', 100, 0, 50, ONE.room, ONE.room, ONE.tau, ONE.chi_h, ONE.chi_c, 2000]);
CASES.push(['s18_4_pos100', 100, 55, 50, ONE.room, ONE.room, ONE.tau, ONE.chi_h, ONE.chi_c, 2000]);

(async () => {
  const outDir = process.argv[2] || 'course_simulator_out';
  fs.mkdirSync(outDir, { recursive: true });
  const browser = await chromium.launch();
  const page = await browser.newPage();
  await page.goto(URL, { waitUntil: 'load' });
  const summary = [];
  for (const [label, P, ua, ub, room, init, tau, ch, cc, total] of CASES) {
    const res = await page.evaluate(({ P, ua, ub, room, init, tau, ch, cc, total }) => {
      const set = (id, v) => { const e = document.getElementById(id); e.value = String(v); };
      set('modelMode', 'lumped'); set('lumpRoom', room); set('lumpInitial', init);
      set('lumpTau', tau); set('lumpHeat', ch); set('lumpCool', cc);
      set('pwmMode', 'sine'); set('pwmOffset', ua); set('pwmAmplitude', ub); set('pwmPeriod', P);
      document.getElementById('btnReset').click();          // Apply / Reset
      sim.paused = true;
      const out = [];
      const n = Math.round(total / 0.1);
      for (let k = 1; k <= n; k++) {
        advance(0.1);
        out.push([sim.time, pwmAt(sim.time), sim.Tc, analyticalAt(sim.time)]);
      }
      return { out, lump: { ...lump }, drive: { ...drive } };
    }, { P, ua, ub, room, init, tau, ch, cc, total });
    const lines = res.out.map(r => [r[0].toFixed(2), r[1].toFixed(4), r[2].toFixed(6),
      Number.isFinite(r[3]) ? r[3].toFixed(6) : ''].join(','));
    fs.writeFileSync(path.join(outDir, `course_lump_${label}.csv`),
      'time_s,pwm,TEC_C,analytical_C\n' + lines.join('\n') + '\n');
    summary.push({ label, lump: res.lump, drive: res.drive });
    console.log(label, JSON.stringify(res.lump), `u_a ${res.drive.offset} u_b ${res.drive.amplitude} P ${res.drive.period} manual ${res.drive.manual} wave ${res.drive.wave}`);
  }
  fs.writeFileSync(path.join(outDir, 'course_lump_meta.json'), JSON.stringify({ url: URL, fetched: new Date().toISOString(), cases: summary }, null, 2) + '\n');
  await browser.close();
})();
