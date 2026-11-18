// Screenshots and the page's own "Save CSV" files from the course six-sensor simulator,
// set up through its controls exactly as a student would, for the S18.3, S18.4, S19,
// and Stage 6 simulations.
//
// For each case: open the page (one-lump cases with ?s18=3), set the controls, press
// Apply / Reset, pause the wall clock, advance the page's own model to the requested
// time, let it redraw, then screenshot the window and click Save CSV.
//
// Usage (from the repository root):
//   node tools/course_simulator_screens.js docs/figures/module_08/course_ui data/module_08/course_simulator/ui_csv
const fs = require('fs');
const path = require('path');
const { chromium } = require('playwright');
const BASE = 'https://sethfraden.github.io/Phys39F26-course/simulators/angstrom-six-sensor/';
const G = 2.55;                    // PWM counts per % duty
const ONE = { tau: 59.1, heat: 0.4115, cool: 0.1732, room: 21.73 };   // measured one-lump (S18.3/S18.4)

const LUMP = [];
for (const P of [100, 200, 300]) for (const cc of [0.23, 0.46])
  LUMP.push({ name: `s18_3_P${P}_chic${cc}`, P, ua: 0, ub: 50, room: 22, init: 22, tau: 80, heat: 0.46, cool: cc,
              t: 12 * P, residual: true });
LUMP.push({ name: 's18_3_P300_equal_startup', P: 300, ua: 0, ub: 50, room: 22, init: 22, tau: 80, heat: 0.46, cool: 0.46,
            t: 600, residual: false });
LUMP.push({ name: 's18_4_zero100', P: 100, ua: 0, ub: 50, room: ONE.room, init: ONE.room, tau: ONE.tau, heat: ONE.heat,
            cool: ONE.cool, t: 2000, residual: true });
LUMP.push({ name: 's18_4_pos100', P: 100, ua: 55, ub: 50, room: ONE.room, init: ONE.room, tau: ONE.tau, heat: ONE.heat,
            cool: ONE.cool, t: 2000, residual: true });

const COUPLED = [
  { name: 's19_p100', P: 100, kp: 6, ki: 0, t: 14400, views: ['thermal'] },
  { name: 's19_pi100', P: 100, kp: 6, ki: 0.04, t: 14400, views: ['thermal'] },
  { name: 's19_pi100_adj', P: 100, kp: 30, ki: 0.3, t: 14400, views: ['thermal'] },
  { name: 'stage6_pi300', P: 300, kp: 30, ki: 0.3, t: 18000, views: ['normal', 'thermal'] },
];

async function shoot(page, file) {
  await page.waitForTimeout(700);            // let requestAnimationFrame redraw
  await page.screenshot({ path: file });
}
async function saveCsv(page, file) {
  const [dl] = await Promise.all([page.waitForEvent('download'), page.click('#btnCSV')]);
  await dl.saveAs(file);
}

(async () => {
  const figDir = process.argv[2] || 'course_ui', csvDir = process.argv[3] || 'course_ui_csv';
  fs.mkdirSync(figDir, { recursive: true }); fs.mkdirSync(csvDir, { recursive: true });
  const browser = await chromium.launch();
  const ctx = await browser.newContext({ viewport: { width: 1680, height: 1050 }, deviceScaleFactor: 2,
                                         acceptDownloads: true });
  const page = await ctx.newPage();

  for (const c of LUMP) {
    await page.goto(BASE + '?s18=3', { waitUntil: 'load' });
    await page.evaluate((c) => {
      const set = (id, v) => { const e = document.getElementById(id); e.value = String(v); };
      set('modelMode', 'lumped'); set('lumpRoom', c.room); set('lumpInitial', c.init); set('lumpTau', c.tau);
      set('lumpHeat', c.heat); set('lumpCool', c.cool);
      set('pwmMode', 'sine'); set('pwmOffset', c.ua); set('pwmAmplitude', c.ub); set('pwmPeriod', c.P);
      document.getElementById('btnReset').click();
      sim.paused = true; document.getElementById('btnPause').textContent = 'Run';
      for (let k = 0; k < c.t; k++) advance(1.0);
      document.getElementById('showResidual').checked = c.residual;
    }, c);
    await shoot(page, path.join(figDir, `${c.name}.png`));
    await saveCsv(page, path.join(csvDir, `${c.name}.csv`));
    const readout = await page.textContent('#settledReadout');
    console.log(c.name, '|', readout.trim());
  }

  for (const c of COUPLED) {
    await page.goto(BASE, { waitUntil: 'load' });
    await page.evaluate((c) => {
      const set = (id, v, ev = 'change') => { const e = document.getElementById(id); e.value = String(v);
                                              e.dispatchEvent(new Event(ev)); };
      const chk = (id, on) => { const e = document.getElementById(id);
                                if (e.checked !== on) { e.checked = on; e.dispatchEvent(new Event('change')); } };
      set('modelMode', 'coupled'); set('material', 'aluminum');
      set('diam', 12.7, 'input'); set('len', 285, 'input');
      chk('manualChk', false);
      set('setpointMode', 'sine'); set('driveTa', 40); set('driveTb', 16); set('drivePeriod', c.P);
      chk('pOn', true); chk('iOn', c.ki > 0); chk('dOn', false);
      set('kp', +(c.kp / 2.55).toFixed(4)); set('ki', +(c.ki / 2.55).toFixed(4));
      document.getElementById('btnReset').click();
      sim.paused = true; document.getElementById('btnPause').textContent = 'Run';
      for (let k = 0; k < c.t; k++) advance(1.0);
    }, c);
    for (const v of c.views) {
      await page.click(v === 'normal' ? '#btnNormal' : '#btnThermal');
      await shoot(page, path.join(figDir, `${c.name}_${v === 'normal' ? 'visual' : 'thermal'}.png`));
    }
    await saveCsv(page, path.join(csvDir, `${c.name}.csv`));
    const st = await page.evaluate(() => ({ kp: pid.kp, ki: pid.ki, iOn: pid.iOn, L: sim.L, t: sim.time }));
    console.log(c.name, JSON.stringify(st));
  }
  await browser.close();
})();
