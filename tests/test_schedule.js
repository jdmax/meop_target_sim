// Reduced-model checks of the discharge schedule in the JavaScript port, he3_model.js.
// Run from the repository root:  node tests/test_schedule.js
// tests/test_schedule.py runs the same checks on the Python model and compares the two ports.
"use strict";
const path = require("path");
const assert = require("assert");
const M = require(path.join(__dirname, "..", "he3_model.js"));

const H = 3600;
// In the MEOP-limited limit the cells share one P, with dP/dt = k (P_inf - P), k = alpha/T_b + beta and
// k P_inf = alpha M_ss/T_b (alpha = N_p/N, beta = (N_t/N)/T_1t). These levels give a chosen (k, P_inf).
const NP = 1, NT = 120, ALPHA = NP / (NP + NT), BETA = 0.1 / H, T1T = (NT / (NP + NT)) / BETA;
const level = (label, k_h, P) => { const k = k_h / H; return { label, M_ss: 100 * k * P / (k - BETA), T_b_s: ALPHA / (k - BETA), T1p_s: 1e9 }; };
const reducedTime = (k_h, P, P0, P1) => Math.log((P - P0) / (P - P1)) / k_h;
const SCHED = { sched_mode: "switch", single_level: 0, early_level: 1, late_level: 0, trigger: "auto", P_switch: 20, t_switch_h: 1,
  hyst_pct: 0, dwell_min: 0, goal_mode: "abs", goal_pct: 30, goal_frac: 90, on_h: 24 };
function run(levels, goal, xlim, over) {
  const G = NP * Math.max(...levels.map(l => 1 / l.T_b_s)) / xlim;
  return M.schedule(NP, NT, G, T1T, levels, Object.assign({}, SCHED, { goal_pct: 100 * goal }, over));
}

const CASES = [
  // [k_L, P_L], [k_H, P_H], goal, P*, low-only h, switched h (the hand calculations)
  [[0.3, 0.40], [0.5, 0.32], 0.35, 0.20, 6.93, 6.58],
  [[0.3, 0.40], [1.1, 0.364], 0.38, 0.35, 10.0, 6.0],
];
for (const [[kL, PL], [kH, PH], goal, Pstar, tLow, tSw] of CASES) {
  const levels = [level("Low", kL, PL), level("High", kH, PH)];
  const r = run(levels, goal, 1e-6);
  const Pred = (kH * PH - kL * PL) / (kH - kL);
  assert(Math.abs(r.sched.P_switch - Pred) < 1e-12);
  assert.strictEqual(+Pred.toFixed(2), Pstar);
  const exactLow = reducedTime(kL, PL, 0, goal), exactSw = reducedTime(kH, PH, 0, Pred) + reducedTime(kL, PL, Pred, goal);
  const low = r.levels[0].t_goal_s / H, sw = r.sched.t_goal_s / H;
  assert(Math.abs(low / exactLow - 1) < 1e-4, `${low} vs ${exactLow}`);
  assert(Math.abs(sw / exactSw - 1) < 1e-4, `${sw} vs ${exactSw}`);
  assert(Math.abs(low - tLow) < 0.05 && Math.abs(sw - tSw) < 0.05);
  const env = run(levels, goal, 1e-3, { sched_mode: "envelope" }), two = run(levels, goal, 1e-3, {});
  assert(Math.abs(env.sched.t_goal_s / two.sched.t_goal_s - 1) < 1e-12, "envelope = automatic switch");
  console.log(`reduced  P* ${Pred.toFixed(4)}  low ${low.toFixed(3)} h (exact ${exactLow.toFixed(3)})  switched ${sw.toFixed(3)} h (exact ${exactSw.toFixed(3)})`);
}
console.log("ok       JavaScript schedule checks");
