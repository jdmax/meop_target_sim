// he3_model.js -- line-for-line JavaScript port of he3_model.py.
// Used for the instant first render and as a fallback when Python can't load.
(function (root) {
  "use strict";
  const KB = 1.380649e-23, NA = 6.02214076e23, M3 = 3.01603e-3, M_ATOM = M3 / NA;
  const RG = 8.314462, CP = 2.5 * RG / M3, GRAV = 9.81;
  const ETA_300 = 17.2e-6, ETA_EXP = 0.68, RHO_D_OVER_ETA = 1.32;
  const eta = T => ETA_300 * Math.pow(T / 300, ETA_EXP);
  const k_gas = T => 3.75 * RG / M3 * eta(T);
  const n_num = (T, p) => p / (KB * T);
  const rho = (T, p) => p * M_ATOM / (KB * T);
  const nD = T => RHO_D_OVER_ETA * eta(T) / M_ATOM;

  const GLASS_T = [4, 6, 10, 20, 30, 50, 77, 100, 150, 200, 250, 300];
  const GLASS_K = [0.10, 0.11, 0.12, 0.15, 0.20, 0.32, 0.45, 0.55, 0.75, 0.92, 1.03, 1.12];
  function k_glass(T) {
    if (T <= GLASS_T[0]) return GLASS_K[0];
    if (T >= GLASS_T[GLASS_T.length - 1]) return GLASS_K[GLASS_K.length - 1];
    const x = Math.log(T);
    for (let i = 1; i < GLASS_T.length; i++) {
      const x1 = Math.log(GLASS_T[i]);
      if (x <= x1) {
        const x0 = Math.log(GLASS_T[i - 1]), f = (x - x0) / (x1 - x0);
        return Math.exp(Math.log(GLASS_K[i - 1]) + f * (Math.log(GLASS_K[i]) - Math.log(GLASS_K[i - 1])));
      }
    }
    return GLASS_K[GLASS_K.length - 1];
  }

  function trap(y, x) { let s = 0; for (let i = 1; i < x.length; i++) s += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]); return s; }
  function cumtrap(y, x) { const o = [0]; let s = 0; for (let i = 1; i < x.length; i++) { s += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]); o.push(s); } return o; }

  const GRID_N = 120;
  function t_grid(lo, hi) {
    const span = hi - lo, pts = new Set([lo, hi]);
    for (let i = 0; i < GRID_N; i++) {
      const f = Math.pow(10, -9 + 9 * i / (GRID_N - 1));
      pts.add(lo + f * span); pts.add(hi - f * span);
      pts.add(lo * Math.pow(hi / lo, i / (GRID_N - 1)));
    }
    return [...pts].filter(t => t >= lo && t <= hi).sort((a, b) => a - b);
  }
  const cache = new Map();
  function grid_K(lo, hi, Aw, Ag) {
    const key = lo.toFixed(9) + "|" + hi.toFixed(9) + "|" + Aw + "|" + Ag;
    let g = cache.get(key);
    if (!g) {
      const T = t_grid(lo, hi), K = T.map(t => Aw * k_glass(t) + Ag * k_gas(t));
      g = [T, K, trap(K, T)];
      if (cache.size > 200) cache.clear();
      cache.set(key, g);
    }
    return g;
  }
  function Tube(d, wall) { const r = d / 2; return { d, A_gas: Math.PI * r * r, A_wall: Math.PI * ((r + wall) ** 2 - r * r) }; }

  function solve_free(L, lo, hi, mdot, tube) {
    const [T, K, IK] = grid_K(lo, hi, tube.A_wall, tube.A_gas);
    const a = mdot * CP, n = T.length;
    if (Math.abs(a) < 1e-30) { const F = -IK / L; return [T, K.map(k => k / (-F)), F, false]; }
    function LdL(F) {
      let s = 0, s2 = 0, g0 = K[0] / (a * T[0] - F);
      for (let i = 1; i < n; i++) {
        const g1 = K[i] / (a * T[i] - F), dx = 0.5 * (T[i] - T[i - 1]);
        s += (g0 + g1) * dx; s2 += (g0 * g0 / K[i - 1] + g1 * g1 / K[i]) * dx; g0 = g1;
      }
      return [s, s2];
    }
    const Fmax = Math.min(a * lo, a * hi), scale = IK / L;
    const Fhi = Fmax - 1e-12 * Math.max(scale, Math.abs(Fmax));
    let F, clamped = false;
    if (LdL(Fhi)[0] < L) { F = Fhi; clamped = true; }
    else {
      let Flo = Fmax - scale;
      while (LdL(Flo)[0] > L) Flo = Fmax - 4 * (Fmax - Flo);
      let lb = Flo, hb = Fhi; F = Flo;
      for (let it = 0; it < 100; it++) {
        const [v, dv] = LdL(F), r = v - L;
        if (Math.abs(r) < 1e-10 * L) break;
        if (r < 0) lb = F; else hb = F;
        let Fn = dv > 0 ? F - r / dv : 0.5 * (lb + hb);
        if (!(lb < Fn && Fn < hb)) Fn = 0.5 * (lb + hb);
        F = Fn;
      }
    }
    const w = T.map((t, i) => K[i] / (a * t - F)), norm = L / trap(w, T);
    return [T, w.map(x => x * norm), F, clamped];
  }

  function free_segment(L, T0, T1, mdot_s, tube) {
    if (Math.abs(T1 - T0) < 1e-6) return { iso: true, T: T0, F: mdot_s * CP * T0, clamped: false };
    if (T1 > T0) { const [T, w, F, cl] = solve_free(L, T0, T1, mdot_s, tube); return { iso: false, T, w, s: cumtrap(w, T), F, clamped: cl }; }
    const [T, w, F, cl] = solve_free(L, T1, T0, -mdot_s, tube);
    return { iso: false, T, w, s: cumtrap(w, T).map(x => L - x), F: -F, clamped: cl };
  }

  function solve_leg(leg, tube, p, mdot, want) {
    const nodes = leg.nodes, frees = leg.free, a = mdot * CP;
    let I_rho = 0, I_nu = 0, I_rd = 0, I_n = 0, s0 = 0, clamped = false;
    const prof = [], Fs = [];
    nodes.forEach(([Tn, cl], i) => {
      if (cl > 0) { I_rho += rho(Tn, p) * cl; I_nu += eta(Tn) / rho(Tn, p) * cl; I_rd += cl / nD(Tn); I_n += n_num(Tn, p) * cl; }
      if (want) { prof.push([s0, Tn]); prof.push([s0 + cl, Tn]); }
      s0 += cl;
      if (i < frees.length) {
        const L = frees[i], seg = free_segment(L, Tn, nodes[i + 1][0], mdot, tube);
        if (seg.iso) {
          const t = seg.T;
          I_rho += rho(t, p) * L; I_nu += eta(t) / rho(t, p) * L; I_rd += L / nD(t); I_n += n_num(t, p) * L;
        } else {
          const T = seg.T, w = seg.w;
          I_rho += trap(T.map((t, j) => rho(t, p) * w[j]), T);
          I_nu += trap(T.map((t, j) => eta(t) / rho(t, p) * w[j]), T);
          I_rd += trap(T.map((t, j) => w[j] / nD(t)), T);
          I_n += trap(T.map((t, j) => n_num(t, p) * w[j]), T);
          clamped = clamped || seg.clamped;
          if (want) {
            const pts = seg.s.map((s, j) => [s, T[j]]).sort((u, v) => u[0] - v[0] || u[1] - v[1]);
            const step = Math.max(1, Math.floor(pts.length / 60));
            for (let j = 0; j < pts.length; j += step) prof.push([s0 + pts[j][0], pts[j][1]]);
            prof.push([s0 + pts[pts.length - 1][0], pts[pts.length - 1][1]]);
          }
        }
        Fs.push(seg.F); s0 += L;
      }
    });
    const heats = nodes.map(([Tn, , kind], i) => {
      const Fin = i === 0 ? a * Tn : Fs[i - 1], Fout = i === nodes.length - 1 ? a * Tn : Fs[i];
      return [kind, Tn, Fin - Fout];
    });
    const out = { I_rho, R_fric: 128 / (Math.PI * tube.d ** 4) * I_nu, R_diff: I_rd / tube.A_gas, N_tube: tube.A_gas * I_n, heats, length: s0, clamped };
    if (want) out.profile = prof.sort((u, v) => u[0] - v[0] || u[1] - v[1]).map(([s, t]) => [+(s * 100).toFixed(4), +t.toFixed(3)]);
    return out;
  }

  function thermosiphon(A_, B_, tube, p) {
    function resid(md) {
      const A = solve_leg(A_, tube, p, md), B = solve_leg(B_, tube, p, -md);
      const drive = GRAV * (B.I_rho - A.I_rho);
      return [md * (A.R_fric + B.R_fric) - drive, drive, A.R_fric + B.R_fric];
    }
    const [r0, d0, R0] = resid(0);
    if (Math.abs(d0) < 1e-15) return 0;
    const sgn = d0 > 0 ? 1 : -1;
    let hi = sgn * Math.abs(d0) / R0, rhi = resid(hi)[0];
    while (rhi * sgn < 0) { hi *= 2; rhi = resid(hi)[0]; }
    let lo = 0, rlo = r0, md = hi;
    for (let it = 0; it < 80; it++) {
      md = (lo * rhi - hi * rlo) / (rhi - rlo);
      const r = resid(md)[0];
      if (Math.abs(r) < 1e-9 * Math.abs(d0) || Math.abs(hi - lo) < 1e-9 * Math.abs(md)) break;
      if ((r > 0) === (rhi > 0)) { hi = md; rhi = r; rlo *= 0.5; } else { lo = md; rlo = r; rhi *= 0.5; }
    }
    return md;
  }

  function build_legs(P) {
    const cm = 0.01, Tc = P.Tc, Th = P.Th, H = P.length_cm * cm;
    const useInt = P.intercept && (P.mode === "single" || P.mode === "loop");
    const Ti = P.T_int, si = P.s_int_cm * cm, errs = [];
    function cold_leg(clamp) {
      if (useInt) {
        if (!(clamp + 0.2 * cm < si && si < H - 0.2 * cm)) {
          errs.push("Intercept height must sit between the cold clamp and the pumping cell.");
          return { nodes: [[Tc, clamp, "cold"], [Th, 0, "warm"]], free: [H - clamp] };
        }
        return { nodes: [[Tc, clamp, "cold"], [Ti, 0, "intercept"], [Th, 0, "warm"]], free: [si - clamp, H - si] };
      }
      return { nodes: [[Tc, clamp, "cold"], [Th, 0, "warm"]], free: [H - clamp] };
    }
    if (P.mode === "single") return [[cold_leg(0)], errs];
    let he = P.engine_cm * cm; const gap = P.gap_cm * cm, Te = P.T_engine;
    if (H - gap - he < 0.5 * cm) {
      errs.push("Legs are too short for the engine section; lengthen them or shorten the engine.");
      he = Math.max(0, H - gap - 0.5 * cm);
    }
    const riser = { nodes: [[Tc, 0, "cold"], [Te, he, "engine"], [Th, 0, "warm"]], free: [gap, H - gap - he] };
    return [[riser, cold_leg(gap + he)], errs];
  }

  function expm2(A, t) {
    const a = A[0][0], b = A[0][1], c = A[1][0], d = A[1][1], m = 0.5 * (a + d), det = a * d - b * c, q2 = m * m - det;
    const q = q2 > 0 ? Math.sqrt(q2) : 0;
    if (q * t < 1e-8) { const e = Math.exp(m * t); return [[e * (1 + (a - m) * t), e * b * t], [e * c * t, e * (1 + (d - m) * t)]]; }
    const l1 = m + q, l2 = m - q, e1 = Math.exp(l1 * t), e2 = Math.exp(l2 * t), k = 1 / (l1 - l2);
    return [[k * (e1 * (a - l2) - e2 * (a - l1)), k * (e1 - e2) * b], [k * (e1 - e2) * c, k * (e1 * (d - l2) - e2 * (d - l1))]];
  }
  function ss(A, b) {
    const det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
    return [-(A[1][1] * b[0] - A[0][1] * b[1]) / det, -(-A[1][0] * b[0] + A[0][0] * b[1]) / det];
  }
  function evolve(rates, P0, t_on, t_end, npts) {
    const t = []; for (let i = 0; i < npts; i++) t.push(t_end * i / (npts - 1));
    if (t_on > 0 && t_on < t_end && !t.includes(t_on)) { t.push(t_on); t.sort((a, b) => a - b); }
    const [A_on, b_on] = rates.on, [A_off, b_off] = rates.off, ss_on = ss(A_on, b_on), ss_off = ss(A_off, b_off);
    const op = [], ot = []; let xsw = null;
    for (const ti of t) {
      let A, xs, x0, dt;
      if (ti <= t_on) { A = A_on; xs = ss_on; x0 = P0; dt = ti; }
      else {
        if (!xsw) {
          const E = expm2(A_on, t_on), d0 = [P0[0] - ss_on[0], P0[1] - ss_on[1]];
          xsw = [ss_on[0] + E[0][0] * d0[0] + E[0][1] * d0[1], ss_on[1] + E[1][0] * d0[0] + E[1][1] * d0[1]];
        }
        A = A_off; xs = ss_off; x0 = xsw; dt = ti - t_on;
      }
      const E = expm2(A, dt), d0 = [x0[0] - xs[0], x0[1] - xs[1]];
      op.push(xs[0] + E[0][0] * d0[0] + E[0][1] * d0[1]); ot.push(xs[1] + E[1][0] * d0[0] + E[1][1] * d0[1]);
    }
    return [t, op, ot, ss_on];
  }
  function eig_times(A) {
    const a = A[0][0], b = A[0][1], c = A[1][0], d = A[1][1], m = 0.5 * (a + d);
    const q = Math.sqrt(Math.max(m * m - (a * d - b * c), 0));
    return [-1 / (m + q), -1 / (m - q)];
  }

  function compute(P) {
    const p = P.p_mbar * 100, Tc = P.Tc, Th = P.Th, tube = Tube(P.id_mm * 1e-3, P.wall_mm * 1e-3);
    const [legs, errs] = build_legs(P);
    const mdot = P.mode === "single" ? 0 : thermosiphon(legs[0], legs[1], tube, p);
    const sol = legs.map((leg, i) => solve_leg(leg, tube, p, i === 0 ? mdot : -mdot, true));
    const G_diff = sol.reduce((s, x) => s + 1 / x.R_diff, 0), G_conv = Math.abs(mdot) / M_ATOM, G = G_diff + G_conv;
    const Np = n_num(Th, p) * P.V_pump_cc * 1e-6, Nt = n_num(Tc, p) * P.V_target_cc * 1e-6;
    const heat = { cold: 0, intercept: 0, engine: 0, warm: 0 };
    sol.forEach(s => s.heats.forEach(([k, , q]) => { heat[k] += q; }));
    const g_op = 1 / P.T_op_s, g_rp = 1 / P.T1p_on_s, g_rp_off = 1 / (P.T1p_off_h * 3600), g_t = 1 / (P.T1t_h * 3600);
    const kp = G / Np, kt = G / Nt, Plim = P.P_lim / 100;
    const mat = (gop, grp) => [[[-(gop + grp + kp), kp], [kt, -(kt + g_t)]], [gop * Plim, 0]];
    const rates = { on: mat(g_op, g_rp), off: mat(0, g_rp_off) };
    const t_end = P.run_h * 3600, t_on = Math.min(P.on_h, P.run_h) * 3600;
    const [t, pp, pt, ss_on] = evolve(rates, [0, 0], t_on, t_end, 361);
    let t90 = null; const target = 0.9 * ss_on[1];
    const pt_at = tt => evolve(rates, [0, 0], 1e30, tt, 2)[2][1];
    if (ss_on[1] > 0) {
      let lo = 0, hi = 1;
      while (pt_at(hi) < target && hi < 1e9) hi *= 2;
      if (hi < 1e9) { for (let i = 0; i < 60; i++) { const mid = 0.5 * (lo + hi); if (pt_at(mid) < target) lo = mid; else hi = mid; } t90 = hi; }
    }
    const warnings = errs.slice(), Re = 4 * Math.abs(mdot) / (Math.PI * tube.d * eta(Tc));
    if (Re > 50) warnings.push(`Reynolds number ${Re.toFixed(0)} at the cold end: gas may lag the wall temperature, so the 1-D model overestimates the flow.`);
    if (heat.cold > 2.0) warnings.push("Cold-stage load exceeds the PT420 rating of about 2 W at 4.2 K.");
    if (sol.some(s => s.clamped)) warnings.push("Flow is strong enough to flatten a free section's temperature profile; treat these results as rough.");
    return {
      ok: true, engine: "javascript", mdot_mg_h: mdot * 3.6e9, Re, G, G_diff, G_conv, Np, Nt,
      N_tubes: sol.reduce((s, x) => s + x.N_tube, 0), tau_x_target_s: Nt / G, tau_x_pump_s: Np / G,
      heat_mW: Object.fromEntries(Object.entries(heat).map(([k, v]) => [k, v * 1e3])),
      legs: sol.map((s, i) => ({ profile: s.profile, length_cm: s.length * 100, nodes: legs[i].nodes.map(n => [n[0], n[1] * 100, n[2]]), free_cm: legs[i].free.map(f => f * 100), G_diff: 1 / s.R_diff })),
      ss_pump: ss_on[0], ss_target: ss_on[1], times_s: eig_times(rates.on[0]), t90_target_s: t90,
      sealed_Mss: g_op * Plim / (g_op + g_rp), sealed_Tb_s: 1 / (g_op + g_rp), production_atoms_s: g_op * Plim * Np,
      t_h: t.map(x => x / 3600), P_pump: pp, P_target: pt, t_on_h: t_on / 3600, warnings,
    };
  }
  function run(json) { try { return JSON.stringify(compute(JSON.parse(json))); } catch (e) { return JSON.stringify({ ok: false, error: String(e) }); } }
  root.He3ModelJS = { run, compute };
  if (typeof module !== "undefined") module.exports = { run, compute };
})(typeof window !== "undefined" ? window : globalThis);
