// he3_model.js -- line-for-line JavaScript port of he3_model.py.
// Used for the instant first render and as a fallback when Python can't load.
(function (root) {
  "use strict";
  const KB = 1.380649e-23, NA = 6.02214076e23, M3 = 3.01603e-3, M_ATOM = M3 / NA;
  const RG = 8.314462, CP = 2.5 * RG / M3, GRAV = 9.81, SIGMA = 5.670374419e-8;
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

  function two_cell(Np, Nt, G, tau_p, tau_t) {
    const N = Np + Nt, t_ex = Np * Nt / (N * G);
    const A = [[-(1 / tau_p + Nt / N / t_ex), Nt / N / t_ex], [Np / N / t_ex, -(1 / tau_t + Np / N / t_ex)]];
    const dg = 1 / tau_p - 1 / tau_t, mean = 0.5 * (1 / tau_p + 1 / tau_t + 1 / t_ex);
    const root = 0.5 * Math.sqrt(1 / t_ex ** 2 + dg * dg + 2 / t_ex * (Nt - Np) / N * dg);
    return { N, Np_N: Np / N, Nt_N: Nt / N, t_ex, tau_p, tau_t, tau_s: 1 / (mean + root), tau_l: 1 / (mean - root), A };
  }
  function amplitudes(sys, Mp0, Mt0) {
    const A = sys.A, rs = 1 / sys.tau_s, rl = 1 / sys.tau_l;
    const dMp = A[0][0] * Mp0 + A[0][1] * Mt0, dMt = A[1][0] * Mp0 + A[1][1] * Mt0;
    const a_l = (dMp + rs * Mp0) / (rs - rl), a_lt = (dMt + rs * Mt0) / (rs - rl);
    return { a_s: Mp0 - a_l, a_l, a_st: a_lt - Mt0, a_lt };
  }
  function relax_curves(sys, amp, t) {
    const rs = 1 / sys.tau_s, rl = 1 / sys.tau_l;
    return [t.map(x => amp.a_s * Math.exp(-rs * x) + amp.a_l * Math.exp(-rl * x)),
            t.map(x => amp.a_lt * Math.exp(-rl * x) - amp.a_st * Math.exp(-rs * x))];
  }

  function disk_view(r, L) {
    const R = r / L, X = 2 + 1 / (R * R);
    return 2 / (X + Math.sqrt(X * X - 4));
  }
  function stage_loads(P, heat, tube, n_legs) {
    const Th = P.Th, Tc = P.Tc, T1 = P.T_shield, q_int = heat.intercept, q_eng = heat.engine;
    const rad1 = SIGMA * P.eps_shield * P.A_shield_m2 * (Th ** 4 - T1 ** 4);
    const rad2 = SIGMA * P.eps_cold * P.A_cold_cm2 * 1e-4 * (T1 ** 4 - Tc ** 4);
    const bore = n_legs * SIGMA * tube.A_gas * disk_view(tube.d / 2, P.length_cm * 0.01) * (Th ** 4 - Tc ** 4);
    const first = { T: T1, intercept: Math.max(0, q_int), radiation: rad1, other: P.other1_W };
    first.total = first.intercept + rad1 + first.other;
    const second = { T: Tc, conduction: heat.cold, engine: Math.max(0, q_eng), radiation: rad2, bore, other: P.other2_mW * 1e-3 };
    second.total = second.conduction + second.engine + rad2 + bore + second.other;
    return { first, second, intercept_heater: Math.max(0, -q_int), engine_heater: Math.max(0, -q_eng), warm_supply: -heat.warm };
  }

  // PT420-RM reference capacity curve (Cryomech, 50 and 60 Hz), read off the published chart.
  const PT420_Q1 = [0, 25, 50, 75, 100], PT420_Q2 = [0, 0.5, 1.8, 5, 8, 10];
  const PT420_T1 = [[28.6, 33.5, 40.2, 49.4, 62.7], [29.0, 34.2, 40.9, 50.3, 65.5], [29.7, 36.1, 42.9, 52.4, 68.0],
                    [32.5, 39.3, 47.7, 60.7, 71.2], [35.2, 42.9, 54.5, 64.3, 74.0], [37.1, 46.0, 56.1, 65.1, 76.7]];
  const PT420_T2 = [[2.31, 2.28, 2.34, 2.63, 2.70], [3.02, 2.99, 2.99, 3.11, 2.89], [4.43, 4.34, 4.17, 4.24, 4.08],
                    [7.35, 6.87, 6.52, 6.81, 7.39], [9.76, 9.44, 9.51, 9.63, 9.76], [10.69, 10.63, 10.73, 10.60, 11.53]];
  function seg(xs, x) {
    let i = 0;
    while (i < xs.length - 2 && x > xs[i + 1]) i++;
    return [i, (x - xs[i]) / (xs[i + 1] - xs[i])];
  }
  function pt_map(Q1, Q2, sc) {
    const [i, u] = seg(PT420_Q1, Q1 / sc), [j, v] = seg(PT420_Q2, Q2 / sc);
    const bil = M => {
      const a = M[j][i] + u * (M[j][i + 1] - M[j][i]), b = M[j + 1][i] + u * (M[j + 1][i + 1] - M[j + 1][i]);
      return a + v * (b - a);
    };
    return [bil(PT420_T1), bil(PT420_T2)];
  }
  function pt_load2(Q1, T2, sc) {
    const [i, u] = seg(PT420_Q1, Q1 / sc), col = PT420_T2.map(r => r[i] + u * (r[i + 1] - r[i]));
    const [j, v] = seg(col, T2);
    return sc * (PT420_Q2[j] + v * (PT420_Q2[j + 1] - PT420_Q2[j]));
  }
  function pt_load1(Q2, T1, sc) {
    const [j, v] = seg(PT420_Q2, Q2 / sc), row = PT420_T1[j].map((a, k) => a + v * (PT420_T1[j + 1][k] - a));
    const [i, u] = seg(row, T1);
    return sc * (PT420_Q1[i] + u * (PT420_Q1[i + 1] - PT420_Q1[i]));
  }

  function thermal(P) {
    const p = P.p_mbar * 100, tube = Tube(P.id_mm * 1e-3, P.wall_mm * 1e-3);
    const [legs, errs] = build_legs(P);
    const mdot = P.mode === "single" ? 0 : thermosiphon(legs[0], legs[1], tube, p);
    const sol = legs.map((leg, i) => solve_leg(leg, tube, p, i === 0 ? mdot : -mdot, true));
    const heat = { cold: 0, intercept: 0, engine: 0, warm: 0 };
    sol.forEach(s => s.heats.forEach(([k, , q]) => { heat[k] += q; }));
    return [tube, legs, errs, mdot, sol, heat];
  }
  function cooler_solve(P) {
    const sc = P.cooler_pct / 100, Tset = P.Tc;
    if (!P.cooler_map) {
      const th = thermal(P), st = stage_loads(P, th[5], th[0], th[1].length);
      const Q1 = st.first.total, Q2 = st.second.total, [T1f, T2f] = pt_map(Q1, Q2, sc);
      st.first.capacity = pt_load1(Q2, P.T_shield, sc);
      st.second.capacity = pt_load2(Q1, Tset, sc);
      const cooler = { map: false, T1: P.T_shield, T2: Tset, set: Tset, trim: 0,
        held: Q1 <= st.first.capacity && Q2 <= st.second.capacity, T1_free: T1f, T2_free: T2f, Q1, Q2,
        off_map: Q1 / sc > PT420_Q1[PT420_Q1.length - 1] || Q2 / sc > PT420_Q2[PT420_Q2.length - 1] };
      return [P, th, st, cooler];
    }
    let Tc = Tset, T1 = 40, prev = null, th, st, Q1, Q2, T1f, T2f, cap2, trim, T2;
    for (let n = 0; n < 30; n++) {
      th = thermal(Object.assign({}, P, { Tc }));
      for (let m = 0; m < 100; m++) {
        st = stage_loads(Object.assign({}, P, { Tc, T_shield: T1 }), th[5], th[0], th[1].length);
        Q1 = st.first.total; Q2 = st.second.total;
        [T1f, T2f] = pt_map(Q1, Q2, sc);
        cap2 = pt_load2(Q1, Tset, sc);
        trim = T2f < Tset ? cap2 - Q2 : 0;
        let T1n;
        [T1n, T2] = trim > 0 ? pt_map(Q1, Q2 + trim, sc) : [T1f, T2f];
        if (Math.abs(T1n - T1) < 1e-9) break;
        T1 = T1n;
      }
      const r = T2 - Tc;
      if (Math.abs(r) < 1e-7) break;
      const Tn = prev && r !== prev[1] ? Tc - r * (Tc - prev[0]) / (r - prev[1]) : T2;   // secant on T2(Tc) - Tc
      prev = [Tc, r];
      Tc = Math.max(Tset, Tn);
    }
    st.first.capacity = null;
    st.second.capacity = cap2;
    const cooler = { map: true, T1, T2: Tc, set: Tset, trim, held: T2f <= Tset, T1_free: T1f, T2_free: T2f, Q1, Q2,
      off_map: Q1 / sc > PT420_Q1[PT420_Q1.length - 1] || (Q2 + trim) / sc > PT420_Q2[PT420_Q2.length - 1] };
    return [Object.assign({}, P, { Tc, T_shield: T1 }), th, st, cooler];
  }

  function compute(P0) {
    const [P, th, stages, cooler] = cooler_solve(P0);
    const [tube, legs, errs, mdot, sol, heat] = th;
    const p = P.p_mbar * 100, Tc = P.Tc, Th = P.Th;
    const G_diff = sol.reduce((s, x) => s + 1 / x.R_diff, 0), G_conv = Math.abs(mdot) / M_ATOM, G = G_diff + G_conv;
    const Np = n_num(Th, p) * P.V_pump_cc * 1e-6, Nt = n_num(Tc, p) * P.V_target_cc * 1e-6;
    const Mss = P.M_ss / 100, g_b = 1 / P.T_b_s, g_t = 1 / (P.T1t_h * 3600);
    const kp = G / Np, kt = G / Nt;
    const A_on = [[-(g_b + kp), kp], [kt, -(kt + g_t)]];
    const rates = { on: [A_on, [g_b * Mss, 0]], off: [A_on, [g_b * Mss, 0]] };
    const t_on = P.on_h * 3600;
    const [t, pp, pt, ss_on] = evolve(rates, [0, 0], t_on, t_on, 361);
    let t90 = null; const target = 0.9 * ss_on[1];
    const pt_at = tt => evolve(rates, [0, 0], 1e30, tt, 2)[2][1];
    if (ss_on[1] > 0) {
      let lo = 0, hi = 1;
      while (pt_at(hi) < target && hi < 1e9) hi *= 2;
      if (hi < 1e9) { for (let i = 0; i < 60; i++) { const mid = 0.5 * (lo + hi); if (pt_at(mid) < target) lo = mid; else hi = mid; } t90 = hi; }
    }
    const sys = two_cell(Np, Nt, G, P.tau_p_s, P.T1t_h * 3600);
    const Mp0 = pp[pp.length - 1], Mt0 = pt[pt.length - 1];
    const amp_zero = amplitudes(sys, 0, Mt0), amp_shut = amplitudes(sys, Mp0, Mt0);
    const t_win = Math.max(60, 8 * sys.tau_s), tz = [];
    for (let i = 0; i <= 360; i++) tz.push(t_win * i / 360);
    const [zp, zt] = relax_curves(sys, amp_zero, tz), [sp, st] = relax_curves(sys, amp_shut, tz);
    const rs = 1 / sys.tau_s, rl = 1 / sys.tau_l, a_s = amp_zero.a_s, a_l = amp_zero.a_l;
    const t_peak = a_s < 0 && 0 < a_l ? Math.log(-rs * a_s / (rl * a_l)) / (rs - rl) : 0;
    const { A: _A, ...paper } = sys;
    Object.assign(paper, { Mp0, Mt0, zeroed: amp_zero, shuttered: amp_shut, t_peak,
      P_peak: a_s * Math.exp(-rs * t_peak) + a_l * Math.exp(-rl * t_peak), late_ratio: a_l / amp_zero.a_lt });
    const warnings = errs.slice(), Re = 4 * Math.abs(mdot) / (Math.PI * tube.d * eta(Tc));
    if (Re > 50) warnings.push(`Reynolds number ${Re.toFixed(0)} at the cold end: gas may lag the wall temperature, so the 1-D model overestimates the flow.`);
    const g3 = v => String(+v.toPrecision(3));
    if (cooler.map) {
      if (!cooler.held) {
        const cap2 = stages.second.capacity;
        if (cap2 > 0) warnings.push(`The PT420 can't hold the target at ${g3(cooler.set)} K: it gives ${g3(cap2)} W there with ${g3(cooler.Q1)} W on the first stage, short of the ${g3(cooler.Q2)} W load, so the target runs at ${g3(Tc)} K.`);
        else warnings.push(`The PT420 can't reach ${g3(cooler.set)} K with ${g3(cooler.Q1)} W on the first stage, so the target runs at ${g3(Tc)} K.`);
      }
    } else {
      [["Second", stages.second, cooler.T2_free], ["First", stages.first, cooler.T1_free]].forEach(([name, stage, T_free]) => {
        if (stage.total > stage.capacity) {
          if (stage.capacity > 0) warnings.push(`${name}-stage load of ${g3(stage.total)} W exceeds the PT420's ${g3(stage.capacity)} W at ${g3(stage.T)} K; with these loads it settles at ${g3(T_free)} K.`);
          else warnings.push(`The PT420 can't reach ${g3(stage.T)} K on the ${name.toLowerCase()} stage with these loads; it settles at ${g3(T_free)} K.`);
        }
      });
    }
    if (cooler.off_map) {
      const sc = P.cooler_pct / 100;
      warnings.push(`The loads run past the edge of the PT420 chart (${g3(PT420_Q1[PT420_Q1.length - 1] * sc)} W on the first stage, ${g3(PT420_Q2[PT420_Q2.length - 1] * sc)} W on the second), so the stage temperatures are extrapolated.`);
    }
    if (sol.some(s => s.clamped)) warnings.push("Flow is strong enough to flatten a free section's temperature profile; treat these results as rough.");
    return {
      ok: true, engine: "javascript", mdot_mg_h: mdot * 3.6e9, Re, G, G_diff, G_conv, Np, Nt,
      N_tubes: sol.reduce((s, x) => s + x.N_tube, 0), tau_x_target_s: Nt / G, tau_x_pump_s: Np / G,
      heat_mW: Object.fromEntries(Object.entries(heat).map(([k, v]) => [k, v * 1e3])), stages, Tc, T_shield: P.T_shield, cooler,
      legs: sol.map((s, i) => ({ profile: s.profile, length_cm: s.length * 100, nodes: legs[i].nodes.map(n => [n[0], n[1] * 100, n[2]]), free_cm: legs[i].free.map(f => f * 100), G_diff: 1 / s.R_diff })),
      ss_pump: ss_on[0], ss_target: ss_on[1], times_s: eig_times(A_on), t90_target_s: t90,
      production_atoms_s: g_b * Mss * Np,
      loss_exchange: (Nt / G) * g_t, loss_capacity: (Nt / Np) * P.T_b_s * g_t,
      t_h: t.map(x => x / 3600), P_pump: pp, P_target: pt, paper,
      zero_t_s: tz, zero_P_pump: zp, zero_P_target: zt, shut_P_pump: sp, shut_P_target: st, warnings,
    };
  }
  function run(json) { try { return JSON.stringify(compute(JSON.parse(json))); } catch (e) { return JSON.stringify({ ok: false, error: String(e) }); } }
  const PT420 = { Q1: PT420_Q1, Q2: PT420_Q2, T1: PT420_T1, T2: PT420_T2, map: pt_map, load1: pt_load1, load2: pt_load2 };
  root.He3ModelJS = { run, compute, PT420 };
  if (typeof module !== "undefined") module.exports = { run, compute };
})(typeof window !== "undefined" ? window : globalThis);
