"""
he3_model.py -- two-cell polarized 3He exchange model (pure Python, Pyodide-ready)

Warm MEOP pumping cell (Th) <-> glass leg(s) <-> cold target cell (Tc), with the
cells filled to the operating pressure at temperature, then valved off from the supply.
Two-cell target concept: Maxwell & Milner, arXiv:2106.06606.

Transport (per leg, 1-D along the tube):
  * glass wall + gas conduction, plus the enthalpy carried by any net flow.
    In each free (unclamped) section the energy flux F = -K dT/ds + mdot*cp*T is
    constant, so ds/dT = K(T) / (mdot*cp*T - F); F is found by Newton iteration.
  * clamped sections (cell ends, cold-end "engine", intercepts) are isothermal;
    the heat each clamp removes is the jump in F across it.
  * diffusion of polarization at uniform pressure: G = A / int ds/(n D);
    n*D = 1.32 eta(T)/m is independent of pressure.
  * thermosiphon (loop mode): buoyancy g*(int rho ds)_down - (int rho ds)_up
    balanced by Poiseuille friction 128 nu mdot / (pi d^4), with the tube
    temperature profiles re-solved for each trial flow.

Polarization (exact solution of the linear 2x2 system, piecewise in time):
  dPp/dt = (Mss - Pp)/Tb - (G/Np)(Pp - Pt)          laser and discharge on
  dPp/dt = -Pp/tau_p     - (G/Np)(Pp - Pt)          laser shuttered, discharge on
  dPt/dt = (G/Nt)(Pp - Pt) - Pt/tau_t
Mss and Tb are the steady polarization and build-up time measured in a sealed
cell at the operating pressure and field (Nikiel et al.: dM/dt = Gop (Mss - M)).
In steady state the target reaches Pt = Mss / (1 + a + b), with
  a = tau_x / T1t            exchange too slow        (tau_x = Nt/G)
  b = (Nt/Np) (Tb / T1t)     MEOP can't keep up with target relaxation.

Pump-cell probe (Maxwell, Epstein & Milner, arXiv:1412.6167): after the build-up
the laser is shuttered and the pumping cell is zeroed; the discharge stays on so
the optical polarimeter can watch polarized gas return from the target. With the
laser off the system is eq. 2 of that paper,
  dMp/dt = -(1/tau_p + (Nt/N)/t_ex) Mp + (Nt/N)/t_ex Mt
  dMt/dt = (Np/N)/t_ex Mp - (1/tau_t + (Np/N)/t_ex) Mt
with N = Np + Nt and t_ex = Np Nt / (N G), and its solution (eq. 3)
  Mp(t) = a_s exp(-t/tau_s) + a_l exp(-t/tau_l)
  Mt(t) = a_lt exp(-t/tau_l) - a_st exp(-t/tau_s)
with tau_s, tau_l from eq. 4. Once the short term has died, Mp/Mt = a_l/a_lt.

Entry point: run(params_json) -> result_json. Edit the property functions freely.
"""
import math
import json

# ------------------------------------------------------------------ constants
KB = 1.380649e-23
NA = 6.02214076e23
M3 = 3.01603e-3                 # kg/mol
M_ATOM = M3 / NA                # kg/atom
RG = 8.314462
CP = 2.5 * RG / M3              # J/kg/K, monatomic ideal gas
GRAV = 9.81

# ------------------------------------------------------------------ 3He gas
ETA_300 = 17.2e-6               # Pa s at 300 K (4He value * sqrt(3/4))
ETA_EXP = 0.68                  # eta ~ T^0.68; likely low by tens of % at 5 K
RHO_D_OVER_ETA = 1.32           # Chapman-Enskog self-diffusion ratio for He


def eta(T):
    return ETA_300 * (T / 300.0) ** ETA_EXP


def k_gas(T):                   # Eucken relation, monatomic
    return 3.75 * RG / M3 * eta(T)


def n_num(T, p):
    return p / (KB * T)


def rho(T, p):
    return p * M_ATOM / (KB * T)


def nD(T):                      # n*D, pressure independent [1/(m s)]
    return RHO_D_OVER_ETA * eta(T) / M_ATOM


# ------------------------------------------------------------------ glass
# Approximate borosilicate thermal conductivity [W/m/K]; swap in NIST data.
GLASS_T = [4.0, 6.0, 10.0, 20.0, 30.0, 50.0, 77.0, 100.0, 150.0, 200.0, 250.0, 300.0]
GLASS_K = [0.10, 0.11, 0.12, 0.15, 0.20, 0.32, 0.45, 0.55, 0.75, 0.92, 1.03, 1.12]


def k_glass(T):
    if T <= GLASS_T[0]:
        return GLASS_K[0]
    if T >= GLASS_T[-1]:
        return GLASS_K[-1]
    x = math.log(T)
    for i in range(1, len(GLASS_T)):
        x1 = math.log(GLASS_T[i])
        if x <= x1:
            x0 = math.log(GLASS_T[i - 1])
            f = (x - x0) / (x1 - x0)
            return math.exp(math.log(GLASS_K[i - 1]) + f * (math.log(GLASS_K[i]) - math.log(GLASS_K[i - 1])))
    return GLASS_K[-1]


# ------------------------------------------------------------------ numerics
def trap(y, x):
    s = 0.0
    for i in range(1, len(x)):
        s += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1])
    return s


def cumtrap(y, x):
    out = [0.0]
    s = 0.0
    for i in range(1, len(x)):
        s += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1])
        out.append(s)
    return out


GRID_N = 120


def t_grid(lo, hi):
    """Temperature grid dense near both ends (where profiles can be singular)."""
    span = hi - lo
    pts = set([lo, hi])
    for i in range(GRID_N):
        f = 10.0 ** (-9.0 + 9.0 * i / (GRID_N - 1))
        pts.add(lo + f * span)
        pts.add(hi - f * span)
        pts.add(lo * (hi / lo) ** (i / (GRID_N - 1)))
    return sorted(t for t in pts if lo <= t <= hi)


_cache = {}


def grid_K(lo, hi, A_wall, A_gas):
    key = (round(lo, 9), round(hi, 9), A_wall, A_gas)
    g = _cache.get(key)
    if g is None:
        T = t_grid(lo, hi)
        K = [A_wall * k_glass(t) + A_gas * k_gas(t) for t in T]
        g = (T, K, trap(K, T))
        if len(_cache) > 200:
            _cache.clear()
        _cache[key] = g
    return g


class Tube:
    def __init__(self, d_in, wall):
        self.d = d_in
        r = d_in / 2.0
        self.A_gas = math.pi * r * r
        self.A_wall = math.pi * ((r + wall) ** 2 - r * r)


def solve_free(L, lo, hi, mdot, tube):
    """Free section from lo (x=0) to hi (x=L), flow mdot along +x.
    Returns T grid, ds/dT weights, energy flux F along +x, clamped flag."""
    T, K, IK = grid_K(lo, hi, tube.A_wall, tube.A_gas)
    a = mdot * CP
    n = len(T)
    if abs(a) < 1e-30:
        F = -IK / L
        return T, [k / (-F) for k in K], F, False

    def LdL(F):
        s = s2 = 0.0
        g0 = K[0] / (a * T[0] - F)
        for i in range(1, n):
            g1 = K[i] / (a * T[i] - F)
            dx = 0.5 * (T[i] - T[i - 1])
            s += (g0 + g1) * dx
            s2 += (g0 * g0 / K[i - 1] + g1 * g1 / K[i]) * dx    # dL/dF
            g0 = g1
        return s, s2

    Fmax = min(a * lo, a * hi)
    scale = IK / L
    Fhi = Fmax - 1e-12 * max(scale, abs(Fmax))
    clamped = False
    if LdL(Fhi)[0] < L:
        F = Fhi
        clamped = True
    else:
        Flo = Fmax - scale
        while LdL(Flo)[0] > L:
            Flo = Fmax - 4.0 * (Fmax - Flo)
        lo_b, hi_b = Flo, Fhi
        F = Flo
        for _ in range(100):
            v, dv = LdL(F)
            r = v - L
            if abs(r) < 1e-10 * L:
                break
            if r < 0:
                lo_b = F
            else:
                hi_b = F
            Fn = F - r / dv if dv > 0 else 0.5 * (lo_b + hi_b)
            if not (lo_b < Fn < hi_b):
                Fn = 0.5 * (lo_b + hi_b)
            F = Fn
    w = [K[i] / (a * T[i] - F) for i in range(n)]
    norm = L / trap(w, T)
    return T, [x * norm for x in w], F, clamped


def free_segment(L, T0, T1, mdot_s, tube):
    """Section running from T0 at its start (s=0) to T1 at its end (s=L).
    mdot_s > 0 flows toward the end. Returns T, weights w=ds/dT, s(T), F along +s."""
    if abs(T1 - T0) < 1e-6:
        return dict(iso=True, T=T0, F=mdot_s * CP * T0, clamped=False)
    if T1 > T0:
        T, w, F, cl = solve_free(L, T0, T1, mdot_s, tube)
        return dict(iso=False, T=T, w=w, s=cumtrap(w, T), F=F, clamped=cl)
    T, w, F, cl = solve_free(L, T1, T0, -mdot_s, tube)
    x = cumtrap(w, T)
    return dict(iso=False, T=T, w=w, s=[L - xi for xi in x], F=-F, clamped=cl)


def solve_leg(leg, tube, p, mdot, want_profile=False):
    """leg: dict(nodes=[(T, clamp_length, kind), ...], free=[lengths]).
    Node 0 sits at the target cell, the last node at the pumping cell.
    mdot > 0 = upward flow (target -> pumping cell)."""
    nodes, frees = leg["nodes"], leg["free"]
    a = mdot * CP
    I_rho = I_nu = I_rd = I_n = 0.0
    prof = []
    s0 = 0.0
    Fs = []
    clamped = False
    for i, (Tn, cl, kind) in enumerate(nodes):
        if cl > 0:
            I_rho += rho(Tn, p) * cl
            I_nu += eta(Tn) / rho(Tn, p) * cl
            I_rd += cl / nD(Tn)
            I_n += n_num(Tn, p) * cl
        if want_profile:
            prof.append((s0, Tn))
            prof.append((s0 + cl, Tn))
        s0 += cl
        if i < len(frees):
            L = frees[i]
            seg = free_segment(L, Tn, nodes[i + 1][0], mdot, tube)
            if seg["iso"]:
                t = seg["T"]
                I_rho += rho(t, p) * L
                I_nu += eta(t) / rho(t, p) * L
                I_rd += L / nD(t)
                I_n += n_num(t, p) * L
            else:
                T, w = seg["T"], seg["w"]
                I_rho += trap([rho(t, p) * wi for t, wi in zip(T, w)], T)
                I_nu += trap([eta(t) / rho(t, p) * wi for t, wi in zip(T, w)], T)
                I_rd += trap([wi / nD(t) for t, wi in zip(T, w)], T)
                I_n += trap([n_num(t, p) * wi for t, wi in zip(T, w)], T)
                clamped = clamped or seg["clamped"]
                if want_profile:
                    pts = sorted(zip(seg["s"], T))
                    m = len(pts)
                    step = max(1, m // 60)
                    for j in range(0, m, step):
                        prof.append((s0 + pts[j][0], pts[j][1]))
                    prof.append((s0 + pts[-1][0], pts[-1][1]))
            Fs.append(seg["F"])
            s0 += L
    heats = []
    for i, (Tn, cl, kind) in enumerate(nodes):
        Fin = a * Tn if i == 0 else Fs[i - 1]
        Fout = a * Tn if i == len(nodes) - 1 else Fs[i]
        heats.append((kind, Tn, Fin - Fout))
    A = tube.A_gas
    out = dict(
        I_rho=I_rho,
        R_fric=128.0 / (math.pi * tube.d ** 4) * I_nu,
        R_diff=I_rd / A,
        N_tube=A * I_n,
        heats=heats,
        length=s0,
        clamped=clamped,
    )
    if want_profile:
        prof.sort()
        out["profile"] = [[round(s * 100, 4), round(t, 3)] for s, t in prof]
    return out


def thermosiphon(legA, legB, tube, p):
    """legA = riser, legB = downcomer, same total height, pumping cell on top.
    Returns mdot (> 0: up A, down B)."""

    def resid(md):
        A = solve_leg(legA, tube, p, md)
        B = solve_leg(legB, tube, p, -md)
        drive = GRAV * (B["I_rho"] - A["I_rho"])
        return md * (A["R_fric"] + B["R_fric"]) - drive, drive, A["R_fric"] + B["R_fric"]

    r0, d0, R0 = resid(0.0)
    if abs(d0) < 1e-15:
        return 0.0
    sgn = 1.0 if d0 > 0 else -1.0
    hi = sgn * abs(d0) / R0
    rhi = resid(hi)[0]
    while rhi * sgn < 0:
        hi *= 2.0
        rhi = resid(hi)[0]
    lo, rlo = 0.0, r0
    md = hi
    for _ in range(80):                      # Illinois false position
        md = (lo * rhi - hi * rlo) / (rhi - rlo)
        r = resid(md)[0]
        if abs(r) < 1e-9 * abs(d0) or abs(hi - lo) < 1e-9 * abs(md):
            break
        if (r > 0) == (rhi > 0):
            hi, rhi = md, r
            rlo *= 0.5
        else:
            lo, rlo = md, r
            rhi *= 0.5
    return md


# ------------------------------------------------------------------ geometry
def build_legs(P):
    cm = 0.01
    Tc, Th = P["Tc"], P["Th"]
    H = P["length_cm"] * cm
    use_int = P["intercept"] and P["mode"] in ("single", "loop")
    Ti, si = P["T_int"], P["s_int_cm"] * cm
    errs = []

    def cold_leg(clamp):
        if use_int:
            if not (clamp + 0.2 * cm < si < H - 0.2 * cm):
                errs.append("Intercept height must sit between the cold clamp and the pumping cell.")
                return dict(nodes=[(Tc, clamp, "cold"), (Th, 0.0, "warm")], free=[H - clamp])
            return dict(nodes=[(Tc, clamp, "cold"), (Ti, 0.0, "intercept"), (Th, 0.0, "warm")],
                        free=[si - clamp, H - si])
        return dict(nodes=[(Tc, clamp, "cold"), (Th, 0.0, "warm")], free=[H - clamp])

    if P["mode"] == "single":
        return [cold_leg(0.0)], errs
    he, gap, Te = P["engine_cm"] * cm, P["gap_cm"] * cm, P["T_engine"]
    if H - gap - he < 0.5 * cm:
        errs.append("Legs are too short for the engine section; lengthen them or shorten the engine.")
        he = max(0.0, H - gap - 0.5 * cm)
    riser = dict(nodes=[(Tc, 0.0, "cold"), (Te, he, "engine"), (Th, 0.0, "warm")],
                 free=[gap, H - gap - he])
    return [riser, cold_leg(gap + he)], errs


# ------------------------------------------------------------------ evolution
def expm2(A, t):
    a, b, c, d = A[0][0], A[0][1], A[1][0], A[1][1]
    m = 0.5 * (a + d)
    det = a * d - b * c
    q2 = m * m - det
    q = math.sqrt(q2) if q2 > 0 else 0.0
    if q * t < 1e-8:
        e = math.exp(m * t)
        return [[e * (1 + (a - m) * t), e * b * t], [e * c * t, e * (1 + (d - m) * t)]]
    l1, l2 = m + q, m - q
    e1, e2 = math.exp(l1 * t), math.exp(l2 * t)
    k = 1.0 / (l1 - l2)
    return [[k * (e1 * (a - l2) - e2 * (a - l1)), k * (e1 - e2) * b],
            [k * (e1 - e2) * c, k * (e1 * (d - l2) - e2 * (d - l1))]]


def evolve(rates, P0, t_on, t_end, npts):
    """Piecewise-exact evolution. rates[k] = (A, b) for discharge on / off."""
    def ss(A, b):
        det = A[0][0] * A[1][1] - A[0][1] * A[1][0]
        return [-(A[1][1] * b[0] - A[0][1] * b[1]) / det, -(-A[1][0] * b[0] + A[0][0] * b[1]) / det]

    t = [t_end * i / (npts - 1) for i in range(npts)]
    if 0 < t_on < t_end and t_on not in t:
        t.append(t_on)
        t.sort()
    out_p, out_t = [], []
    A_on, b_on = rates["on"]
    A_off, b_off = rates["off"]
    ss_on, ss_off = ss(A_on, b_on), ss(A_off, b_off)
    x_switch = None
    for ti in t:
        if ti <= t_on:
            A, xs, x0, dt = A_on, ss_on, P0, ti
        else:
            if x_switch is None:
                E = expm2(A_on, t_on)
                d0 = [P0[0] - ss_on[0], P0[1] - ss_on[1]]
                x_switch = [ss_on[0] + E[0][0] * d0[0] + E[0][1] * d0[1],
                            ss_on[1] + E[1][0] * d0[0] + E[1][1] * d0[1]]
            A, xs, x0, dt = A_off, ss_off, x_switch, ti - t_on
        E = expm2(A, dt)
        d0 = [x0[0] - xs[0], x0[1] - xs[1]]
        out_p.append(xs[0] + E[0][0] * d0[0] + E[0][1] * d0[1])
        out_t.append(xs[1] + E[1][0] * d0[0] + E[1][1] * d0[1])
    return t, out_p, out_t, ss_on


def eig_times(A):
    a, b, c, d = A[0][0], A[0][1], A[1][0], A[1][1]
    m = 0.5 * (a + d)
    q = math.sqrt(max(m * m - (a * d - b * c), 0.0))
    return [-1.0 / (m + q), -1.0 / (m - q)]      # slow, fast


def two_cell(Np, Nt, G, tau_p, tau_t):
    """Laser-off rate equations in the terms of arXiv:1412.6167, eqs. 2 and 4."""
    N = Np + Nt
    t_ex = Np * Nt / (N * G)
    A = [[-(1.0 / tau_p + Nt / N / t_ex), Nt / N / t_ex],
         [Np / N / t_ex, -(1.0 / tau_t + Np / N / t_ex)]]
    dg = 1.0 / tau_p - 1.0 / tau_t
    mean = 0.5 * (1.0 / tau_p + 1.0 / tau_t + 1.0 / t_ex)
    root = 0.5 * math.sqrt(1.0 / t_ex ** 2 + dg * dg + 2.0 / t_ex * (Nt - Np) / N * dg)
    return dict(N=N, Np_N=Np / N, Nt_N=Nt / N, t_ex=t_ex, tau_p=tau_p, tau_t=tau_t,
                tau_s=1.0 / (mean + root), tau_l=1.0 / (mean - root), A=A)


def amplitudes(sys2, Mp0, Mt0):
    """a_s, a_l, a_st, a_lt of eq. 3 for starting polarizations Mp0, Mt0,
    from matching the values and slopes at t = 0."""
    A = sys2["A"]
    rs, rl = 1.0 / sys2["tau_s"], 1.0 / sys2["tau_l"]
    dMp = A[0][0] * Mp0 + A[0][1] * Mt0
    dMt = A[1][0] * Mp0 + A[1][1] * Mt0
    a_l = (dMp + rs * Mp0) / (rs - rl)
    a_lt = (dMt + rs * Mt0) / (rs - rl)
    return dict(a_s=Mp0 - a_l, a_l=a_l, a_st=a_lt - Mt0, a_lt=a_lt)


def relax_curves(sys2, amp, t):
    rs, rl = 1.0 / sys2["tau_s"], 1.0 / sys2["tau_l"]
    Mp = [amp["a_s"] * math.exp(-rs * x) + amp["a_l"] * math.exp(-rl * x) for x in t]
    Mt = [amp["a_lt"] * math.exp(-rl * x) - amp["a_st"] * math.exp(-rs * x) for x in t]
    return Mp, Mt


# ------------------------------------------------------------------ entry
def compute(P):
    p = P["p_mbar"] * 100.0
    Tc, Th = P["Tc"], P["Th"]
    tube = Tube(P["id_mm"] * 1e-3, P["wall_mm"] * 1e-3)
    legs, errs = build_legs(P)

    if P["mode"] == "single":
        mdot = 0.0
    else:
        mdot = thermosiphon(legs[0], legs[1], tube, p)
    sol = []
    for i, leg in enumerate(legs):
        md = mdot if i == 0 else -mdot
        sol.append(solve_leg(leg, tube, p, md, want_profile=True))

    G_diff = sum(1.0 / s["R_diff"] for s in sol)
    G_conv = abs(mdot) / M_ATOM
    G = G_diff + G_conv

    Np = n_num(Th, p) * P["V_pump_cc"] * 1e-6
    Nt = n_num(Tc, p) * P["V_target_cc"] * 1e-6

    heat = dict(cold=0.0, intercept=0.0, engine=0.0, warm=0.0)
    for s in sol:
        for kind, Tn, q in s["heats"]:
            heat[kind] += q

    Mss = P["M_ss"] / 100.0
    g_b = 1.0 / P["T_b_s"]
    g_t = 1.0 / (P["T1t_h"] * 3600.0)
    kp, kt = G / Np, G / Nt

    # laser and discharge on: build-up from zero
    A_on = [[-(g_b + kp), kp], [kt, -(kt + g_t)]]
    rates = dict(on=(A_on, [g_b * Mss, 0.0]), off=(A_on, [g_b * Mss, 0.0]))
    t_on = P["on_h"] * 3600.0
    t, pp, pt, ss_on = evolve(rates, [0.0, 0.0], t_on, t_on, 361)

    # 90% rise time of the target cell under pumping
    t90 = None
    target = 0.9 * ss_on[1]
    lo_t, hi_t = 0.0, 1.0

    def pt_at(tt):
        return evolve(rates, [0.0, 0.0], 1e30, tt, 2)[2][-1]
    if ss_on[1] > 0:
        while pt_at(hi_t) < target and hi_t < 1e9:
            hi_t *= 2.0
        if hi_t < 1e9:
            for _ in range(60):
                mid = 0.5 * (lo_t + hi_t)
                if pt_at(mid) < target:
                    lo_t = mid
                else:
                    hi_t = mid
            t90 = hi_t

    # laser shuttered at the end of the build-up, with and without zeroing the pumping cell
    sys2 = two_cell(Np, Nt, G, P["tau_p_s"], P["T1t_h"] * 3600.0)
    Mp0, Mt0 = pp[-1], pt[-1]
    amp_zero = amplitudes(sys2, 0.0, Mt0)
    amp_shut = amplitudes(sys2, Mp0, Mt0)
    t_win = max(60.0, 8.0 * sys2["tau_s"])
    tz = [t_win * i / 360 for i in range(361)]
    zp, zt = relax_curves(sys2, amp_zero, tz)
    sp, st = relax_curves(sys2, amp_shut, tz)
    rs, rl = 1.0 / sys2["tau_s"], 1.0 / sys2["tau_l"]
    a_s, a_l = amp_zero["a_s"], amp_zero["a_l"]
    t_peak = math.log(-rs * a_s / (rl * a_l)) / (rs - rl) if a_s < 0 < a_l else 0.0
    paper = {k: v for k, v in sys2.items() if k != "A"}
    paper.update(Mp0=Mp0, Mt0=Mt0, zeroed=amp_zero, shuttered=amp_shut,
                 t_peak=t_peak, P_peak=a_s * math.exp(-rs * t_peak) + a_l * math.exp(-rl * t_peak),
                 late_ratio=a_l / amp_zero["a_lt"])

    warnings = list(errs)
    Re = 4.0 * abs(mdot) / (math.pi * tube.d * eta(Tc))
    if Re > 50:
        warnings.append("Reynolds number %.0f at the cold end: gas may lag the wall temperature, so the 1-D model overestimates the flow." % Re)
    if heat["cold"] > 2.0:
        warnings.append("Cold-stage load exceeds the PT420 rating of about 2 W at 4.2 K.")
    if any(s["clamped"] for s in sol):
        warnings.append("Flow is strong enough to flatten a free section's temperature profile; treat these results as rough.")

    return dict(
        ok=True,
        engine="python",
        mdot_mg_h=mdot * 3.6e9,
        Re=Re,
        G=G, G_diff=G_diff, G_conv=G_conv,
        Np=Np, Nt=Nt, N_tubes=sum(s["N_tube"] for s in sol),
        tau_x_target_s=Nt / G, tau_x_pump_s=Np / G,
        heat_mW={k: v * 1e3 for k, v in heat.items()},
        legs=[dict(profile=s["profile"], length_cm=s["length"] * 100,
                   nodes=[[n[0], n[1] * 100, n[2]] for n in leg["nodes"]],
                   free_cm=[f * 100 for f in leg["free"]],
                   G_diff=1.0 / s["R_diff"])
              for s, leg in zip(sol, legs)],
        ss_pump=ss_on[0], ss_target=ss_on[1],
        times_s=eig_times(A_on),
        t90_target_s=t90,
        production_atoms_s=g_b * Mss * Np,
        loss_exchange=(Nt / G) * g_t, loss_capacity=(Nt / Np) * P["T_b_s"] * g_t,
        t_h=[x / 3600.0 for x in t], P_pump=pp, P_target=pt,
        paper=paper,
        zero_t_s=tz, zero_P_pump=zp, zero_P_target=zt, shut_P_pump=sp, shut_P_target=st,
        warnings=warnings,
    )


def run(params_json):
    try:
        return json.dumps(compute(json.loads(params_json)))
    except Exception as e:  # report, don't crash the page
        return json.dumps(dict(ok=False, error="%s: %s" % (type(e).__name__, e)))


DEFAULTS = dict(
    mode="loop", p_mbar=100.0, Tc=5.0, Th=300.0,
    id_mm=8.0, wall_mm=1.0, length_cm=10.0,
    engine_cm=2.0, gap_cm=0.5, T_engine=25.0,
    intercept=True, T_int=45.0, s_int_cm=4.0,
    V_pump_cc=49.0, V_target_cc=100.0,
    M_ss=50.0, T_b_s=150.0, tau_p_s=600.0, T1t_h=10.0,
    on_h=24.0,
)

if __name__ == "__main__" and __import__("sys").platform != "emscripten":
    import time
    for mode in ("single", "loop"):
        P = dict(DEFAULTS, mode=mode)
        t0 = time.time()
        r = compute(P)
        print(mode, "%.2fs" % (time.time() - t0),
              {k: r[k] for k in ("mdot_mg_h", "Re", "G", "G_diff", "tau_x_target_s", "ss_pump", "ss_target",
                                  "times_s", "t90_target_s", "heat_mW", "warnings")})
        print("  paper", {k: v for k, v in r["paper"].items()})
