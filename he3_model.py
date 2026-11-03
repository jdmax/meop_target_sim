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
    the heat each clamp removes is the jump in F across it. The engine can be off
    (the downcomer intercept alone then drives the loop) or held by the first stage
    instead of a heater, in which case it sits at the first-stage temperature.
  * diffusion of polarization at uniform pressure: G = A / int ds/(n D);
    n*D = 1.32 eta(T)/m is independent of pressure.
  * thermosiphon (loop mode): buoyancy g*(int rho ds)_down - (int rho ds)_up
    balanced by Poiseuille friction 128 nu mdot / (pi d^4), with the tube
    temperature profiles re-solved for each trial flow.

Cooler stages (stage_loads): the intercept, 300 K radiation onto the shield and other
loads go to the first stage; the cold clamps and target cell, an engine clamp that
removes heat, shield radiation onto the cold surfaces, direct radiation down the
bores and other loads go to the second stage. Radiation is sigma eps A (Th^4 - Tc^4).
With the PT420 map on (cooler_solve), the first stage floats at the temperature the
Cryomech capacity chart gives for the two loads, and the second stage is held at the
set point Tc by a trim heater, or floats above it if the loads alone keep it warmer.

Polarization (exact solution of the linear 2x2 system, piecewise in time):
  dPp/dt = (Mss - Pp)/Tb - (G/Np)(Pp - Pt)          laser and discharge on
  dPp/dt = -Pp/tau_p     - (G/Np)(Pp - Pt)          laser shuttered, discharge on
  dPt/dt = (G/Nt)(Pp - Pt) - Pt/tau_t
Mss, Tb and tau_p here are the pumping cell's effective values, from a sealed cell's fit
M(t) = M_ss (1 - exp(-t/T_b)) at the operating pressure and field (Nikiel et al.), with
T1p its relaxation with the discharge on and the laser off. The fit is read as fluxes:
  Q0 = N_s M_ss/T_b        optical production at P = 0 [atoms/s],  N_s = n(Th, p) V_seal
  R_dis = N_s/T1p          discharge relaxation [atoms/s per unit P]
  P_lim = M_ss/(Gop T_b),  Gop = 1/T_b - 1/T1p.
At high field the discharge, and the pump light matched to it, fill a thin layer at the
walls, and diffusion spreads polarization through the rest of the pumping cell in well
under a second. So both fluxes are set by the discharge, not by the pumping cell's volume. The
light is taken to keep up with the discharge, so both scale with its size relative to the sealed cell, s:
  Np dPp/dt = Q0' (1 - Pp/P_lim) - R' Pp - G (Pp - Pt),   Q0' = s Q0,  R' = s R_dis + Np/T1_dark
which is the form above with 1/Tb = (Q0'/P_lim + R')/Np = (s N_s/T_b + Np/T1_dark)/Np,
Mss = Q0'/(Q0'/P_lim + R'), tau_p = Np/R'. With no dark relaxation, Mss = M_ss whatever Np is. pump_scaling = "volume"
keeps the older reading, the fit as rates (Mss = M_ss, Tb = T_b), where production grows with Np.
In steady state the target reaches Pt = Mss / (1 + a + b), with
  a = tau_x / T1t            exchange too slow        (tau_x = Nt/G)
  b = (Nt/Np) (Tb / T1t)     MEOP can't keep up with target relaxation
    = Nt Mss / (Q0' T1t),    independent of Np when Mss is.

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
SIGMA = 5.670374419e-8          # W/m^2/K^4
HNU = 6.62607015e-34 * 2.99792458e8 / 1083.0e-9     # J, a 1083 nm pump photon

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
    """Riser first, then the downcomer, each from the target cell (s = 0) to the pumping cell (s = H).
    The riser carries the engine clamp unless it is off; the intercept sits on the downcomer."""
    cm = 0.01
    Tc, Th = P["Tc"], P["Th"]
    H = P["length_cm"] * cm
    Ti, si = P["T_int"], P["s_int_cm"] * cm
    loop = P["mode"] == "loop"
    eng = loop and P["engine"]
    he, gap, Te = P["engine_cm"] * cm, P["gap_cm"] * cm, P["T_engine"]
    errs = []
    if eng and H - gap - he < 0.5 * cm:
        errs.append("Legs are too short for the engine section; lengthen them or shorten the engine.")
        he = max(0.0, H - gap - 0.5 * cm)

    def leg(clamp, engine, intercept):
        nodes, free, s = [(Tc, clamp, "cold")], [], clamp
        if engine:
            nodes.append((Te, he, "engine"))
            free.append(gap)
            s += gap + he
        if intercept:
            if s + 0.2 * cm < si < H - 0.2 * cm:
                nodes.append((Ti, 0.0, "intercept"))
                free.append(si - s)
                s = si
            else:
                errs.append("Intercept height must sit between the %s and the pumping cell." % ("engine" if engine else "cold clamp"))
        nodes.append((Th, 0.0, "warm"))
        free.append(H - s)
        return dict(nodes=nodes, free=free)

    if not loop:
        return [leg(0.0, False, P["intercept"])], errs
    return [leg(0.0, eng, False), leg(gap + he if eng else 0.0, False, P["intercept"])], errs


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


def ss(A, b):
    """Steady state of dx/dt = A x + b."""
    det = A[0][0] * A[1][1] - A[0][1] * A[1][0]
    return [-(A[1][1] * b[0] - A[0][1] * b[1]) / det, -(-A[1][0] * b[0] + A[0][0] * b[1]) / det]


def evolve(rates, P0, t_on, t_end, npts):
    """Piecewise-exact evolution. rates[k] = (A, b) for discharge on / off."""
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


# ------------------------------------------------------------------ gas temperature lag (diagnostic)
NU_LAMINAR = 3.66               # fully developed laminar flow, wall at fixed temperature
LAG_STEPS = 2000


def gas_lag(sol, legs, mdot, Tc):
    """How far the flowing gas lags the tube wall. The model takes the gas at the wall temperature; here the gas
    instead exchanges heat with the wall profile at q' = Nu pi k_gas (Tg - Tw), so along the flow
    dTg/dx = -(Tg - Tw) / lam, lam = |mdot| cp / (Nu pi k_gas(film)), starting from the cell it leaves.
    The wall profile is held fixed, so this is an estimate, not part of the solution."""
    if mdot == 0.0:
        return None
    mcp = abs(mdot) * CP
    lam_of = lambda T: mcp / (NU_LAMINAR * math.pi * k_gas(T))
    out = []
    for i, (s, leg) in enumerate(zip(sol, legs)):
        prof, L = s["profile"], s["length"] * 100.0          # profile: [s cm, T K] from the target cell
        n, m = LAG_STEPS, len(s["profile"])
        xs = [L * k / n for k in range(n + 1)]
        tw, j = [], 1
        for x in xs:                                         # wall temperature on the grid, as the page samples it
            while j < m - 1 and x > prof[j][0]:
                j += 1
            (s0, t0), (s1, t1) = prof[j - 1], prof[j]
            tw.append(t1 if s1 <= s0 else t0 + (t1 - t0) * (x - s0) / (s1 - s0))
        up = (i == 0) == (mdot > 0)                          # gas moves toward the pumping cell
        order = range(n + 1) if up else range(n, -1, -1)
        tg = [0.0] * (n + 1)
        g, dx = None, L / n * 0.01
        for k in order:
            g = tw[k] if g is None else tw[k] + (g - tw[k]) * math.exp(-dx / lam_of(0.5 * (g + tw[k])))
            tg[k] = g
        kmax = max(range(n + 1), key=lambda k: abs(tg[k] - tw[k]))
        marks, x = [], 0.0
        for ni, (Tn, cl, kind) in enumerate(leg["nodes"]):  # the gas where it meets each clamp
            for xm in ((x, x + cl * 100) if cl > 0 else (x,)):
                k = min(n, max(0, int(xm / L * n + 0.5)))
                marks.append([kind, xs[k], tw[k], tg[k]])
            x += cl * 100 + (leg["free"][ni] * 100 if ni < len(leg["free"]) else 0.0)
        end = n if up else 0
        out.append(dict(up=up, max_K=tg[kmax] - tw[kmax], s_max_cm=xs[kmax], T_wall_max=tw[kmax],
                        arrive_K=tg[end], arrive_wall_K=tw[end], marks=marks,
                        profile=[[xs[k], tw[k], tg[k]] for k in range(0, n + 1, n // 100)]))
    down = [d for d in out if not d["up"]]
    past = sum(mcp * (mk[3] - mk[2]) for d in down for mk in d["marks"] if mk[0] == "intercept")
    return dict(mcp_W_K=mcp, lam_cold_mm=lam_of(Tc) * 1e3, legs=out, past_intercept_W=past)


# ------------------------------------------------------------------ cooler stages
def disk_view(r, L):
    """View factor between two equal coaxial disks of radius r a distance L apart:
    the direct line of sight from a tube's warm end to its cold end."""
    R = r / L
    X = 2.0 + 1.0 / (R * R)
    return 2.0 / (X + math.sqrt(X * X - 4.0))     # = (X - sqrt(X^2 - 4))/2, without the cancellation


def stage_loads(P, heat, tube, n_legs):
    """Split the heat into the two cooler stages [W].
    First stage (T_shield): the intercept clamp, 300 K radiation onto the shield, other loads.
    Second stage (Tc): the cold clamps and target cell, an engine clamp that removes heat,
    shield radiation onto the cold surfaces, direct radiation down the tube bores, other loads.
    A clamp that has to be heated (negative heat) is reported as a heater, not a cooler load.
    An engine held by the first stage puts its heat on that stage, with either sign: when the riser
    section needs heat, the stage supplies it and its load drops by that much."""
    Th, Tc, T1 = P["Th"], P["Tc"], P["T_shield"]
    q_int, q_eng = heat["intercept"], heat["engine"]
    on1 = engine_on_stage(P)
    rad1 = SIGMA * P["eps_shield"] * P["A_shield_m2"] * (Th ** 4 - T1 ** 4)
    rad2 = SIGMA * P["eps_cold"] * P["A_cold_cm2"] * 1e-4 * (T1 ** 4 - Tc ** 4)
    bore = n_legs * SIGMA * tube.A_gas * disk_view(tube.d / 2.0, P["length_cm"] * 0.01) * (Th ** 4 - Tc ** 4)
    first = dict(T=T1, intercept=max(0.0, q_int), engine=q_eng if on1 else 0.0, radiation=rad1, other=P["other1_W"])
    first["total"] = first["intercept"] + first["engine"] + rad1 + first["other"]
    second = dict(T=Tc, conduction=heat["cold"], engine=0.0 if on1 else max(0.0, q_eng), radiation=rad2, bore=bore,
                  other=P["other2_mW"] * 1e-3)
    second["total"] = second["conduction"] + second["engine"] + rad2 + bore + second["other"]
    return dict(first=first, second=second,
                intercept_heater=max(0.0, -q_int), engine_heater=0.0 if on1 else max(0.0, -q_eng), warm_supply=-heat["warm"])


def engine_on_stage(P):
    """True when the riser engine is held by the first stage, at its temperature, instead of a heater."""
    return P["mode"] == "loop" and P["engine"] and P["engine_stage"]


# PT420-RM reference capacity curve (Cryomech, 50 and 60 Hz), read off the published chart:
# stage temperatures [K] for first-stage loads PT420_Q1 (columns) and second-stage loads PT420_Q2 (rows) [W].
PT420_Q1 = [0.0, 25.0, 50.0, 75.0, 100.0]
PT420_Q2 = [0.0, 0.5, 1.8, 5.0, 8.0, 10.0]
PT420_T1 = [[28.6, 33.5, 40.2, 49.4, 62.7],
            [29.0, 34.2, 40.9, 50.3, 65.5],
            [29.7, 36.1, 42.9, 52.4, 68.0],
            [32.5, 39.3, 47.7, 60.7, 71.2],
            [35.2, 42.9, 54.5, 64.3, 74.0],
            [37.1, 46.0, 56.1, 65.1, 76.7]]
PT420_T2 = [[2.31, 2.28, 2.34, 2.63, 2.70],
            [3.02, 2.99, 2.99, 3.11, 2.89],
            [4.43, 4.34, 4.17, 4.24, 4.08],
            [7.35, 6.87, 6.52, 6.81, 7.39],
            [9.76, 9.44, 9.51, 9.63, 9.76],
            [10.69, 10.63, 10.73, 10.60, 11.53]]


def seg(xs, x):
    """Interval i of the increasing list xs that holds x (the end intervals extend outward), and the fraction along it."""
    i = 0
    while i < len(xs) - 2 and x > xs[i + 1]:
        i += 1
    return i, (x - xs[i]) / (xs[i + 1] - xs[i])


def pt_map(Q1, Q2, sc):
    """Stage temperatures (T1, T2) with loads Q1, Q2 [W] on a cooler at fraction sc of the reference
    capacity: bilinear in the chart grid, extrapolated past its edges."""
    i, u = seg(PT420_Q1, Q1 / sc)
    j, v = seg(PT420_Q2, Q2 / sc)

    def bil(M):
        a = M[j][i] + u * (M[j][i + 1] - M[j][i])
        b = M[j + 1][i] + u * (M[j + 1][i + 1] - M[j + 1][i])
        return a + v * (b - a)
    return bil(PT420_T1), bil(PT420_T2)


def pt_load2(Q1, T2, sc):
    """Second-stage load that puts the second stage at T2 with Q1 on the first stage (negative if out of reach)."""
    i, u = seg(PT420_Q1, Q1 / sc)
    col = [r[i] + u * (r[i + 1] - r[i]) for r in PT420_T2]     # rises with Q2 at every Q1
    j, v = seg(col, T2)
    return sc * (PT420_Q2[j] + v * (PT420_Q2[j + 1] - PT420_Q2[j]))


def pt_load1(Q2, T1, sc):
    """First-stage load that puts the first stage at T1 with Q2 on the second stage."""
    j, v = seg(PT420_Q2, Q2 / sc)
    row = [a + v * (b - a) for a, b in zip(PT420_T1[j], PT420_T1[j + 1])]
    i, u = seg(row, T1)
    return sc * (PT420_Q1[i] + u * (PT420_Q1[i + 1] - PT420_Q1[i]))


# ------------------------------------------------------------------ entry
def thermal(P):
    """Tube profiles, loop flow and clamp heats with the target cell at P["Tc"]. Also the clamp heats of the same
    tubes with no flow (conduction alone): the difference is the heat the circulating gas carries between clamps."""
    p = P["p_mbar"] * 100.0
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
    heat = dict(cold=0.0, intercept=0.0, engine=0.0, warm=0.0)
    for s in sol:
        for kind, Tn, q in s["heats"]:
            heat[kind] += q
    still = dict(heat)
    if mdot != 0.0:
        still = dict(cold=0.0, intercept=0.0, engine=0.0, warm=0.0)
        for leg in legs:
            for kind, Tn, q in solve_leg(leg, tube, p, 0.0)["heats"]:
                still[kind] += q
    return tube, legs, errs, mdot, sol, heat, still


def cooler_solve(P):
    """Stage temperatures and loads on the PT420.
    Map on: the first stage floats on the map; the second stage is held at the set point P["Tc"] by
    a trim heater, or floats above it when the loads alone keep it warmer. The tube solution depends
    on Tc, so it is redone, with secant steps, until Tc settles; the loads change little with Tc.
    Map off: Tc and T_shield are taken as given, and the map gives each stage's capacity there.
    An engine held by the first stage sits at T_shield, so with the map on the tube solution is
    also redone until the first-stage temperature settles.
    Returns (P with the stage and engine temperatures, the tube solution, the stage loads, the cooler state)."""
    sc = P["cooler_pct"] / 100.0
    Tset = P["Tc"]
    on1 = engine_on_stage(P)
    if not P["cooler_map"]:
        if on1:
            P = dict(P, T_engine=P["T_shield"])
        th = thermal(P)
        st = stage_loads(P, th[5], th[0], len(th[1]))
        Q1, Q2 = st["first"]["total"], st["second"]["total"]
        T1f, T2f = pt_map(Q1, Q2, sc)
        st["first"]["capacity"] = pt_load1(Q2, P["T_shield"], sc)
        st["second"]["capacity"] = pt_load2(Q1, Tset, sc)
        cooler = dict(map=False, T1=P["T_shield"], T2=Tset, set=Tset, trim=0.0,
                      held=Q1 <= st["first"]["capacity"] and Q2 <= st["second"]["capacity"],
                      T1_free=T1f, T2_free=T2f, Q1=Q1, Q2=Q2,
                      off_map=Q1 / sc > PT420_Q1[-1] or Q2 / sc > PT420_Q2[-1])
        return P, th, st, cooler
    Tc, T1, prev = Tset, 40.0, None
    for _ in range(30):
        Te = T1 if on1 else P["T_engine"]
        th = thermal(dict(P, Tc=Tc, T_engine=Te))
        for _ in range(100):
            st = stage_loads(dict(P, Tc=Tc, T_shield=T1), th[5], th[0], len(th[1]))
            Q1, Q2 = st["first"]["total"], st["second"]["total"]
            T1f, T2f = pt_map(Q1, Q2, sc)
            cap2 = pt_load2(Q1, Tset, sc)
            trim = cap2 - Q2 if T2f < Tset else 0.0
            T1n, T2 = pt_map(Q1, Q2 + trim, sc) if trim > 0 else (T1f, T2f)
            if abs(T1n - T1) < 1e-9:
                break
            T1 = T1n
        r = T2 - Tc
        if abs(r) < 1e-7 and (not on1 or abs(T1 - Te) < 1e-7):
            break
        Tn = Tc - r * (Tc - prev[0]) / (r - prev[1]) if prev and r != prev[1] else T2   # secant on T2(Tc) - Tc
        prev = (Tc, r)
        Tc = max(Tset, Tn)
    st["first"]["capacity"] = None
    st["second"]["capacity"] = cap2
    cooler = dict(map=True, T1=T1, T2=Tc, set=Tset, trim=trim, held=T2f <= Tset,
                  T1_free=T1f, T2_free=T2f, Q1=Q1, Q2=Q2,
                  off_map=Q1 / sc > PT420_Q1[-1] or (Q2 + trim) / sc > PT420_Q2[-1])
    return dict(P, Tc=Tc, T_shield=T1, T_engine=Te), th, st, cooler


def system(P):
    """Cooler, tubes and exchange: the stage temperatures, G, the atoms in each cell, and the
    warnings these raise. Everything the rate equations take from the transfer model."""
    P, th, stages, cooler = cooler_solve(P)
    tube, legs, errs, mdot, sol, heat, still = th
    p = P["p_mbar"] * 100.0
    Tc, Th = P["Tc"], P["Th"]

    G_diff = sum(1.0 / s["R_diff"] for s in sol)
    G_conv = abs(mdot) / M_ATOM
    G = G_diff + G_conv

    Np = n_num(Th, p) * P["V_pump_cc"] * 1e-6
    Nt = n_num(Tc, p) * P["V_target_cc"] * 1e-6

    warnings = list(errs)
    Re = 4.0 * abs(mdot) / (math.pi * tube.d * eta(Tc))
    lag = gas_lag(sol, legs, mdot, Tc)
    if lag:
        dn = [d for d in lag["legs"] if not d["up"]][0]
        icp = [mk for mk in dn["marks"] if mk[0] == "intercept"]
        if max(abs(d["max_K"]) for d in lag["legs"]) > 10.0 or dn["arrive_K"] - dn["arrive_wall_K"] > max(1.0, 0.2 * Tc):
            warnings.append("The gas lags the tube walls near the cold end: coming down, it %sreaches the target cell at %.3g K (target %.3g K). The model takes the gas at the wall temperature, so it %slikely overestimates the loop flow; see the gas temperature on the Heat loads tab."
                            % ("passes the intercept at %.3g K (wall %.3g K) and " % (icp[0][3], icp[0][2]) if icp else "", dn["arrive_K"], Tc,
                               "puts too much of the heat on the intercept and " if icp else ""))
    if cooler["map"]:
        if not cooler["held"]:
            cap2 = stages["second"]["capacity"]
            if cap2 > 0:
                warnings.append("The PT420 can't hold the target at %.3g K: it gives %.3g W there with %.3g W on the first stage, short of the %.3g W load, so the target runs at %.3g K."
                                % (cooler["set"], cap2, cooler["Q1"], cooler["Q2"], Tc))
            else:
                warnings.append("The PT420 can't reach %.3g K with %.3g W on the first stage, so the target runs at %.3g K."
                                % (cooler["set"], cooler["Q1"], Tc))
    else:
        for name, stage, T_free in (("Second", stages["second"], cooler["T2_free"]), ("First", stages["first"], cooler["T1_free"])):
            if stage["total"] > stage["capacity"]:
                if stage["capacity"] > 0:
                    warnings.append("%s-stage load of %.3g W exceeds the PT420's %.3g W at %.3g K; with these loads it settles at %.3g K."
                                    % (name, stage["total"], stage["capacity"], stage["T"], T_free))
                else:
                    warnings.append("The PT420 can't reach %.3g K on the %s stage with these loads; it settles at %.3g K."
                                    % (stage["T"], name.lower(), T_free))
    if cooler["off_map"]:
        sc = P["cooler_pct"] / 100.0
        warnings.append("The loads run past the edge of the PT420 chart (%.3g W on the first stage, %.3g W on the second), so the stage temperatures are extrapolated."
                        % (PT420_Q1[-1] * sc, PT420_Q2[-1] * sc))
    if P["mode"] == "loop" and not P["engine"] and not P["intercept"]:
        warnings.append("With no engine and no intercept the two legs are identical, so nothing drives the loop; the exchange is by diffusion alone.")
    if P["mode"] == "loop" and mdot < 0:
        warnings.append("The loop runs backwards: the riser ends up colder and denser than the downcomer, so gas flows down the riser and up the downcomer.")
    if any(s["clamped"] for s in sol):
        warnings.append("Flow is strong enough to flatten a free section's temperature profile; treat these results as rough.")
    return dict(P=P, tube=tube, legs=legs, mdot=mdot, sol=sol, heat=heat, still=still, stages=stages, cooler=cooler,
                G=G, G_diff=G_diff, G_conv=G_conv, Np=Np, Nt=Nt, Re=Re, lag=lag, warnings=warnings)


def pumping_effective(P, Np, M_ss, T_b_s, T1p_s):
    """A sealed-cell fit (M_ss as a fraction, T_b, and T1p with the discharge on and the laser off) as fluxes,
    and the effective M_ss, T_b and T1p of a pumping cell of Np atoms that the rate equations take.
    In flux mode the optical production Q0' and the relaxation R' are set by the discharge, with the light taken to
    keep up with it, so both scale by s; in volume mode the fit is used as it stands, so production grows with Np."""
    g_op = 1.0 / T_b_s - 1.0 / T1p_s
    out = dict(mode=P["pump_scaling"], g_op=g_op, P_lim=M_ss / (g_op * T_b_s) if g_op > 0 else None, eta=None)
    if P["pump_scaling"] == "volume":
        out.update(Ns=None, Q0=None, R_dis=None, s=None, Q0p=Np * M_ss / T_b_s, Rp=Np / T1p_s,
                   M_ss_eff=M_ss, T_b_eff=T_b_s, T1p_eff=T1p_s)
        return out
    Ns = n_num(P["Th"], P["p_mbar"] * 100.0) * P["V_seal_cc"] * 1e-6
    Q0, R_dis = Ns * M_ss / T_b_s, Ns / T1p_s
    if P["Pabs_on"]:                    # s from the absorbed pump power, which keeps up with the discharge
        s = P["Pabs_mW"] / P["Pabs_seal_mW"]
        out["eta"] = Q0 / (P["Pabs_seal_mW"] * 1e-3 / HNU)
    else:
        s = P["s_dis"]
    dark = Np / (P["T1_dark_h"] * 3600.0) if P["T1_dark_on"] else 0.0
    Q0p, Rp = s * Q0, s * R_dis + dark
    loss = s * Ns / T_b_s + dark        # Q0'/P_lim + R' = s Ns (Gop + 1/T1p) + dark
    out.update(Ns=Ns, Q0=Q0, R_dis=R_dis, s=s, Q0p=Q0p, Rp=Rp,
               M_ss_eff=Q0p / loss, T_b_eff=Np / loss, T1p_eff=Np / Rp)
    return out


def pumping_warnings(P, pe):
    """What pumping_effective()'s result says about the inputs, in flux mode."""
    if pe["mode"] == "volume":
        return []
    out = []
    if pe["eta"] is not None and pe["eta"] > 2.0:
        out.append("The sealed-cell fit gives %.3g polarized atoms per absorbed photon, N_s M_ss/T_b over P_abs/(h nu), more than the 2 a photon can carry: the absorbed power or the sealed-cell volume is off."
                   % pe["eta"])
    if pe["T_b_eff"] > P["on_h"] * 3600.0:
        out.append("The pumping cell's build-up time T_b' = %.3g h is longer than the %.3g h run: the discharge can't fill a pumping cell this large in that time."
                   % (pe["T_b_eff"] / 3600.0, P["on_h"]))
    return out


def compute(P):
    S = system(P)
    P, legs, mdot, sol, heat, still, stages, cooler = S["P"], S["legs"], S["mdot"], S["sol"], S["heat"], S["still"], S["stages"], S["cooler"]
    Tc = P["Tc"]
    G, G_diff, G_conv, Np, Nt, Re, lag = S["G"], S["G_diff"], S["G_conv"], S["Np"], S["Nt"], S["Re"], S["lag"]

    pe = pumping_effective(P, Np, P["M_ss"] / 100.0, P["T_b_s"], P["tau_p_s"])
    warnings = S["warnings"] + pumping_warnings(P, pe)
    if pe["g_op"] <= 0:
        warnings.append("tau_p (%s s) must be longer than the sealed-cell build-up time T_b (%s s), since 1/T_b = Gop + 1/tau_p."
                        % (fmt4(P["tau_p_s"]), fmt4(P["T_b_s"])))
    Mss = pe["M_ss_eff"]
    g_b = 1.0 / pe["T_b_eff"]
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
    sys2 = two_cell(Np, Nt, G, pe["T1p_eff"], P["T1t_h"] * 3600.0)
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

    return dict(
        ok=True,
        engine="python",
        mdot_mg_h=mdot * 3.6e9,
        Re=Re,
        G=G, G_diff=G_diff, G_conv=G_conv,
        Np=Np, Nt=Nt, N_tubes=sum(s["N_tube"] for s in sol),
        tau_x_target_s=Nt / G, tau_x_pump_s=Np / G,
        heat_mW={k: v * 1e3 for k, v in heat.items()},
        heat_still_mW={k: v * 1e3 for k, v in still.items()},
        stages=stages,
        Tc=Tc, T_shield=P["T_shield"], T_engine=P["T_engine"], cooler=cooler,
        legs=[dict(profile=s["profile"], length_cm=s["length"] * 100,
                   nodes=[[n[0], n[1] * 100, n[2]] for n in leg["nodes"]],
                   free_cm=[f * 100 for f in leg["free"]],
                   G_diff=1.0 / s["R_diff"])
              for s, leg in zip(sol, legs)],
        ss_pump=ss_on[0], ss_target=ss_on[1],
        times_s=eig_times(A_on),
        t90_target_s=t90,
        production_atoms_s=pe["Q0p"], pump_eff=pe,
        loss_exchange=(Nt / G) * g_t, loss_capacity=(Nt / Np) * pe["T_b_eff"] * g_t,
        t_h=[x / 3600.0 for x in t], P_pump=pp, P_target=pt,
        paper=paper,
        lag=lag,
        zero_t_s=tz, zero_P_pump=zp, zero_P_target=zt, shut_P_pump=sp, shut_P_target=st,
        warnings=warnings,
    )


def run(params_json):
    try:
        return json.dumps(compute(json.loads(params_json)))
    except Exception as e:  # report, don't crash the page
        return json.dumps(dict(ok=False, error="%s: %s" % (type(e).__name__, e)))


# ------------------------------------------------------------------ discharge schedule
# Each discharge level d is a sealed-cell fit (M_ss,d, T_b,d, T1p,d); pumping_effective() turns it into the
# pumping cell's (M_ss,d', T_b,d'), so its net production per pumping-cell atom is the straight line
# f_d(P) = (M_ss,d' - P)/T_b,d' = a_d - s_d P. In volume mode T1p only splits it: f_d = Gop (Plim - P) - P/T1p
# with Gop = 1/T_b - 1/T1p and Gop Plim = M_ss/T_b, so it leaves f_d, and the dynamics, unchanged.
# Within a level the rate equations are those of compute(), solved exactly; levels change at
# events, each located by scanning the exact solution and bisecting the first bracket.

SYSTEM_KEYS = ("mode", "p_mbar", "Tc", "Th", "id_mm", "wall_mm", "length_cm",
               "engine", "engine_stage", "engine_cm", "gap_cm", "T_engine",
               "intercept", "T_int", "s_int_cm", "V_pump_cc", "V_target_cc",
               "T_shield", "A_shield_m2", "eps_shield", "other1_W",
               "A_cold_cm2", "eps_cold", "other2_mW", "cooler_map", "cooler_pct")
SCAN_N = 2000            # steps across the horizon when looking for an event
HORIZON_TAU = 20.0       # the horizon covers this many slow time constants of the slowest level
TAIL_TAU = 40.0          # a level left on for good is searched this many slow time constants out
MAX_SWITCHES = 200
_sys_memo = [None, None]


def system_cached(P):
    """system() for the inputs that set it; the schedule page reruns often with only the levels changed."""
    sp = {k: P[k] for k in SYSTEM_KEYS}
    key = json.dumps(sp, sort_keys=True)
    if _sys_memo[0] != key:
        _sys_memo[1] = system(sp)
        _sys_memo[0] = key
    return _sys_memo[1]


def advance(A, xs, x0, dt):
    """Exact solution of dx/dt = A (x - xs), from x0, after dt."""
    E = expm2(A, dt)
    d0 = [x0[0] - xs[0], x0[1] - xs[1]]
    return [xs[0] + E[0][0] * d0[0] + E[0][1] * d0[1], xs[1] + E[1][0] * d0[0] + E[1][1] * d0[1]]


def first_time(g, t0, t1, dt):
    """First time in [t0, t1] where g(t) > 0: scan in steps of about dt, then bisect the first bracket."""
    if g(t0) > 0:
        return t0
    n = max(1, int(math.ceil((t1 - t0) / dt)))
    lo = t0
    for k in range(1, n + 1):
        hi = t0 + (t1 - t0) * k / n
        if g(hi) > 0:
            for _ in range(200):
                mid = 0.5 * (lo + hi)
                if mid <= lo or mid >= hi:
                    break
                if g(mid) > 0:
                    hi = mid
                else:
                    lo = mid
            return hi
        lo = hi
    return None


def fmt4(v):
    """A number to four decimals, trailing zeros dropped (written the same way by the JS port)."""
    return ("%.4f" % v).rstrip("0").rstrip(".")


def level_rates(G, Np, Nt, g_t, pe):
    """A and b of compute()'s rate equations with a discharge level on; pe is its pumping_effective()."""
    g_b = 1.0 / pe["T_b_eff"]
    kp, kt = G / Np, G / Nt
    return [[-(g_b + kp), kp], [kt, -(kt + g_t)]], [g_b * pe["M_ss_eff"], 0.0]


def crossover(a, s, i, j):
    """Polarization where f_i and f_j cross, or None for parallel lines."""
    if s[i] == s[j]:
        return None
    return (a[i] - a[j]) / (s[i] - s[j])


def envelope(a, s, P_max):
    """Upper envelope of the lines f_d = a_d - s_d P over [0, P_max]: [[P_from, P_to, level], ...]."""
    n = len(a)
    c = 0
    for d in range(1, n):
        if a[d] > a[c] or (a[d] == a[c] and s[d] < s[c]):
            c = d
    pieces, P0 = [], 0.0
    while True:
        best, nd = None, None
        for d in range(n):
            if d == c or s[d] >= s[c]:
                continue
            x = crossover(a, s, c, d)
            if x > P0 and (best is None or x < best or (x == best and s[d] < s[nd])):
                best, nd = x, d
        if best is None or best >= P_max:
            pieces.append([P0, P_max, c])
            return pieces
        pieces.append([P0, best, c])
        P0, c = best, nd


def schedule(Np, Nt, G, T1t_s, levels, S):
    """Build-up from zero under a discharge schedule. levels: [{label, M_ss %, T_b_s, T1p_s}], sealed-cell fits;
    S: sched_mode ("single", "switch" or "envelope") and its settings, the goal, on_h and the pumping-cell scaling."""
    if not levels:
        raise ValueError("no discharge levels")
    n = len(levels)
    g_t = 1.0 / T1t_s
    N = Np + Nt
    alpha, beta = Np / N, (Nt / N) * g_t
    E = [pumping_effective(S, Np, lv["M_ss"] / 100.0, lv["T_b_s"], lv["T1p_s"]) for lv in levels]
    L = []
    for pe in E:
        A, b = level_rates(G, Np, Nt, g_t, pe)
        L.append(dict(A=A, xs=ss(A, b), times=eig_times(A)))
    a = [level_rates(G, Np, Nt, g_t, pe)[1][0] for pe in E]           # f_d(0) = M_ss'/T_b'
    s = [1.0 / pe["T_b_eff"] for pe in E]                             # -df_d/dP = 1/T_b'
    span = S["on_h"] * 3600.0
    T_h = max(span, HORIZON_TAU * max(l["times"][0] for l in L))
    dt = T_h / SCAN_N
    warnings = []

    def state(segs, t):
        k = len(segs) - 1
        while k > 0 and segs[k]["t0"] > t:
            k -= 1
        sg = segs[k]
        l = L[sg["level"]]
        return advance(l["A"], l["xs"], sg["x0"], t - sg["t0"])

    def seg(t0, d, x0):
        return dict(t0=t0, level=d, x0=x0)

    def pp_on(sg):
        l = L[sg["level"]]
        return lambda t: advance(l["A"], l["xs"], sg["x0"], t - sg["t0"])

    def goal_time(segs, goal):
        """First time the target reaches the goal, or None."""
        for k, sg in enumerate(segs):
            t1 = segs[k + 1]["t0"] if k + 1 < len(segs) else sg["t0"] + TAIL_TAU * L[sg["level"]]["times"][0]
            x_at = pp_on(sg)
            t = first_time(lambda tt: x_at(tt)[1] - goal, sg["t0"], t1, dt)
            if t is not None:
                return t
        return None

    # the schedule
    mode = S["sched_mode"]
    P_sw = None
    clamp = lambda i: min(n - 1, max(0, int(i)))
    if mode == "single":
        segs = [seg(0.0, clamp(S["single_level"]), [0.0, 0.0])]
    elif mode == "switch":
        e, l = clamp(S["early_level"]), clamp(S["late_level"])
        segs = [seg(0.0, e, [0.0, 0.0])]
        t_sw = None
        if e != l:
            if S["trigger"] == "t":
                t_sw = S["t_switch_h"] * 3600.0
            else:
                if S["trigger"] == "P":
                    P_sw = S["P_switch"] / 100.0
                elif a[e] - a[l] <= 0:
                    P_sw = 0.0                          # the late level is better from the start
                elif s[e] > s[l]:
                    P_sw = crossover(a, s, e, l)
                if P_sw is not None:
                    x_at = pp_on(segs[0])
                    t_sw = first_time(lambda t: x_at(t)[0] - P_sw, 0.0, T_h, dt) if P_sw > 0 else 0.0
        if t_sw is not None:
            segs = [seg(0.0, l, [0.0, 0.0])] if t_sw <= 0 else segs + [seg(t_sw, l, state(segs, t_sw))]
    elif mode == "envelope":
        h = S["hyst_pct"] / 100.0
        dwell = S["dwell_min"] * 60.0
        c = envelope(a, s, 1.0)[0][2]
        segs = [seg(0.0, c, [0.0, 0.0])]

        def gain(c, Pp):
            """Best margin of another level over c at Pp, less the hysteresis, and that level."""
            best, bd = None, None
            for d in range(n):
                if d == c:
                    continue
                v = (a[d] - s[d] * Pp) - (a[c] - s[c] * Pp) - h * abs(s[d] - s[c])
                if best is None or v > best:
                    best, bd = v, d
            return best, bd
        while n > 1:
            sg = segs[-1]
            t_a = sg["t0"] + dwell
            if t_a >= T_h:
                break
            x_at = pp_on(sg)
            t_sw = first_time(lambda t: gain(sg["level"], x_at(t)[0])[0], t_a, T_h, dt)
            if t_sw is None:
                break
            if len(segs) > MAX_SWITCHES:
                warnings.append("The envelope schedule switched more than %d times; add hysteresis or a minimum dwell." % MAX_SWITCHES)
                break
            x = x_at(t_sw)
            segs.append(seg(t_sw, gain(sg["level"], x[0])[1], x))
    else:
        raise ValueError("unknown schedule mode %r" % mode)

    # each level on its own, and the goal
    best_ss = max(l["xs"][1] for l in L)
    goal = S["goal_frac"] / 100.0 * best_ss if S["goal_mode"] == "frac" else S["goal_pct"] / 100.0
    t_grid = [span * i / 360 for i in range(361)]
    out_levels, single_p, single_t = [], [], []
    for d, (lv, l, pe) in enumerate(zip(levels, L, E)):
        one = [seg(0.0, d, [0.0, 0.0])]
        xs = [state(one, t) for t in t_grid]
        single_p.append([x[0] for x in xs])
        single_t.append([x[1] for x in xs])
        g_op = pe["g_op"]
        k_sys = alpha * s[d] + beta
        out_levels.append(dict(
            label=lv["label"], M_ss=lv["M_ss"], T_b_s=lv["T_b_s"], T1p_s=lv["T1p_s"], pump_eff=pe,
            f0=a[d], slope=s[d], gamma_op=g_op, P_lim=pe["P_lim"],
            x_lim=Np * s[d] / G, k_sys=k_sys, P_sys=alpha * a[d] / k_sys,
            ss_pump=l["xs"][0], ss_target=l["xs"][1], times_s=l["times"],
            t_goal_s=goal_time(one, goal)))
        if g_op <= 0:
            warnings.append("%s: T1p (%s s) must be longer than the sealed-cell build-up time T_b (%s s), since 1/T_b = Gop + 1/T1p."
                            % (lv["label"], fmt4(lv["T1p_s"]), fmt4(lv["T_b_s"])))
        warnings.extend("%s: %s" % (lv["label"], w) for w in pumping_warnings(S, pe))
        if Np * s[d] / G > 0.2:
            warnings.append("%s: N_p/(T_b G) = %.2f. The exchange %s the build-up, so the pumping cell runs ahead of the target "
                            "and a stronger discharge gains less than its f_d suggests."
                            % (lv["label"], Np * s[d] / G, "limits" if Np * s[d] / G > 1 else "is close to limiting"))

    # the schedule's trace, with each switch on the grid
    switches = [{"t_s": sg["t0"], "from": segs[k]["level"], "to": sg["level"], "P_pump": sg["x0"][0], "P_target": sg["x0"][1]}
                for k, sg in enumerate(segs[1:])]
    tt = sorted(set(t_grid) | set(w["t_s"] for w in switches if 0 < w["t_s"] < span))
    xx = [state(segs, t) for t in tt]
    t_goal = goal_time(segs, goal)
    fin = L[segs[-1]["level"]]["xs"]
    singles = [x["t_goal_s"] for x in out_levels if x["t_goal_s"] is not None]
    best_t = min(singles) if singles else None
    if goal > best_ss:
        warnings.append("No level on its own reaches the goal of %.1f%%: the best steady state is %.1f%%." % (goal * 100, best_ss * 100))
    elif t_goal is None:
        warnings.append("The schedule never reaches the goal of %.1f%%; it settles at %.1f%%." % (goal * 100, fin[1] * 100))

    P_max = max(pe["M_ss_eff"] for pe in E)
    pieces = envelope(a, s, P_max)
    kinks = [(pieces[k][1], set((pieces[k][2], pieces[k + 1][2]))) for k in range(len(pieces) - 1)]
    crossings = []
    for i in range(n):
        for j in range(i + 1, n):
            x = crossover(a, s, i, j)
            if x is not None and 0 < x < P_max:
                crossings.append(dict(P=x, i=i, j=j, f=a[i] - s[i] * x,
                                      on_envelope=any(xk == x and pair == set((i, j)) for xk, pair in kinks)))

    return dict(
        levels=out_levels, goal=goal, best_ss_target=best_ss, horizon_s=T_h,
        loss_slope=(Nt / Np) * g_t,
        t_h=[t / 3600.0 for t in t_grid], single_P_pump=single_p, single_P_target=single_t,
        sched=dict(mode=mode, t_h=[t / 3600.0 for t in tt], P_pump=[x[0] for x in xx], P_target=[x[1] for x in xx],
                   segments=[dict(t0_s=sg["t0"], t1_s=segs[k + 1]["t0"] if k + 1 < len(segs) else None, level=sg["level"])
                             for k, sg in enumerate(segs)],
                   switches=switches, P_switch=P_sw, t_goal_s=t_goal, final_level=segs[-1]["level"],
                   ss_pump=fin[0], ss_target=fin[1], best_single_s=best_t,
                   saved=(best_t - t_goal) / best_t if best_t and t_goal is not None else None),
        envelope=dict(P_max=P_max, pieces=pieces, crossings=crossings),
        warnings=warnings,
    )


def compute_schedule(P):
    S = system_cached(P)
    out = schedule(S["Np"], S["Nt"], S["G"], P["T1t_h"] * 3600.0, P["levels"], P)
    out.update(ok=True, engine="python", G=S["G"], G_diff=S["G_diff"], G_conv=S["G_conv"], Np=S["Np"], Nt=S["Nt"],
               Tc=S["P"]["Tc"], T_shield=S["P"]["T_shield"], mdot_mg_h=S["mdot"] * 3.6e9,
               tau_x_target_s=S["Nt"] / S["G"], tau_x_pump_s=S["Np"] / S["G"],
               warnings=S["warnings"] + out["warnings"])
    return out


def run_schedule(params_json):
    try:
        return json.dumps(compute_schedule(json.loads(params_json)))
    except Exception as e:  # report, don't crash the page
        return json.dumps(dict(ok=False, error="%s: %s" % (type(e).__name__, e)))


DEFAULTS = dict(
    mode="loop", p_mbar=100.0, Tc=5.0, Th=300.0,
    id_mm=8.0, wall_mm=1.0, length_cm=8.0,
    engine=True, engine_stage=False, engine_cm=2.0, gap_cm=0.5, T_engine=25.0,
    intercept=True, T_int=45.0, s_int_cm=4.0,
    V_pump_cc=49.0, V_target_cc=100.0,
    M_ss=50.0, T_b_s=150.0, tau_p_s=600.0, T1t_h=10.0,
    pump_scaling="flux", V_seal_cc=49.0, s_dis=1.0, T1_dark_on=False, T1_dark_h=10.0,
    Pabs_on=False, Pabs_seal_mW=200.0, Pabs_mW=200.0,
    on_h=24.0,
    T_shield=45.0, A_shield_m2=0.5, eps_shield=0.003, other1_W=0.0,
    A_cold_cm2=300.0, eps_cold=0.05, other2_mW=0.0,
    cooler_map=True, cooler_pct=100.0,
)

# Placeholder discharge levels, to be replaced with measured sealed-cell values: "Low" is the main
# page's default pumping, "High" a stronger discharge that pumps faster to a lower polarization.
SCHEDULE_DEFAULTS = dict(
    levels=[dict(label="Low", M_ss=50.0, T_b_s=150.0, T1p_s=600.0),
            dict(label="High", M_ss=40.0, T_b_s=60.0, T1p_s=200.0)],
    sched_mode="switch", single_level=0, early_level=1, late_level=0,
    trigger="auto", P_switch=20.0, t_switch_h=1.0, hyst_pct=0.0, dwell_min=0.0,
    goal_mode="frac", goal_pct=30.0, goal_frac=90.0,
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
        print("  stages", r["stages"])
        print("  cooler", r["cooler"])
    for sm in ("single", "switch", "envelope"):
        r = compute_schedule(dict(DEFAULTS, **dict(SCHEDULE_DEFAULTS, sched_mode=sm)))
        print("schedule", sm, {k: r["sched"][k] for k in ("t_goal_s", "best_single_s", "saved", "P_switch", "switches")},
              "goal", r["goal"], [(x["label"], x["t_goal_s"], x["ss_target"]) for x in r["levels"]])
