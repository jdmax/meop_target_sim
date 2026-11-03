"""
Checks of the discharge-schedule model in he3_model.py, and of its JavaScript port.

Run from the repository root with either
    python tests/test_schedule.py
    python -m pytest tests/test_schedule.py
The parity checks need Node.js on the PATH and are skipped without it.
"""
import json
import math
import os
import shutil
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import he3_model as m  # noqa: E402

H = 3600.0


# ------------------------------------------------------------------ reduced model
# In the MEOP-limited limit (N_p/T_b << G) the cells share one polarization P, and
#   N dP/dt = N_p (M_ss - P)/T_b - N_t P/T_1t,   i.e.   dP/dt = k (P_inf - P)
# with k = alpha/T_b + beta and k P_inf = alpha M_ss/T_b, where alpha = N_p/N and beta = (N_t/N)/T_1t.
# These build sealed-cell levels that give a chosen whole-system pair (k, P_inf).
NP, NT = 1.0, 120.0
ALPHA = NP / (NP + NT)
BETA = 0.1 / H                          # loss rate shared by every level
T1T = (NT / (NP + NT)) / BETA


def level(label, k_h, P_inf):
    k = k_h / H
    return dict(label=label, M_ss=100.0 * k * P_inf / (k - BETA), T_b_s=ALPHA / (k - BETA), T1p_s=1e9)


def reduced_time(k_h, P_inf, P0, P1):
    """Hours for dP/dt = k (P_inf - P) to go from P0 to P1."""
    return math.log((P_inf - P0) / (P_inf - P1)) / k_h


def run_reduced(levels, goal, xlim=1e-6, **S):
    """The full model with G large enough that N_p/(T_b G) = xlim for the fastest level."""
    G = NP * max(1.0 / lv["T_b_s"] for lv in levels) / xlim
    sched = dict(m.SCHEDULE_DEFAULTS, on_h=24.0, goal_mode="abs", goal_pct=100.0 * goal, early_level=1, late_level=0,
                 pump_scaling="volume")      # the levels are whole-system rates for NP abstract atoms
    sched.update(S)
    return m.schedule(NP, NT, G, T1T, levels, sched)


CASES = [
    # (k_L, P_L), (k_H, P_H), goal, P*, low-only h, switched h (the hand calculations)
    ((0.3, 0.40), (0.5, 0.32), 0.35, 0.20, 6.93, 6.58),
    ((0.3, 0.40), (1.1, 0.364), 0.38, 0.35, 10.0, 6.0),
]


def test_reduced_model():
    for (kL, PL), (kH, PH), goal, P_star, t_low, t_sw in CASES:
        levels = [level("Low", kL, PL), level("High", kH, PH)]
        r = run_reduced(levels, goal)
        # the crossover from the f_d lines is the reduced-model formula, whatever the losses
        P_red = (kH * PH - kL * PL) / (kH - kL)
        assert abs(r["sched"]["P_switch"] - P_red) < 1e-12, (r["sched"]["P_switch"], P_red)
        assert round(P_red, 2) == P_star
        # whole-system rates and limits the page reports
        assert abs(r["levels"][0]["k_sys"] * H - kL) < 1e-12 and abs(r["levels"][0]["P_sys"] - PL) < 1e-12
        assert abs(r["levels"][1]["k_sys"] * H - kH) < 1e-12 and abs(r["levels"][1]["P_sys"] - PH) < 1e-12
        # times: exact reduced model, and the hand values to the digits quoted
        exact_low = reduced_time(kL, PL, 0.0, goal)
        exact_sw = reduced_time(kH, PH, 0.0, P_red) + reduced_time(kL, PL, P_red, goal)
        got_low = r["levels"][0]["t_goal_s"] / H
        got_sw = r["sched"]["t_goal_s"] / H
        assert abs(got_low / exact_low - 1) < 1e-4, (got_low, exact_low)
        assert abs(got_sw / exact_sw - 1) < 1e-4, (got_sw, exact_sw)
        assert abs(got_low - t_low) < 0.05 and abs(got_sw - t_sw) < 0.05, (got_low, got_sw)
        assert len(r["sched"]["switches"]) == 1 and abs(r["sched"]["switches"][0]["P_pump"] - P_red) < 1e-12
        # the times converge on the reduced model as the exchange gets faster
        errs = []
        for xlim in (1e-2, 1e-3, 1e-4):
            rr = run_reduced(levels, goal, xlim=xlim)
            errs.append(abs(rr["sched"]["t_goal_s"] / H / exact_sw - 1))
        assert errs[0] > errs[1] > errs[2], errs
        print("reduced  P* %.4f  low %.3f h (exact %.3f)  switched %.3f h (exact %.3f)  saved %.1f%%"
              % (P_red, got_low, exact_low, got_sw, exact_sw, 100 * (got_low - got_sw) / got_low))


def test_envelope_matches_switch():
    """With two levels, no hysteresis and no dwell, following the envelope is the automatic switch."""
    for (kL, PL), (kH, PH), goal, *_ in CASES:
        levels = [level("Low", kL, PL), level("High", kH, PH)]
        a = run_reduced(levels, goal, xlim=1e-3, sched_mode="switch")
        b = run_reduced(levels, goal, xlim=1e-3, sched_mode="envelope")
        assert abs(a["sched"]["t_goal_s"] / b["sched"]["t_goal_s"] - 1) < 1e-12
        assert abs(a["sched"]["switches"][0]["t_s"] / b["sched"]["switches"][0]["t_s"] - 1) < 1e-12


def test_envelope_options():
    """Hysteresis delays the switch past the crossover; the dwell holds the first level at least that long."""
    (kL, PL), (kH, PH), goal, *_ = CASES[0]
    levels = [level("Low", kL, PL), level("High", kH, PH)]
    P_red = (kH * PH - kL * PL) / (kH - kL)
    r = run_reduced(levels, goal, xlim=1e-3, sched_mode="envelope", hyst_pct=2.0)
    assert abs(r["sched"]["switches"][0]["P_pump"] - (P_red + 0.02)) < 1e-12
    t0 = run_reduced(levels, goal, xlim=1e-3, sched_mode="envelope")["sched"]["switches"][0]["t_s"]
    r = run_reduced(levels, goal, xlim=1e-3, sched_mode="envelope", dwell_min=3 * t0 / 60)
    assert abs(r["sched"]["switches"][0]["t_s"] - 3 * t0) < 1e-6 * t0


def test_single_level_reproduces_main_page():
    """One level with the main page's M_ss and T_b gives compute()'s build-up trace exactly."""
    for over in ({}, dict(mode="single"), dict(T1t_h=48.0, M_ss=70.0, T_b_s=40.0, on_h=6.0)):
        P = dict(m.DEFAULTS, **over)
        main = m.compute(dict(P))
        lv = dict(label="Main", M_ss=P["M_ss"], T_b_s=P["T_b_s"], T1p_s=P["tau_p_s"])
        sp = dict(P, **dict(m.SCHEDULE_DEFAULTS, levels=[lv], sched_mode="single", single_level=0, goal_mode="frac", goal_frac=90.0))
        r = m.compute_schedule(sp)
        assert r["sched"]["t_h"] == main["t_h"]
        assert r["sched"]["P_pump"] == main["P_pump"] and r["sched"]["P_target"] == main["P_target"]
        assert r["single_P_target"][0] == main["P_target"]
        assert r["levels"][0]["ss_target"] == main["ss_target"] and r["levels"][0]["ss_pump"] == main["ss_pump"]
        # goal = 90% of this level's steady state is the main page's 90% rise time
        assert abs(r["sched"]["t_goal_s"] / main["t90_target_s"] - 1) < 1e-12, (r["sched"]["t_goal_s"], main["t90_target_s"])
        assert r["G"] == main["G"] and r["Nt"] == main["Nt"]


# ------------------------------------------------------------------ Python / JavaScript parity
PARITY_CASES = [
    dict(),
    dict(sched_mode="single", single_level=1),
    dict(sched_mode="switch", trigger="P", P_switch=25.0),
    dict(sched_mode="switch", trigger="t", t_switch_h=1.5),
    dict(sched_mode="envelope"),
    dict(sched_mode="envelope", hyst_pct=1.0, dwell_min=20.0, goal_mode="abs", goal_pct=31.0),
    dict(levels=[dict(label="Low", M_ss=60.0, T_b_s=300.0, T1p_s=900.0), dict(label="Mid", M_ss=50.0, T_b_s=120.0, T1p_s=400.0),
                 dict(label="High", M_ss=35.0, T_b_s=40.0, T1p_s=150.0)],
         sched_mode="envelope", T1t_h=30.0),
    dict(levels=[dict(label="Low", M_ss=60.0, T_b_s=300.0, T1p_s=900.0), dict(label="Mid", M_ss=50.0, T_b_s=120.0, T1p_s=400.0),
                 dict(label="High", M_ss=35.0, T_b_s=40.0, T1p_s=30.0)],
         sched_mode="envelope", hyst_pct=1.5, dwell_min=45.0, T1t_h=30.0),
    dict(levels=[dict(label="Low", M_ss=60.0, T_b_s=300.0, T1p_s=900.0), dict(label="High", M_ss=35.0, T_b_s=20.0, T1p_s=60.0)],
         mode="single", id_mm=3.0, T1t_h=40.0),
]


def close(a, b, path, tol):
    if isinstance(a, dict):
        assert isinstance(b, dict) and set(a) == set(b), (path, sorted(a), sorted(b))
        for k in a:
            if k != "engine":
                close(a[k], b[k], path + "." + k, tol)
    elif isinstance(a, list):
        assert isinstance(b, list) and len(a) == len(b), (path, len(a), len(b))
        for i, (x, y) in enumerate(zip(a, b)):
            close(x, y, "%s[%d]" % (path, i), tol)
    elif isinstance(a, bool) or a is None or isinstance(a, str):
        assert a == b, (path, a, b)
    else:
        assert abs(a - b) <= tol * max(abs(a), abs(b), 1e-300) or abs(a - b) < 1e-300, (path, a, b)


def run_js(params, fn="run_schedule"):
    node = shutil.which("node")
    if not node:
        return None
    js = ("const M=require(process.argv[1]);let s='';process.stdin.on('data',d=>s+=d)"
          ".on('end',()=>process.stdout.write(M.%s(s)));" % fn)
    p = subprocess.run([node, "-e", js, os.path.join(ROOT, "he3_model.js")], input=json.dumps(params),
                       capture_output=True, text=True, check=True)
    return json.loads(p.stdout)


def test_js_parity():
    if not shutil.which("node"):
        print("parity   skipped: node is not on the PATH")
        return
    worst = 0.0
    for over in PARITY_CASES:
        params = dict(m.DEFAULTS, **dict(m.SCHEDULE_DEFAULTS, **over))
        py = json.loads(m.run_schedule(json.dumps(params)))
        js = run_js(params)
        assert py["ok"] and js["ok"], (py.get("error"), js.get("error"))
        close(py, js, "result", 1e-12)
        # the worst relative difference, for the report
        stack = [(py, js)]
        while stack:
            a, b = stack.pop()
            if isinstance(a, dict):
                stack.extend((a[k], b[k]) for k in a if k != "engine")
            elif isinstance(a, list):
                stack.extend(zip(a, b))
            elif isinstance(a, float) and a != b:
                worst = max(worst, abs(a - b) / max(abs(a), abs(b)))
    print("parity   %d cases, every number agrees; worst relative difference %.2g" % (len(PARITY_CASES), worst))


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("ok      ", name)
