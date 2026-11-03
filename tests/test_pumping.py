"""
Checks of the pumping term in he3_model.py: the sealed-cell fit read as fluxes (pump_scaling = "flux"),
against the older reading as rates (pump_scaling = "volume"), and of the JavaScript port.

Run from the repository root with either
    python tests/test_pumping.py
    python -m pytest tests/test_pumping.py
The parity check needs Node.js on the PATH and is skipped without it.
"""
import json
import os
import shutil
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import he3_model as m  # noqa: E402
from test_schedule import close, run_js  # noqa: E402

FLUX = dict(m.DEFAULTS, pump_scaling="flux")
VOLUME = dict(m.DEFAULTS, pump_scaling="volume")
SCHED = dict(m.SCHEDULE_DEFAULTS, sched_mode="envelope")


def strip(r):
    """A result without what only one of the two modes reports."""
    r = dict(r)
    r.pop("pump_eff", None)
    if "levels" in r:
        r["levels"] = [{k: v for k, v in lv.items() if k != "pump_eff"} for lv in r["levels"]]
    return r


def test_flux_reproduces_volume_mode():
    """With the sealed cell the size of the pumping cell, both readings of the fit agree; in volume mode the
    model is the flux model with V_seal = V_pump, whatever V_pump is."""
    for over in ({}, dict(mode="single"), dict(T1t_h=48.0, M_ss=70.0, T_b_s=40.0, on_h=6.0)):
        close(strip(m.compute(dict(VOLUME, **over))), strip(m.compute(dict(FLUX, **over))), "compute", 1e-9)
        close(strip(m.compute_schedule(dict(VOLUME, **dict(SCHED, **over)))),
              strip(m.compute_schedule(dict(FLUX, **dict(SCHED, **over)))), "schedule", 1e-9)
    for V in (20.0, 200.0):
        vol = m.compute(dict(VOLUME, V_pump_cc=V))
        close(strip(vol), strip(m.compute(dict(FLUX, V_pump_cc=V, V_seal_cc=V))), "compute V=%g" % V, 1e-9)
        assert abs(vol["production_atoms_s"] / (vol["Np"] * m.DEFAULTS["M_ss"] / 100.0 / m.DEFAULTS["T_b_s"]) - 1) < 1e-12
        assert vol["pump_eff"]["T_b_eff"] == m.DEFAULTS["T_b_s"] and vol["pump_eff"]["T1p_eff"] == m.DEFAULTS["tau_p_s"]


def test_volume_independence():
    """Without dark relaxation the steady state and the production don't depend on the pumping cell's size;
    only the dynamics do."""
    base = m.compute(dict(FLUX))
    t90 = []
    for V in (20.0, 49.0, 100.0, 400.0):
        r = m.compute(dict(FLUX, V_pump_cc=V))
        for k in ("ss_target", "ss_pump", "production_atoms_s", "loss_capacity"):
            assert abs(r[k] / base[k] - 1) < 1e-9, (V, k, r[k], base[k])
        assert abs(r["pump_eff"]["M_ss_eff"] - m.DEFAULTS["M_ss"] / 100.0) < 1e-12
        t90.append(r["t90_target_s"])
    assert all(b > a for a, b in zip(t90, t90[1:])), t90        # a bigger cell takes longer to fill
    # in volume mode a bigger pumping cell pumps harder
    assert m.compute(dict(VOLUME, V_pump_cc=400.0))["ss_target"] > base["ss_target"] * 1.3


def test_scaling():
    Np = m.compute(dict(FLUX))["Np"]
    eff = lambda **kw: m.pumping_effective(dict(FLUX, **kw), Np, 0.5, 150.0, 600.0)
    one = eff()
    # a bigger discharge, with the light keeping up, pumps faster to the same M_ss
    for k in (0.25, 3.0, 8.0):
        e = eff(s_dis=k)
        assert abs(e["Q0p"] / (k * one["Q0p"]) - 1) < 1e-12
        assert abs(e["T_b_eff"] * k / one["T_b_eff"] - 1) < 1e-12
        assert abs(e["T1p_eff"] * k / one["T1p_eff"] - 1) < 1e-12        # T1p' is the discharge relaxation alone
        assert abs(e["M_ss_eff"] / 0.5 - 1) < 1e-12
    # T1p < T_b leaves M_ss and T_b as measured
    e = m.pumping_effective(dict(FLUX, s_dis=3.0), Np, 0.5, 150.0, 100.0)
    assert e["P_lim"] is None and abs(e["M_ss_eff"] / 0.5 - 1) < 1e-12


def test_dark_relaxation():
    """Gas outside the discharge that relaxes on its own makes a bigger pumping cell cost polarization."""
    ss = [m.compute(dict(FLUX, V_pump_cc=V, T1_dark_on=True, T1_dark_h=0.5))["ss_target"] for V in (20.0, 49.0, 100.0, 400.0)]
    assert all(b < a for a, b in zip(ss, ss[1:])), ss
    assert ss[0] < m.compute(dict(FLUX, V_pump_cc=20.0))["ss_target"]


def test_photon_efficiency():
    r = m.compute(dict(FLUX, Pabs_on=True, Pabs_seal_mW=200.0, Pabs_mW=500.0))
    pe = r["pump_eff"]
    photons = 0.2 / m.HNU
    assert abs(photons / (0.2 * 5.45e18) - 1) < 2e-3                  # 1 W is 5.45e18 photons/s at 1083 nm
    assert abs(pe["eta"] - pe["Q0"] / photons) < 1e-12 * pe["eta"]
    assert abs(pe["s"] - 2.5) < 1e-15                                 # the powers set s, not the slider
    assert not any("per absorbed photon" in w for w in r["warnings"])
    r = m.compute(dict(FLUX, Pabs_on=True, Pabs_seal_mW=1.0, Pabs_mW=1.0))
    assert r["pump_eff"]["eta"] > 2 and any("per absorbed photon" in w for w in r["warnings"])


def test_warnings():
    r = m.compute(dict(FLUX, V_pump_cc=500.0, s_dis=0.1, on_h=1.0))
    assert r["pump_eff"]["T_b_eff"] > 3600.0 and any("can't fill a pumping cell" in w for w in r["warnings"])
    assert not any("can't fill" in w for w in m.compute(dict(VOLUME, V_pump_cc=500.0, on_h=0.01))["warnings"])
    r = m.compute(dict(FLUX, tau_p_s=100.0, s_dis=10.0))              # T1p < T_b: no room for optical pumping
    assert r["ok"] and any("must be longer than the sealed-cell build-up time" in w for w in r["warnings"])


PARITY_CASES = [
    dict(V_pump_cc=250.0, s_dis=2.0, T1_dark_on=True, T1_dark_h=3.0),
    dict(V_pump_cc=15.0, s_dis=0.3, mode="single", on_h=2.0),
    dict(V_pump_cc=120.0, V_seal_cc=30.0, Pabs_on=True, Pabs_seal_mW=1.0, Pabs_mW=3.0, s_dis=1.5),
    dict(tau_p_s=100.0, s_dis=4.0),
    dict(pump_scaling="volume", V_pump_cc=300.0),
]


def test_js_parity():
    if not shutil.which("node"):
        print("parity   skipped: node is not on the PATH")
        return
    for over in PARITY_CASES:
        P = dict(FLUX, **over)
        py, js = json.loads(m.run(json.dumps(P))), run_js(P, "run")
        assert py["ok"] and js["ok"], (py.get("error"), js.get("error"))
        close(py, js, "compute", 1e-12)
        P = dict(P, **SCHED)
        py, js = json.loads(m.run_schedule(json.dumps(P))), run_js(P)
        assert py["ok"] and js["ok"], (py.get("error"), js.get("error"))
        close(py, js, "schedule", 1e-12)
    print("parity   %d flux-mode cases agree, main page and schedule" % len(PARITY_CASES))


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("ok      ", name)
