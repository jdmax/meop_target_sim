"""Temperature-dependent material properties, 1-300 K.

Units: mm, MPa, K.  Thermal strain is (L(T) - L(293 K)) / L(293 K), i.e. negative on cooling.

Sources / anchors
-----------------
Al 6061-T6 and Al 1100 : NIST cryogenic materials fits (6061-T6 page) for linear expansion and
                          Young's modulus.  1100 contracts within ~1 % of 6061; the same curve is used.
Pyrex 7740             : total contraction 293->4 K = 5.6e-4 (Ekin, Table A6.4); CTE at 293 K ~3.0e-6/K.
                          Curve shape: single-Einstein alpha(T) fitted to those two anchors.
Stycast 2850FT/24LV    : contraction 293->77 K = 4.78e-3 (FNAL TS-SSC-90-060), 293->4 K = 5.08e-3 (Ekin).
                          Curve shape: single-Einstein alpha(T) fitted to those two anchors.
Pyrex k(T)             : smooth table, integral 4-300 K ~ 200 W/m (consistent with literature).
"""
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import brentq

TREF = 293.0
_Tg = np.linspace(0.2, 420.0, 20000)


def _einstein(x):
    x = np.clip(x, 1e-6, 300.0)
    ex = np.exp(x)
    return x * x * ex / (ex - 1.0) ** 2


def _shape(theta):
    c = _einstein(theta / _Tg)
    cum = cumulative_trapezoid(c, _Tg, initial=0.0)
    return c, cum - np.interp(TREF, _Tg, cum)  # integral from TREF to T (<0 below TREF)


def _fit_pyrex(alpha293=3.0e-6, total=5.6e-4):
    def f(th):
        c, s = _shape(th)
        return np.interp(TREF, _Tg, c) / (-np.interp(4.0, _Tg, s)) - alpha293 / total
    th = brentq(f, 20.0, 2000.0)
    c, s = _shape(th)
    A = total / (-np.interp(4.0, _Tg, s))
    return th, A * s, A * c


def _fit_epoxy(e77=4.78e-3, e4=5.08e-3):
    def f(th):
        _, s = _shape(th)
        return np.interp(77.0, _Tg, s) / np.interp(4.0, _Tg, s) - e77 / e4
    th = brentq(f, 5.0, 3000.0)
    c, s = _shape(th)
    A = e4 / (-np.interp(4.0, _Tg, s))
    return th, A * s, A * c


THETA_PYREX, _EPS_PYREX, _ALPHA_PYREX = _fit_pyrex()
THETA_EPOXY, _EPS_EPOXY, _ALPHA_EPOXY = _fit_epoxy()


def eps_al(T):
    T = np.asarray(T, float)
    Tc = np.maximum(T, 18.0)
    y = -4.1277e2 - 3.0389e-1 * Tc + 8.7696e-3 * Tc**2 - 9.9821e-6 * Tc**3
    y = np.where(T < 18.0, -415.45, y)
    return y * 1e-5


def eps_pyrex(T):
    return np.interp(T, _Tg, _EPS_PYREX)


def eps_epoxy(T):
    return np.interp(T, _Tg, _EPS_EPOXY)


def alpha_pyrex(T):
    return np.interp(T, _Tg, _ALPHA_PYREX)


def E_al(T):  # MPa, NIST 6061-T6 fit (2-295 K)
    T = np.clip(np.asarray(T, float), 2.0, 295.0)
    y = 7.771221e1 + 1.030646e-2 * T - 2.924100e-4 * T**2 + 8.993600e-7 * T**3 - 1.070900e-9 * T**4
    return y * 1e3


def E_epoxy(T, E293=7.0e3, E0=13.0e3):  # MPa; ~linear stiffening on cooling (assumption)
    T = np.clip(np.asarray(T, float), 0.0, 293.0)
    return E0 + (E293 - E0) * T / 293.0


E_GLASS, NU_GLASS = 63.0e3, 0.20
NU_EPOXY, NU_AL = 0.33, 0.33

# Pyrex thermal conductivity, W/m/K
_kT = np.array([1.0, 2.0, 4.0, 10.0, 20.0, 40.0, 77.0, 100.0, 150.0, 200.0, 250.0, 300.0, 320.0])
_kv = np.array([0.025, 0.06, 0.10, 0.12, 0.16, 0.28, 0.48, 0.58, 0.75, 0.88, 1.00, 1.10, 1.14])
_KT = np.linspace(1.0, 320.0, 64000)
_Kk = np.exp(np.interp(np.log(_KT), np.log(_kT), np.log(_kv)))
_KINT = cumulative_trapezoid(_Kk, _KT, initial=0.0)  # integral from 1 K, W/m


def k_pyrex(T):
    return np.exp(np.interp(np.log(T), np.log(_kT), np.log(_kv)))


def Kint(T):
    return np.interp(T, _KT, _KINT)


def Kinv(K):
    return np.interp(K, _KINT, _KT)


if __name__ == "__main__":
    print(f"Pyrex Einstein theta = {THETA_PYREX:.1f} K; alpha(293)={alpha_pyrex(293)*1e6:.2f} ppm/K")
    print(f"Epoxy Einstein theta = {THETA_EPOXY:.1f} K; alpha(293)={np.interp(293,_Tg,_ALPHA_EPOXY)*1e6:.1f} ppm/K")
    for T in [4, 20, 40, 77, 100, 150, 200, 250, 293]:
        print(f"T={T:4d}  Pyrex {eps_pyrex(T)*1e5:8.2f}  Epoxy {eps_epoxy(T)*1e5:8.1f}  Al {eps_al(T)*1e5:8.1f} (x1e-5)  "
              f"E_Al={E_al(T)/1e3:.1f} GPa  k={k_pyrex(T):.2f}")
    print("int k dT 4-300 K =", Kint(300) - Kint(4), "W/m")
    # check vs design-reference profile for a 100 mm tube, 5 K -> 300 K
    for d in [5, 17, 55]:
        K = Kint(5) + (Kint(300) - Kint(5)) * d / 100
        print(f"  100 mm tube: {d} mm from cold end -> {Kinv(K):.0f} K")
