"""Stress intensity for circumferential surface cracks in the glass, from the FE stress field.

Crack: ring (edge) crack of depth a, growing radially inward from the outer surface at plane z0,
opened by sigma_zz.  K from Glinka's universal weight function for an edge crack in a half-plane
(M1=0.0719768, M2=0.246984, M3=0.529659; reproduces K = 1.1215 sigma sqrt(pi a) for uniform stress).
A long ring crack is conservative relative to a semi-elliptical (thumbnail) flaw of the same depth.
"""
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from fe_model import R_OUT, Z_BOND

KIC = 0.77          # MPa sqrt(m), Pyrex 7740
M1, M2, M3 = 0.0719768, 0.246984, 0.529659
_v, _wv = np.polynomial.legendre.leggauss(48)
_v = 0.5 * (_v + 1); _wv = 0.5 * _wv


def glass_grid(X, sig):
    z = np.round(X[:, 1], 9); r = np.round(X[:, 0], 9)
    zu, ru = np.unique(z), np.unique(r)
    G = np.full((len(zu), len(ru), sig.shape[1]), np.nan)
    G[np.searchsorted(zu, z), np.searchsorted(ru, r)] = sig
    assert not np.isnan(G).any()
    return zu, ru, G


def K_profile(X, sig, a_mm, z_lo=39.0, z_hi=Z_BOND):
    """Max K (MPa sqrt m) over crack planes z0 in [z_lo, z_hi) for each depth a (mm)."""
    zu, ru, G = glass_grid(X, sig)
    f = RegularGridInterpolator((zu, ru), G[:, :, 1], bounds_error=False, fill_value=None)
    planes = zu[(zu >= z_lo) & (zu < z_hi)]
    Kmax = np.zeros(len(a_mm)); zbest = np.zeros(len(a_mm))
    for i, a in enumerate(a_mm):
        x = a * (1 - _v**2)                          # depth from the outer surface
        rr = R_OUT - x
        pts = np.column_stack([np.repeat(planes, len(x)), np.tile(rr, len(planes))])
        s = f(pts).reshape(len(planes), len(x))
        wfun = 1 + M1 * _v + M2 * _v**2 + M3 * _v**3
        K = (2 / np.sqrt(2 * np.pi)) * 2 * np.sqrt(a * 1e-3) * (s * wfun) @ _wv
        j = np.argmax(K); Kmax[i] = K[j]; zbest[i] = planes[j]
    return Kmax, zbest


def critical_depth(a_mm, Kmax, kic=KIC):
    above = np.where(Kmax >= kic)[0]
    if len(above) == 0:
        return np.inf
    i = above[0]
    if i == 0:
        return a_mm[0]
    # log interpolation
    return float(np.exp(np.interp(np.log(kic), np.log(Kmax[i-1:i+1]), np.log(a_mm[i-1:i+1]))))


if __name__ == "__main__":
    # weight-function self-check: uniform stress
    a = 0.05
    K = (2/np.sqrt(2*np.pi))*2*np.sqrt(a*1e-3)*((1+M1*_v+M2*_v**2+M3*_v**3)@_wv)*100.0
    print("uniform 100 MPa, a=50um: K =", K, " expected", 1.1215*100*np.sqrt(np.pi*a*1e-3))
