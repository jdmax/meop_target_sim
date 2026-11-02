"""Axisymmetric thermoelastic(-plastic) FE model: borosilicate tube bonded into an Al heat sink.

Geometry (z along tube axis, z = 0 at the warm end; mm):
  glass tube   r in [3, 4],      z in [0, 50]      (8 mm OD, 1 mm wall)
  warm zone    z in [0, 10]  held at T_hot = 292 K, mechanically free
  free span    z in [10, 40] steady 1-D conduction, k(T) via Kirchhoff transform
  bond         z in [40, 50] Stycast 2850FT bond line, thickness b, between glass OD and Al
  Al           r in [4+b, 4+b+t(z)], z in [40, 54]   (t(z) defines collar vs. feather)
Everything at z >= 40 sits at the sink temperature.  Stress-free state: 293 K (room-temp cure).

Elements: 9-node Lagrange quads, 3x3 Gauss.  Al optionally J2 elastoplastic with linear isotropic
hardening (radial return, consistent tangent).  Glass and epoxy linear elastic.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.integrate import cumulative_trapezoid
import materials as M

GLASS, EPOXY, AL = 0, 1, 2
R_IN, R_OUT = 3.0, 4.0
Z_HOT, Z_BOND, Z_GLASS_END, Z_AL_END = 10.0, 40.0, 50.0, 54.0


# ---------------------------------------------------------------- mesh utilities
def graded(a, b, ha, hb, hmax, g=0.15):
    """Element edges on [a,b], size ~ha at a, ~hb at b, growing ~(1+g) per element up to hmax."""
    x = a + (b - a) * np.linspace(0, 1, 200001)
    h = np.minimum.reduce([np.full_like(x, hmax), ha + g * (x - a), hb + g * (b - x)])
    s = cumulative_trapezoid(1.0 / h, x, initial=0.0)
    n = max(1, int(np.ceil(s[-1])))
    return np.interp(np.linspace(0, s[-1], n + 1), s, x)


def with_mid(edges):
    out = np.empty(2 * len(edges) - 1)
    out[0::2] = edges
    out[1::2] = 0.5 * (edges[:-1] + edges[1:])
    return out


class Case:
    def __init__(self, name, kind, h0=0.005, bond=0.075, tip=0.1, taper=6.0, wall=1.0,
                 collar_r=12.0, al_yield=None, al_H=1.5e3, n_al=24, label=None, ep_scale=1.0):
        self.name, self.kind, self.h0, self.bond = name, kind, h0, bond
        self.tip, self.taper, self.wall, self.collar_r = tip, taper, wall, collar_r
        self.al_yield, self.al_H, self.n_al = al_yield, al_H, n_al
        self.label = label or name
        self.ep_scale = ep_scale

    def t_al(self, z):
        """Radial thickness of the Al band at axial position z."""
        r0 = R_OUT + self.bond
        full = self.collar_r - r0
        z = np.asarray(z, float)
        if self.kind == "collar":
            return np.full_like(z, full)
        # feather: tip -> wall over the taper, wall to the glass end, then into a thick block
        t = np.interp(z, [Z_BOND, Z_BOND + self.taper, Z_GLASS_END, Z_GLASS_END + 1.0, Z_AL_END],
                      [self.tip, self.wall, self.wall, full, full])
        return t


def build_mesh(case):
    h0 = case.h0
    zb = [0.0, Z_HOT, Z_BOND, Z_BOND + case.taper, Z_GLASS_END, Z_GLASS_END + 1.0, Z_AL_END]
    hz = [1.0, 0.2, h0, 0.1, 4 * h0, 0.2, 0.5]
    zedges = [np.array([0.0])]
    for i in range(len(zb) - 1):
        e = graded(zb[i], zb[i + 1], hz[i], hz[i + 1], 1.0)
        zedges.append(e[1:])
    zedges = np.concatenate(zedges)
    zn = with_mid(zedges)                      # all z node coordinates (shared by every band)

    r_glass = with_mid(graded(R_IN, R_OUT, 0.02, h0, 0.1))
    n_ep = max(4, int(np.ceil(case.bond / max(h0, 0.008))))
    r_ep = with_mid(np.linspace(R_OUT, R_OUT + case.bond, n_ep + 1))
    # Al band: parametric s in [0,1], graded toward the bonded inner face
    s_edges = graded(0.0, 1.0, 0.01, 0.12, 0.12, g=0.25)
    s_al = with_mid(s_edges)

    i_bond = np.searchsorted(zn, Z_BOND - 1e-9)
    i_gend = np.searchsorted(zn, Z_GLASS_END - 1e-9)
    assert abs(zn[i_bond] - Z_BOND) < 1e-9 and abs(zn[i_gend] - Z_GLASS_END) < 1e-9

    coords, elems, mats = [], [], []
    node_count = 0

    def add_band(zsl, rgrid_fn, mat):
        nonlocal node_count
        zz = zn[zsl]
        R = rgrid_fn(zz)                         # (nz, nr) radial node coords
        nz, nr = R.shape
        Z = np.repeat(zz[:, None], nr, axis=1)
        ids = node_count + np.arange(nz * nr).reshape(nz, nr)
        coords.append(np.column_stack([R.ravel(), Z.ravel()]))
        node_count += nz * nr
        ez, er = (nz - 1) // 2, (nr - 1) // 2
        for p in range(ez):
            for q in range(er):
                blk = ids[2 * p:2 * p + 3, 2 * q:2 * q + 3]   # [a(z), b(r)]
                elems.append(blk.ravel())                  # local k = 3*a + b
                mats.append(mat)

    add_band(slice(0, i_gend + 1), lambda zz: np.tile(r_glass, (len(zz), 1)), GLASS)
    add_band(slice(i_bond, i_gend + 1), lambda zz: np.tile(r_ep, (len(zz), 1)), EPOXY)
    r0 = R_OUT + case.bond
    add_band(slice(i_bond, len(zn)), lambda zz: r0 + np.outer(case.t_al(zz), s_al), AL)

    X = np.vstack(coords)
    elems = np.array(elems)
    mats = np.array(mats)
    # merge coincident nodes
    key = np.round(X / 1e-7).astype(np.int64)
    _, uniq_idx, inv = np.unique(key[:, 0] * 10**12 + key[:, 1], return_index=True, return_inverse=True)
    Xu = X[uniq_idx]
    elems = inv[elems]
    if case.kind == "glass_only":           # drop epoxy + Al, compact node numbering
        keep = mats == GLASS
        elems, mats = elems[keep], mats[keep]
        used, newid = np.unique(elems, return_inverse=True)
        Xu, elems = Xu[used], newid.reshape(elems.shape)
    return Xu, elems, mats


# ---------------------------------------------------------------- element machinery
_gp1 = np.array([-np.sqrt(0.6), 0.0, np.sqrt(0.6)])
_gw1 = np.array([5 / 9, 8 / 9, 5 / 9])


def _L(x):
    return np.array([0.5 * x * (x - 1), 1 - x * x, 0.5 * x * (x + 1)])


def _dL(x):
    return np.array([x - 0.5, -2 * x, x + 0.5])


def shape_tables():
    N = np.zeros((9, 9)); dNxi = np.zeros((9, 9)); dNeta = np.zeros((9, 9)); W = np.zeros(9)
    g = 0
    for i, xi in enumerate(_gp1):          # xi -> local a (z direction)
        for j, eta in enumerate(_gp1):     # eta -> local b (r direction)
            La, Lb, dLa, dLb = _L(xi), _L(eta), _dL(xi), _dL(eta)
            for a in range(3):
                for b in range(3):
                    k = 3 * a + b
                    N[g, k] = La[a] * Lb[b]
                    dNxi[g, k] = dLa[a] * Lb[b]
                    dNeta[g, k] = La[a] * dLb[b]
            W[g] = _gw1[i] * _gw1[j]
            g += 1
    return N, dNxi, dNeta, W


class Model:
    def __init__(self, case):
        self.case = case
        X, E, mats = build_mesh(case)
        self.X, self.E, self.mats = X, E, mats
        self.nn, self.ne = len(X), len(E)
        N, dNxi, dNeta, W = shape_tables()
        self.N = N
        xe = X[E]                                            # (ne, 9, 2)  [r, z]
        J = np.empty((self.ne, 9, 2, 2))
        J[:, :, 0, 0] = np.einsum('gk,ek->eg', dNxi, xe[:, :, 0])
        J[:, :, 0, 1] = np.einsum('gk,ek->eg', dNxi, xe[:, :, 1])
        J[:, :, 1, 0] = np.einsum('gk,ek->eg', dNeta, xe[:, :, 0])
        J[:, :, 1, 1] = np.einsum('gk,ek->eg', dNeta, xe[:, :, 1])
        det = J[..., 0, 0] * J[..., 1, 1] - J[..., 0, 1] * J[..., 1, 0]
        assert np.all(np.abs(det) > 0)
        inv = np.empty_like(J)
        inv[..., 0, 0] = J[..., 1, 1] / det; inv[..., 1, 1] = J[..., 0, 0] / det
        inv[..., 0, 1] = -J[..., 0, 1] / det; inv[..., 1, 0] = -J[..., 1, 0] / det
        dNr = inv[..., 0, 0, None] * dNxi[None] + inv[..., 0, 1, None] * dNeta[None]
        dNz = inv[..., 1, 0, None] * dNxi[None] + inv[..., 1, 1, None] * dNeta[None]
        rg = np.einsum('gk,ek->eg', N, xe[:, :, 0])
        zg = np.einsum('gk,ek->eg', N, xe[:, :, 1])
        self.rg, self.zg = rg, zg
        B = np.zeros((self.ne, 9, 4, 18))
        B[:, :, 0, 0::2] = dNr
        B[:, :, 1, 1::2] = dNz
        B[:, :, 2, 0::2] = N[None] / rg[..., None]
        B[:, :, 3, 0::2] = dNz
        B[:, :, 3, 1::2] = dNr
        self.B = B
        self.w = rg * np.abs(det) * W[None]                  # 2*pi omitted
        dof = np.empty((self.ne, 18), dtype=np.int64)
        dof[:, 0::2] = 2 * E; dof[:, 1::2] = 2 * E + 1
        self.dof = dof
        self.Irow = np.repeat(dof, 18, axis=1).ravel()
        self.Icol = np.tile(dof, (1, 18)).ravel()
        self.ndof = 2 * self.nn
        # anchor: axial dof at the outer, far corner of the Al block
        far = np.where(np.abs(X[:, 1] - X[:, 1].max()) < 1e-9)[0]
        anchor = far[np.argmax(X[far, 0])]
        self.fixed = np.array([2 * anchor + 1])
        self.free = np.setdiff1d(np.arange(self.ndof), self.fixed)
        # state
        self.u = np.zeros(self.ndof)
        self.ep = np.zeros((self.ne, 9, 4))
        self.alpha = np.zeros((self.ne, 9))
        self.matg = np.repeat(mats[:, None], 9, axis=1)
        self.Minv = np.linalg.inv(N)                         # gp values -> nodal (extrapolation)

    # ------------------------------------------------------------ fields
    def temperature(self, z, Ts, Th=292.0):
        z = np.asarray(z, float)
        T = np.where(z <= Z_HOT, Th, Ts).astype(float)
        m = (z > Z_HOT) & (z < Z_BOND)
        Ks, Kh = M.Kint(Ts), M.Kint(Th)
        T[m] = M.Kinv(Ks + (Kh - Ks) * (Z_BOND - z[m]) / (Z_BOND - Z_HOT))
        return T

    def props(self, Ts):
        T = self.temperature(self.zg, Ts)
        m = self.matg
        E = np.where(m == GLASS, M.E_GLASS, np.where(m == EPOXY, self.case.ep_scale * M.E_epoxy(T), M.E_al(T)))
        nu = np.where(m == GLASS, M.NU_GLASS, np.where(m == EPOXY, M.NU_EPOXY, M.NU_AL))
        eth = np.where(m == GLASS, M.eps_pyrex(T), np.where(m == EPOXY, M.eps_epoxy(T), M.eps_al(T)))
        return T, E, nu, eth

    # ------------------------------------------------------------ constitutive
    def constitutive(self, eps, E, nu, eth, ep, alpha):
        lam = E * nu / ((1 + nu) * (1 - 2 * nu)); mu = E / (2 * (1 + nu)); kap = lam + 2 * mu / 3
        mvec = np.array([1.0, 1.0, 1.0, 0.0])
        ee = eps - eth[..., None] * mvec - ep
        tr = ee[..., 0] + ee[..., 1] + ee[..., 2]
        sig = lam[..., None] * tr[..., None] * mvec + 2 * mu[..., None] * ee * np.array([1, 1, 1, 0.5])
        C = (lam[..., None, None] * np.outer(mvec, mvec)
             + mu[..., None, None] * np.diag([2.0, 2.0, 2.0, 1.0]))
        C = np.broadcast_to(C, sig.shape[:-1] + (4, 4)).copy()
        ep_new, al_new = ep.copy(), alpha.copy()
        if self.case.al_yield is not None and np.any(self.matg == AL):
            H = self.case.al_H
            p = (sig[..., 0] + sig[..., 1] + sig[..., 2]) / 3
            s = sig - p[..., None] * mvec
            sn = np.sqrt(s[..., 0]**2 + s[..., 1]**2 + s[..., 2]**2 + 2 * s[..., 3]**2)
            q = np.sqrt(1.5) * sn
            f = q - (self.case.al_yield + H * alpha)
            y = (self.matg == AL) & (f > 0)
            if np.any(y):
                dg = f[y] / (3 * mu[y] + H)
                nh = s[y] / sn[y][:, None]
                sig[y] -= (2 * mu[y] * dg * np.sqrt(1.5))[:, None] * nh
                dep = (dg * np.sqrt(1.5))[:, None] * nh
                dep[:, 3] *= 2.0
                ep_new[y] += dep
                al_new[y] += dg
                beta = 1 - 3 * mu[y] * dg / q[y]
                thb = 3 * mu[y] / (3 * mu[y] + H) - (1 - beta)
                Idev = np.diag([1, 1, 1, 0.5]) - np.outer(mvec, mvec) / 3
                C[y] = (kap[y][:, None, None] * np.outer(mvec, mvec)
                        + (2 * mu[y] * beta)[:, None, None] * Idev
                        - (2 * mu[y] * thb)[:, None, None] * np.einsum('ni,nj->nij', nh, nh))
        return sig, C, ep_new, al_new

    # ------------------------------------------------------------ solve one load step
    def step(self, Ts, tol=1e-9, maxit=40, verbose=False):
        T, E, nu, eth = self.props(Ts)
        u = self.u.copy()
        ref = None
        for it in range(maxit):
            eps = np.einsum('egij,ej->egi', self.B, u[self.dof])
            sig, C, ep_n, al_n = self.constitutive(eps, E, nu, eth, self.ep, self.alpha)
            fe = np.einsum('egij,egi,eg->ej', self.B, sig, self.w)
            R = np.zeros(self.ndof); np.add.at(R, self.dof, fe)
            Rf = R[self.free]
            if ref is None:
                # reference force: thermal load from a zero-displacement state
                ref = max(np.linalg.norm(Rf), 1e-12) if it == 0 else ref
            nr = np.linalg.norm(Rf)
            if verbose:
                print(f"    it {it}: |R| = {nr:.3e}  (rel {nr/ref:.2e})")
            if nr <= tol * ref and it > 0:
                break
            CB = np.einsum('egkl,eglj->egkj', C, self.B)
            Ke = np.einsum('egki,egkj,eg->eij', self.B, CB, self.w)
            K = sp.csc_matrix((Ke.ravel(), (self.Irow, self.Icol)), shape=(self.ndof, self.ndof))
            Kff = K[self.free][:, self.free]
            du = spla.spsolve(Kff, -Rf)
            u[self.free] += du
        else:
            print(f"  WARNING: {self.case.name} Ts={Ts:.1f} not converged, rel {nr/ref:.2e}")
        self.u = u
        self.ep, self.alpha = ep_n, al_n
        self.sig, self.T = sig, T
        return it

    # ------------------------------------------------------------ post-processing
    def nodal(self, values_gp, mat):
        """Extrapolate gp values (ne,9,k) to nodes, averaged over elements of one material."""
        sel = self.mats == mat
        vn = np.einsum('kg,egc->ekc', self.Minv, values_gp[sel])
        nodes = self.E[sel]
        acc = np.zeros((self.nn, values_gp.shape[-1])); cnt = np.zeros(self.nn)
        np.add.at(acc, nodes.ravel(), vn.reshape(-1, values_gp.shape[-1]))
        np.add.at(cnt, nodes.ravel(), 1)
        idx = np.where(cnt > 0)[0]
        return idx, acc[idx] / cnt[idx, None]


def principal(s):
    """s: (..., 4) [rr, zz, tt, rz] -> (s1 max, s3 min)."""
    c = 0.5 * (s[..., 0] + s[..., 1])
    rad = np.sqrt((0.5 * (s[..., 0] - s[..., 1]))**2 + s[..., 3]**2)
    a, b, t = c + rad, c - rad, s[..., 2]
    return np.maximum(np.maximum(a, b), t), np.minimum(np.minimum(a, b), t)


def von_mises(s):
    return np.sqrt(0.5 * ((s[..., 0] - s[..., 1])**2 + (s[..., 1] - s[..., 2])**2 + (s[..., 2] - s[..., 0])**2)
                   + 3 * s[..., 3]**2)
