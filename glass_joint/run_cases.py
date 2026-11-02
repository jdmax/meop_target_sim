import numpy as np, time, pickle, sys
from fe_model import *
import materials as M

T_HOT = 292.0
EDGE_D = np.array([0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0])
times = np.arange(0, 10.01, 0.5)                     # h
Tsink = np.maximum(T_HOT - 29.0 * times, 2.0)

def corner_dist(X, zc=Z_BOND, rc=R_OUT):
    return np.hypot(X[:, 0] - rc, X[:, 1] - zc)

def run(case, snaps=Tsink, tlist=times, keep_field=True, verbose=False, tol=1e-8):
    m = Model(case)
    out = dict(case=case.name, label=case.label, nn=m.nn, ne=m.ne, t=[], Ts=[], prof=[], stats=[])
    t0 = time.time()
    for t, Ts in zip(tlist, snaps):
        it = m.step(Ts, tol=tol)
        idx, sg = m.nodal(m.sig, GLASS)
        Xg = m.X[idx]
        s1, s3 = principal(sg)
        d = corner_dist(Xg)
        d50 = np.hypot(Xg[:, 0] - R_OUT, Xg[:, 1] - Z_GLASS_END)
        st = {f"s1_ex{rho}": s1[(d >= rho) & (d50 >= rho)].max() for rho in [0.0, 0.025, 0.05, 0.1, 0.2, 0.5]}
        st["s3_min"] = s3.min()
        # free-span stresses, away from both thermal kinks and the joint (12 < z < 37)
        fs = (Xg[:, 1] > 15) & (Xg[:, 1] < 30)
        st["free_s1"] = s1[fs].max(); st["free_s3"] = s3[fs].min()
        st["free_szz_abs"] = np.abs(sg[fs, 1]).max()
        st["edge_d"] = EDGE_D
        if np.any(m.mats == EPOXY):
            ie, se = m.nodal(m.sig, EPOXY)
            e1, _ = principal(se)
            st["ep_s1"] = e1.max()
            # bond mid-line (non-singular) shear and peel
            xg = np.einsum('gk,ekc->egc', m.N, m.X[m.E])
            sel = m.mats == EPOXY
            x = xg[sel].reshape(-1, 2); sgp = m.sig[sel].reshape(-1, 4)
            mid = np.abs(x[:, 0] - (R_OUT + 0.5 * case.bond)) < 0.25 * case.bond
            st["bond_tau_edge"] = np.abs(sgp[mid & (x[:, 1] < 42), 3]).max()
            st["bond_tau_deep"] = np.abs(sgp[mid & (x[:, 1] > 48), 3]).max()
            st["bond_peel_edge"] = sgp[mid & (x[:, 1] < 42), 0].max()
            al = m.matg == AL
            st["al_vm"] = von_mises(m.sig[al]).max()
            st["al_peeq"] = m.alpha[al].max()
            # plastic strain in the feather body, away from singular corners (mid-thickness, 41 < z < 48)
            xa = xg[m.mats == AL].reshape(-1, 2); aa = m.alpha[m.mats == AL].ravel()
            r0 = R_OUT + case.bond
            body = (xa[:, 1] > 41) & (xa[:, 1] < 48) & (xa[:, 0] > r0 + 0.3 * case.t_al(xa[:, 1])) & (xa[:, 0] < r0 + 0.7 * case.t_al(xa[:, 1]))
            st["al_peeq_body"] = aa[body].max() if body.any() else 0.0
        else:
            for k in ["ep_s1", "bond_tau_edge", "bond_tau_deep", "bond_peel_edge", "al_vm", "al_peeq", "al_peeq_body"]:
                st[k] = np.nan
        # where is the glass s1 max (excluding 0.1 mm around singular corners)?
        msk = (d >= 0.1) & (d50 >= 0.1)
        j = np.argmax(np.where(msk, s1, -np.inf))
        st["s1_loc"] = tuple(Xg[j])
        # surface profiles along z
        outer = np.abs(Xg[:, 0] - R_OUT) < 1e-9
        inner = np.abs(Xg[:, 0] - R_IN) < 1e-9
        prof = {}
        for nm, msk2 in [("outer", outer), ("inner", inner)]:
            o = np.argsort(Xg[msk2, 1])
            prof[nm] = dict(z=Xg[msk2, 1][o], s=sg[msk2][o], s1=s1[msk2][o])
        zo, s1o, so = prof["outer"]["z"], prof["outer"]["s1"], prof["outer"]["s"]
        zi, s1i, si = prof["inner"]["z"], prof["inner"]["s1"], prof["inner"]["s"]
        st["outer_s1_at_d"] = np.interp(Z_BOND - EDGE_D, zo, s1o)
        st["outer_szz_at_d"] = np.interp(Z_BOND - EDGE_D, zo, so[:, 1])
        w = (zo > 36) & (zo < Z_BOND - 0.2)
        st["outer_s1_bend"] = s1o[w].max(); st["inner_s1_bend"] = s1i[(zi > 36) & (zi < 44)].max()
        out.setdefault("glass_sig_t", []).append(sg.astype(np.float32))
        out["t"].append(t); out["Ts"].append(Ts); out["prof"].append(prof); out["stats"].append(st)
        if verbose:
            print(f"  t={t:4.1f} h Ts={Ts:6.1f} | glass outer s1 at d=25/50/100/250um: "
                  + "/".join(f"{v:.0f}" for v in st['outer_s1_at_d'][1:5])
                  + f"  bend-zone outer {st['outer_s1_bend']:.1f} inner {st['inner_s1_bend']:.1f}"
                  f" | free-span s1 {st['free_s1']:.2f} s3 {st['free_s3']:.2f}"
                  f" | bond tau edge {st['bond_tau_edge']:.0f} deep {st['bond_tau_deep']:.0f}"
                  f" | Al peeq body {st['al_peeq_body']:.1e} | its={it}", flush=True)
    if keep_field:
        idx, sg = m.nodal(m.sig, GLASS)
        out["idx_glass"] = idx
        out["field"] = dict(X=m.X, E=m.E, mats=m.mats, idx_glass=idx, sig_glass=sg,
                            sig_gp=m.sig, alpha=m.alpha, T=m.T, u=m.u)
    out["wall"] = time.time() - t0
    return out

if __name__ == "__main__":
    which = sys.argv[1]
    cases = {
      "A": Case("A", "collar", label="A: solid 6061-T6 collar, epoxy-bonded (as originally posed)", al_yield=300.0, al_H=1.0e3),
      "G": Case("G", "glass_only", label="Glass tube alone (gradient only, no joint)"),
      "B": Case("B", "feather", label="B: 1100-O Al feather, 0.1 mm tip, epoxy-bonded (planned joint)", al_yield=40.0, al_H=1.5e3),
      "Bh": Case("Bh", "feather", label="B, work-hardened Al (yield 100 MPa)", al_yield=100.0, al_H=1.5e3),
      "Be": Case("Be", "feather", label="B, Al elastic (no yield)", al_yield=None),
      "Bt": Case("Bt", "feather", label="B, blunt 0.5 mm feather tip", tip=0.5, al_yield=40.0, al_H=1.5e3),
    }
    c = cases[which]
    print("Case", c.label, flush=True)
    res = run(c, verbose=True)
    print(f"  nodes {res['nn']} elems {res['ne']} wall {res['wall']:.0f} s")
    pickle.dump(res, open(f"res_{which}.pkl", "wb"))
