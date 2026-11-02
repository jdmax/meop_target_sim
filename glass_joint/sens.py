import sys, pickle, numpy as np
from fe_model import Case
from run_cases import run
steps = np.array([270, 240, 200, 160, 120, 80, 40, 15, 2.0])
tl = np.arange(len(steps))
V = {
 # mesh convergence
 "A_h10": Case("A_h10","collar", h0=0.010, al_yield=300.0, al_H=1e3),
 "A_h025": Case("A_h025","collar", h0=0.0025, al_yield=300.0, al_H=1e3),
 "B_h10": Case("B_h10","feather", h0=0.010, al_yield=40.0),
 "B_h025": Case("B_h025","feather", h0=0.0025, al_yield=40.0),
 # sensitivity (feather)
 "B_e": Case("B_e","feather", label="Al elastic (no yield)"),
 "B_y100": Case("B_y100","feather", al_yield=100.0, label="Al yield 100 MPa (work-hardened)"),
 "B_tip03": Case("B_tip03","feather", tip=0.03, al_yield=40.0, label="tip 0.03 mm"),
 "B_tip25": Case("B_tip25","feather", tip=0.25, al_yield=40.0, label="tip 0.25 mm"),
 "B_tip50": Case("B_tip50","feather", tip=0.5, al_yield=40.0, label="tip 0.5 mm"),
 "B_bond150": Case("B_bond150","feather", bond=0.15, al_yield=40.0, label="bond line 150 um"),
 "B_bond040": Case("B_bond040","feather", bond=0.04, al_yield=40.0, label="bond line 40 um"),
 "B_ep06": Case("B_ep06","feather", ep_scale=0.6, al_yield=40.0, label="epoxy modulus x0.6"),
 "B_ep15": Case("B_ep15","feather", ep_scale=1.5, al_yield=40.0, label="epoxy modulus x1.5"),
 "B_taper9": Case("B_taper9","feather", taper=9.0, al_yield=40.0, label="taper 9 mm"),
 "B_wall05": Case("B_wall05","feather", wall=0.5, al_yield=40.0, label="feather wall 0.5 mm"),
}
for k in sys.argv[1:]:
    res = run(V[k], snaps=steps, tlist=tl)
    st = res["stats"][-1]
    print(f"{k:10s} outer s1 d=25/50/100/250um " + "/".join(f"{v:.0f}" for v in st['outer_s1_at_d'][1:5])
          + f"  bend {st['outer_s1_bend']:.1f}  tau edge/deep {st['bond_tau_edge']:.0f}/{st['bond_tau_deep']:.0f}"
          + f"  peeq body {st['al_peeq_body']:.1e}  ({res['wall']:.0f}s, {res['nn']} nodes)", flush=True)
    pickle.dump(res, open(f"res_{k}.pkl","wb"))
