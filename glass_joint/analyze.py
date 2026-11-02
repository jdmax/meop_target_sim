import pickle, json, glob, numpy as np, os
from fracture import K_profile, critical_depth, KIC
from fe_model import R_IN, R_OUT, Z_BOND, principal
a_mm = np.logspace(np.log10(0.002), np.log10(0.6), 70)

def load(k): return pickle.load(open(f"res_{k}.pkl","rb"))

def fr(X, sig):
    res = {}
    for tag, c in [("01", 0.1), ("025", 0.25)]:
        K, zb = K_profile(X, sig, a_mm, z_lo=36.0, z_hi=Z_BOND - c + 1e-9)
        ac = critical_depth(a_mm, K)
        res["ac_um_" + tag] = float(ac * 1e3) if np.isfinite(ac) else None
        res["K25_" + tag] = float(np.interp(0.025, a_mm, K))
        res["K50_" + tag] = float(np.interp(0.050, a_mm, K))
        res["Kpeak_" + tag] = float(K.max()); res["a_at_Kpeak_" + tag] = float(a_mm[np.argmax(K)] * 1e3)
        if tag == "01": Kout = K
    res["ac_um"] = res["ac_um_01"]; res["K25"] = res["K25_01"]
    return Kout, None, res

out = {"KIC": KIC, "a_um": (a_mm*1e3).tolist()}
# time sweeps
for k in ["A","B","G"]:
    r = load(k); Xg = r["field"]["X"][r["idx_glass"]]
    ts = dict(t=[float(x) for x in r["t"]], Ts=[float(x) for x in r["Ts"]], ac_um=[], ac_um_025=[], K25=[], K50=[], Kpeak=[], Kpeak_025=[],
              outer_bend=[], s_250=[], s_50=[], tau_edge=[], tau_deep=[], peeq=[], free_s1=[], free_s3=[], glass_s3=[])
    prof = dict(z=None, outer_szz=[], inner_szz=[], outer_s1=[])
    for i, (st, sg) in enumerate(zip(r["stats"], r["glass_sig_t"])):
        K, ac, d = fr(Xg, sg.astype(float))
        ts["ac_um"].append(d["ac_um_01"]); ts["ac_um_025"].append(d["ac_um_025"]); ts["K25"].append(d["K25_01"]); ts["K50"].append(d["K50_01"]); ts["Kpeak"].append(d["Kpeak_01"]); ts["Kpeak_025"].append(d["Kpeak_025"])
        ts["outer_bend"].append(float(st["outer_s1_bend"]))
        ts["s_250"].append(float(st["outer_s1_at_d"][4])); ts["s_50"].append(float(st["outer_s1_at_d"][2]))
        for kk, src in [("tau_edge","bond_tau_edge"),("tau_deep","bond_tau_deep"),("peeq","al_peeq_body")]:
            v = st[src]; ts[kk].append(None if (v is None or np.isnan(v)) else float(v))
        ts["free_s1"].append(float(st["free_s1"])); ts["free_s3"].append(float(st["free_s3"]))
        ts["glass_s3"].append(float(st["s3_min"]))
        p = r["prof"][i]
        zo = p["outer"]["z"]; zi = p["inner"]["z"]
        # downsample profile: keep every node beyond 0.5 mm of the edge sparsely, dense near the joint
        keep = np.ones(len(zo), bool)
        if prof["z"] is None:
            sel = np.unique(np.concatenate([np.where((zo < 35) | (zo > 45))[0][::3], np.where((zo >= 35) & (zo <= 45))[0]]))
            prof["z"] = zo[sel].round(5).tolist(); prof["sel"] = sel
            assert np.allclose(zo, zi)
        sel = prof["sel"]
        prof["outer_szz"].append(p["outer"]["s"][sel, 1].round(4).tolist())
        prof["inner_szz"].append(p["inner"]["s"][sel, 1].round(4).tolist())
        prof["outer_s1"].append(p["outer"]["s1"][sel].round(4).tolist())
        if i == len(r["stats"]) - 1:
            ts["K_final"] = K.tolist()
    prof.pop("sel")
    # temperature profile along the tube for each snapshot
    from fe_model import Model, Case
    ts["label"] = r["label"]
    out[k] = dict(ts=ts, prof=prof)
    print(k, "final: ac=%s um, K25=%.2f, outer bend=%.1f, s@250um=%.1f, s@50um=%.1f, free s1=%.3f s3=%.3f, glass s3 min=%.0f" %
          (ts["ac_um"][-1], ts["K25"][-1], ts["outer_bend"][-1], ts["s_250"][-1], ts["s_50"][-1], ts["free_s1"][-1], ts["free_s3"][-1], ts["glass_s3"][-1]))
    print("   ac(0.1) vs time:", [None if v is None else round(v,1) for v in ts["ac_um"][::2]])
    print("   ac(0.25) vs time:", [None if v is None else round(v,1) for v in ts["ac_um_025"][::2]])

# temperature profiles
from fe_model import Model, Case
m = Model.__new__(Model)
zT = np.linspace(0, 50, 251)
out["Tprof"] = dict(z=zT.tolist(), T=[Model.temperature(m, zT, Ts).round(2).tolist() for Ts in out["A"]["ts"]["Ts"]])

# near-joint maps (final state), glass max principal
maps = {}
for k in ["A","B"]:
    r = load(k); f = r["field"]; X = f["X"][f["idx_glass"]]; s = f["sig_glass"]
    s1, s3 = principal(s)
    z = np.round(X[:,1], 9); rr = np.round(X[:,0], 9)
    zu, ru = np.unique(z), np.unique(rr)
    G = np.full((len(zu), len(ru)), np.nan); G[np.searchsorted(zu,z), np.searchsorted(ru,rr)] = s1
    Gz = np.full((len(zu), len(ru)), np.nan); Gz[np.searchsorted(zu,z), np.searchsorted(ru,rr)] = s[:,1]
    w = (zu >= 37.5) & (zu <= 43.0)
    zs = zu[w]; 
    # thin the z sampling where it is very dense, keep near-edge detail
    idx = np.unique(np.concatenate([np.arange(len(zs))[::2], np.where(np.abs(zs-40)<0.3)[0]]))
    maps[k] = dict(z=zs[idx].round(5).tolist(), r=ru.round(5).tolist(), s1=G[w][idx].T.round(2).tolist(), szz=Gz[w][idx].T.round(2).tolist())
out["maps"] = maps

# sensitivity + convergence (final state)
names = {"B":"Baseline: 1100-O feather, 0.1 mm tip, 75 µm bond, 6 mm taper",
         "B_tip03":"Feather tip 0.03 mm","B_tip25":"Feather tip 0.25 mm","B_tip50":"Feather tip 0.5 mm",
         "B_bond040":"Bond line 40 µm","B_bond150":"Bond line 150 µm",
         "B_ep06":"Epoxy modulus ×0.6","B_ep15":"Epoxy modulus ×1.5",
         "B_taper9":"Taper 9 mm (vs 6)","B_wall05":"Feather wall 0.5 mm (vs 1.0)",
         "B_y100":"Al yield 100 MPa (work-hardened)","B_e":"Al elastic (never yields)",
         "A":"Solid 6061-T6 collar (as originally posed)"}
sens = []
for k, nm in names.items():
    if not os.path.exists(f"res_{k}.pkl"): print("missing", k); continue
    r = load(k); f = r["field"]; X = f["X"][f["idx_glass"]]; st = r["stats"][-1]
    K, ac, d = fr(X, f["sig_glass"])
    row = dict(key=k, name=nm, bend=float(st["outer_s1_bend"]), s250=float(st["outer_s1_at_d"][4]),
               s50=float(st["outer_s1_at_d"][2]), tau_edge=float(st["bond_tau_edge"]), tau_deep=float(st["bond_tau_deep"]),
               peeq=float(st["al_peeq_body"]), **d)
    sens.append(row)
    print(f"{k:10s} bend {row['bend']:6.1f}  s@250 {row['s250']:6.1f}  s@50 {row['s50']:6.1f}  K25 {d['K25']:.2f}  ac(0.1/0.25) {d['ac_um_01']}/{d['ac_um_025']}  tau {row['tau_edge']:.0f}/{row['tau_deep']:.0f}  peeq {row['peeq']:.1e}")
out["sens"] = sens
conv = []
for k, h in [("A_h10",10),("A",5),("A_h025",2.5),("B_h10",10),("B",5),("B_h025",2.5)]:
    r = load(k); f = r["field"]; X = f["X"][f["idx_glass"]]; st = r["stats"][-1]
    K, ac, d = fr(X, f["sig_glass"])
    conv.append(dict(case=k[0], h_um=h, nodes=r["nn"], s25=float(st["outer_s1_at_d"][1]), s100=float(st["outer_s1_at_d"][3]),
                     s250=float(st["outer_s1_at_d"][4]), bend=float(st["outer_s1_bend"]), K25=d["K25"], ac_um=d["ac_um_01"], ac_um_025=d["ac_um_025"]))
    print("conv", conv[-1])
out["conv"] = conv
json.dump(out, open("results.json","w"), separators=(",",":"))
print("json size MB", os.path.getsize("results.json")/1e6)
