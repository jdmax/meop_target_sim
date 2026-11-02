# Transfer Tube Cooldown Stress

Axisymmetric finite-element model of the borosilicate transfer tube (8 mm OD, 1 mm wall) where it is epoxy-bonded into the aluminum target cell, through a cooldown from 292 K to 2 K at −29 K/h. It compares a solid 6061-T6 collar with the planned soft-aluminum (1100-O) feather seal and checks the glass against fracture.

Published at **<https://jdmax.github.io/meop_target_sim/glass_joint/>**.

## Results in brief

| | Solid 6061 collar | 1100-O feather |
|---|---|---|
| Axial tension on the glass surface, 0.25 mm outside the bond edge, at 2 K | 98 MPa | 29 MPa |
| K for a 25 µm surface ring crack, as a fraction of K<sub>IC</sub> = 0.77 MPa√m | 1.6 | 0.64 |
| Worst ring crack (≥ 0.1 mm from the edge), K/K<sub>IC</sub> | 2.5 | 1.05 |
| Fraction of the 2 K surface stress already present at 77 K | 92 % | 91 % |

- The axial temperature gradient alone produces at most 0.85 MPa, at the two kinks in the temperature profile, and about 0.02 MPa through the middle of the free span.
- The load is a shell-bending effect at the end of the bond. Inside the bond the metal and epoxy squeeze the glass radially (0.36 % contraction mismatch). Where the bond ends, the wall bends back to its free radius over about √(Rt) ≈ 1.5 mm, putting the outer surface into axial tension.
- The feather works only because it yields. An elastic feather matches the solid collar, and a work-hardened one (100 MPa yield) raises the edge stress by about 70 %.
- A thinner tip, a thinner feather wall, a longer taper, a thinner bond line and a softer epoxy each lower the edge stress by 10–18 %. See the sensitivity table on the page.

## Model

- 9-node Lagrange quadrilaterals, 3 × 3 Gauss, axisymmetric small strain. Mesh graded to 5 µm at the bond edge (about 26,000 nodes); results agree within 2 % at 10, 5 and 2.5 µm.
- Thermal strain is the integrated contraction from 293 K (room-temperature cure). Each 30-minute snapshot is one Newton load step. The axial thermal time constant of the free span is about 10 minutes, so the tube is treated as quasi-static.
- Aluminum: J2 plasticity with linear isotropic hardening, radial return and a consistent tangent. Glass and epoxy are linear elastic.
- Temperature: the first 10 mm is held at 292 K; the 30 mm free span carries steady 1-D conduction with borosilicate k(T) through a Kirchhoff transform; the 10 mm bond is at the sink temperature.
- Fracture: long circumferential surface crack opened by σ<sub>zz</sub>, with K from Glinka's edge-crack weight function on the FE stress field through the wall. Crack planes within 0.1 mm of the epoxy corner are excluded, because a perfectly sharp elastic corner gives a singular field that a real fillet and yielding epoxy do not.

### Material data

| Material | Data |
|---|---|
| Pyrex 7740 | E 63 GPa, ν 0.20, ΔL/L(293→4 K) = 5.6 × 10⁻⁴ (Ekin), α(293 K) = 3.0 × 10⁻⁶ K⁻¹, K<sub>IC</sub> 0.77 MPa√m |
| 6061-T6, 1100 Al | NIST cryogenic fits for ΔL/L(T) and E(T); yield 300 MPa (6061-T6), 40 MPa with 1.5 GPa hardening (1100-O) |
| Stycast 2850FT | ΔL/L = 4.78 × 10⁻³ at 77 K (FNAL TS-SSC-90-060), 5.08 × 10⁻³ at 4 K (Ekin); E 7 → 13 GPa on cooling (assumed) |

The Pyrex and Stycast curves between those points are single-Einstein shapes fitted to the tabulated values. The 2 K results depend only on the totals.

**Not modeled:** the epoxy fillet, epoxy yielding or crazing, radiation onto the free span, and gas in the bore.

## Files

| File | Purpose |
|---|---|
| `index.html` | The published page (built; do not edit by hand) |
| `page_template.html` | Page source, with a `__DATA__` placeholder for the results |
| `results.json` | Results embedded in the page, written by `analyze.py` |
| `materials.py` | Temperature-dependent properties |
| `fe_model.py` | Mesh, element routines, constitutive update and Newton solver |
| `fracture.py` | Stress-intensity calculation for surface ring cracks |
| `run_cases.py` | Full 30-minute cooldown sweeps: `A` collar, `B` feather, `G` glass alone |
| `sens.py` | Mesh-convergence and feather-design variants, run straight to 2 K |
| `analyze.py` | Fracture metrics and assembly of `results.json` |
| `build_page.py` | Embeds `results.json` in the template to produce `index.html` |

## Reproduce

Needs only NumPy and SciPy.

```
cd glass_joint
python3 run_cases.py A      # about 1–2 min each
python3 run_cases.py B
python3 run_cases.py G
python3 sens.py A_h10 A_h025 B_h10 B_h025 B_e B_y100 B_tip03 B_tip25 B_tip50 \
                B_bond040 B_bond150 B_ep06 B_ep15 B_taper9 B_wall05
python3 analyze.py
python3 build_page.py
```

The runs write `res_*.pkl` files, which are not committed. To model a different joint, add a `Case` in `sens.py`: `tip`, `taper`, `wall`, `bond`, `al_yield`, `ep_scale` and `kind` (`"collar"` or `"feather"`) set the geometry and materials.

## Sources

- [NIST cryogenic material properties, 6061-T6 aluminum](https://trc.nist.gov/cryogenics/materials/6061%20Aluminum/6061_T6Aluminum_rev.htm)
- J. Ekin, *Experimental Techniques for Low-Temperature Measurements*, thermal contraction table ([cryogenic reference tables](https://jwyoo.unist.ac.kr/main/Cryogenic_Reference_Tables.pdf))
- [Fermilab TS-SSC-90-060, thermal contraction of Stycast 2850FT](https://lss.fnal.gov/archive/other/ssc_notes/TS-SSC-90-060.pdf)
- [MIT material property database, Pyrex](https://www.mit.edu/~6.777/matprops/pyrex.htm)
- G. Glinka and G. Shen, "Universal features of weight functions for cracks in mode I," *Eng. Fract. Mech.* 40, 1135 (1991)
