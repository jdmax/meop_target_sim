# ³He Two-Cell Exchange

Interactive model of a polarized ³He target with a warm MEOP pumping cell (~300 K) and a cold target cell (~5 K). The cells are joined by glass transfer tubes, filled to the operating pressure at temperature, and then valved off from the supply.

The physics runs in Python ([Pyodide](https://pyodide.org) 0.26.4, loaded from the jsDelivr CDN) directly in the browser. Nothing is computed on a server.

It is published at **<https://jdmax.github.io/meop_target_sim/>**.

## What it models

- **Tube temperature profiles.** Glass-wall and gas conduction along each leg, including the enthalpy carried by any net gas flow. Clamped sections (cold clamp, riser "engine", intercept) are isothermal; the heat each one removes or supplies is reported.
- **Gas exchange by diffusion.** Polarization diffuses at uniform pressure with conductance `G = A / ∫ ds/(n D)`. Because `n·D ∝ T^0.68` regardless of pressure, the cold end dominates.
- **Thermosiphon loop.** Buoyancy from the riser/downcomer density difference is balanced against laminar friction. The temperature profiles are re-solved for each trial flow. The riser engine, a clamped section near the cold end, can be held at its own temperature by a heater, held by the first stage instead (it then sits at the first-stage temperature, and the stage supplies its heat, so no heater is needed), or turned off, leaving the downcomer intercept to drive the loop alone.
- **Polarization dynamics.** An exact solution of the two-cell rate equations. Pumping comes from a sealed cell measured at the operating pressure and field: a fit to `M(t) = Mss (1 − e^(−t/Tb))`, and `τp`, its relaxation with the discharge on and the laser off. The fit is read as fluxes of polarized atoms, not as rates. With `Ns` atoms in the sealed cell, it gives the optical production at zero polarization, `Q0 = Ns·Mss/Tb` atoms/s, the discharge relaxation `Rdis = Ns/τp`, and the limit `Plim = Mss/(Γop·Tb)`, where `Γop = 1/Tb − 1/τp`. At high field the discharge sits in a thin layer at the walls, and only metastables there absorb the pump light. So production and discharge relaxation scale with the absorbed photon flux (`s_light`) and the discharge volume (`s_dis`), relative to the sealed cell, and not with the pumping-cell volume:

  ```
  Np dPp/dt = Q0'(1 − Pp/Plim) − R' Pp − G(Pp − Pt)     Q0' = s_light·Q0,  R' = s_dis·Rdis (+ Np/T1,dark)
  dPt/dt    = (G/Nt)(Pp − Pt) − Pt/τt
  ```

  The first equation has the sealed-cell form `dPp/dt = (Mss' − Pp)/Tb' − (G/Np)(Pp − Pt)`, with `1/Tb' = (Q0'/Plim + R')/Np`, `Mss' = Q0'/(Q0'/Plim + R')`, and `τp' = Np/R'` once the laser is shuttered. In steady state the target reaches `Pt = Mss' / (1 + a + b)`, where `a = τx/T1t` (exchange too slow, τx = Nt/G) and `b = Nt·Mss'/(Q0'·T1t)` (MEOP can't keep up with target relaxation). With `s_light = s_dis` and no dark relaxation, `Mss' = Mss` and `b` doesn't depend on `Np`, so the steady state doesn't depend on the size of the pumping cell; that changes only the build-up time, `t_ex` and the probe signal. Optionally, `s_light` can come from the absorbed pump power with and without the discharge, and the page then reports the photon efficiency the sealed-cell fit implies, `Q0` per absorbed photon, which cannot exceed 2. Choosing **Pumping scales with: Cell volume** restores the older reading, `Mss' = Mss` and `Tb' = Tb`, in which production grows with the pumping-cell volume.
- **Reading the target from the pumping cell.** After the build-up, the laser is shuttered and the pumping cell is zeroed. The discharge stays on so the optical polarimeter can watch polarized gas come back from the target. The page shows this in the terms of [Maxwell, Epstein & Milner, arXiv:1412.6167](https://arxiv.org/abs/1412.6167): `N = Np + Nt`, the exchange time `t_ex = Np·Nt/(N·G)`, the decay constants `τs` and `τl` (eq. 4), and the fit amplitudes `a_s`, `a_l`, `a_st` and `a_lt` (eq. 3), both with and without zeroing. Once the short term has died out, `Mp/Mt = a_l/a_lt`.

- **Heat loads on the two cooler stages** (the **Heat loads** tab). The first stage, at the shield temperature, takes the intercept, radiation from 300 K onto the shield, and any other loads you enter. The second stage, at the target temperature, takes the conduction into the target cell and cold clamps, an engine clamp that has to remove heat, radiation from the shield onto the cold surfaces, the direct line of sight down each tube bore, and its own other loads. Radiation is `σεA(T_hot⁴ − T_cold⁴)` with an effective emissivity for each pair of surfaces. Clamps that need heating are listed as heaters, not cooler loads. In loop mode, the heat the circulating gas carries is shown separately from conduction: the tubes are solved a second time with no flow, and the difference in each clamp's heat is the gas cooled on its way down (at the intercept and the cold clamps) and rewarmed on its way up (by the engine and the warm end).
- **Gas temperature against the wall (a diagnostic).** The model takes the gas at the wall temperature. The **Heat loads** tab also estimates how far the flowing gas lags the wall: it exchanges heat with the model's wall profile as fully developed laminar flow does (Nusselt number 3.66), catching up over a length `λ = ṁ c_p / (Nu π k_gas)` that is about 1 mm at 300 K and 1–1.5 cm at 5 K. The tab plots wall and gas along both tubes and gives the gas temperature at each clamp, where it reaches the target cell, and the heat it carries past the intercept; a warning appears when the lag is large. The estimate holds the wall profile fixed, so it is not fed back into the solution.
- **The cryocooler (PT420 capacity map).** The stage temperatures come from Cryomech's PT420-RM reference capacity chart, read off at its 30 points and interpolated bilinearly in the two loads. The first stage floats at whatever temperature the two loads give. A trim heater holds the second stage at the target set point; if the loads alone keep it warmer, the target runs at that temperature, and the whole model uses it. The tube solution depends on the target temperature, so it is re-solved with secant steps until the two agree. A performance setting derates the chart; the certified rating, 1.8 W at 4.2 K with 50 W at 45 K, is a little below the reference curve. With the map off, both temperatures are fixed and the map gives each stage's capacity there. The **Heat loads** tab plots the operating point on the chart; hover over it to read the loads that give any pair of stage temperatures.

The **Formulas** tab lays these out as eleven steps in the order the model runs them. Each step gives its equations, what they apply to, how the result is reached, and the numbers they produce for the current settings.

**Not modeled:** the temperature drop between the second stage and the target cell, radiation exchanged with the tube walls (only the direct line of sight down each bore is counted), relaxation in the tubes, bends, 3-D flow in the cells, and quantum corrections to ³He transport below about 20 K. In the pumping cell: radial structure (the model assumes diffusion keeps the cell uniform, and no pumping or discharge relaxation outside the discharge beyond the optional dark `T1`), optical depth or saturation of the ring beam, any dependence of the photon efficiency on polarization beyond the linear sealed-cell form, and line shifts from field inhomogeneity along a longer cell. The borosilicate conductivity table, the default cell volumes and the sealed-cell values are approximate placeholders.

## Transfer tube joint stress

A separate page, **<https://jdmax.github.io/meop_target_sim/glass_joint/>**, models the thermal stress in the borosilicate transfer tube where it is bonded into the aluminum target cell, from 292 K down to 2 K. It compares a solid aluminum collar with the planned 1100-O feather seal and checks the glass against fracture. The model and how to rerun it are described in [glass_joint/README.md](glass_joint/README.md).

## Saving figures, settings and data

**Save / load** at the top of the page opens a panel with a title and notes field, a live preview, and four actions:

- **Download figure (PNG).** Choose any of the cell schematic, the two polarization charts, the PT420 map, a settings table and a results table. The figure is drawn on a white background at 2800 px wide, with the title, notes, date, model version and the page address in its header. With only the schematic and charts selected, it comes out close to 16:9 for a slide.
- **Save settings (JSON).** Every input, plus the key results and both sets of curves. **Load settings…** reads the file back into the page. After loading, the page says whether the current model reproduces the results recorded in the file, so a changed model is easy to spot. Files saved before the PT420 map load with the map off, so they keep their fixed temperatures. Files saved before the choice of pumping scaling load with pumping scaled by the cell volume, as it was then.
- **Save curves (CSV).** Columns `phase,t_s,P_pump,P_target`, with the settings in `#` comment lines at the top. Read it with `pandas.read_csv(path, comment="#")`.

## References

- J. Maxwell and R. Milner, "A Concept for Polarized ³He Targets for High Luminosity Scattering Experiments in High Magnetic Field Environments," [arXiv:2106.06606](https://arxiv.org/abs/2106.06606).
- J.D. Maxwell, C.S. Epstein and R.G. Milner, "Diffusive Transfer of Polarized ³He Gas through Depolarizing Magnetic Gradients," [arXiv:1412.6167](https://arxiv.org/abs/1412.6167).

## Files

| File | Purpose |
|---|---|
| `index.html` | The page: controls, schematic, the Polarization, Heat loads and Formulas tabs |
| `he3_model.py` | The model. The page fetches this file and runs it in Pyodide |
| `he3_model.js` | Line-for-line JavaScript port, used for instant previews while a slider is dragged and as a fallback if Python can't start |
| `glass_joint/` | Transfer-tube joint stress analysis: finite-element model, results and its page (`glass_joint/index.html`) |
| `.github/workflows/pages.yml` | Deploys the simulator files and the joint-stress page to GitHub Pages on each push to `main` |

## Editing the model

Edit `he3_model.py` and push. The page calls `run(params_json)` and expects the same JSON keys back.

After each Python run, the page compares the result against the JavaScript port. While they agree, it uses the JavaScript port for instant slider previews. If they differ, the status line at the bottom of the Formulas tab says so and the page waits for Python instead. To restore instant previews, mirror the change in `he3_model.js`. If the change alters an equation, update its step in the `FX` list in `index.html`. `python tests/test_pumping.py` checks the pumping model and, with Node.js on the PATH, that the two ports agree.

## Run locally

The page reads `he3_model.py` with `fetch`, which browsers block for `file://` pages. Serve the folder instead:

```
python3 -m http.server 8000
```

Then open http://localhost:8000.

## Author
Written in 2026 by J. Maxwell (https://orcid.org/0000-0003-2710-4646).

## License

MIT License - see [LICENSE](LICENSE) file for details.
