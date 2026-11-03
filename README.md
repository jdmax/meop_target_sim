# ³He Two-Cell Exchange

Interactive model of a polarized ³He target with a warm MEOP pumping cell (~300 K) and a cold target cell (~5 K). The cells are joined by glass transfer tubes, filled to the operating pressure at temperature, and then valved off from the supply.

The physics runs in Python ([Pyodide](https://pyodide.org) 0.26.4, loaded from the jsDelivr CDN) directly in the browser. Nothing is computed on a server.

It is published at **<https://jdmax.github.io/meop_target_sim/>**.

## What it models

- **Tube temperature profiles.** Glass-wall and gas conduction along each leg, including the enthalpy carried by any net gas flow. Clamped sections (cold clamp, riser "engine", intercept) are isothermal; the heat each one removes or supplies is reported.
- **Gas exchange by diffusion.** Polarization diffuses at uniform pressure with conductance `G = A / ∫ ds/(n D)`. Because `n·D ∝ T^0.68` regardless of pressure, the cold end dominates.
- **Thermosiphon loop.** Buoyancy from the riser/downcomer density difference is balanced against laminar friction. The temperature profiles are re-solved for each trial flow. The riser engine, a clamped section near the cold end, can be held at its own temperature by a heater, held by the first stage instead (it then sits at the first-stage temperature, and the stage supplies its heat, so no heater is needed), or turned off, leaving the downcomer intercept to drive the loop alone.
- **Polarization dynamics.** An exact solution of the two-cell rate equations. The pumping cell is described by the steady polarization `Mss` and build-up time `Tb` measured in a sealed cell at the operating pressure and field, from a fit to `M(t) = Mss (1 − e^(−t/Tb))`:

  ```
  dPp/dt = (Mss − Pp)/Tb − (G/Np)(Pp − Pt)        laser and discharge on
  dPp/dt = −Pp/τp        − (G/Np)(Pp − Pt)        laser shuttered, discharge on
  dPt/dt = (G/Nt)(Pp − Pt) − Pt/τt
  ```

  In steady state the target reaches `Pt = Mss / (1 + a + b)`, where `a = τx/T1t` (exchange too slow, τx = Nt/G) and `b = (Nt/Np)(Tb/T1t)` (MEOP can't keep up with target relaxation).
- **Reading the target from the pumping cell.** After the build-up, the laser is shuttered and the pumping cell is zeroed. The discharge stays on so the optical polarimeter can watch polarized gas come back from the target. The page shows this in the terms of [Maxwell, Epstein & Milner, arXiv:1412.6167](https://arxiv.org/abs/1412.6167): `N = Np + Nt`, the exchange time `t_ex = Np·Nt/(N·G)`, the decay constants `τs` and `τl` (eq. 4), and the fit amplitudes `a_s`, `a_l`, `a_st` and `a_lt` (eq. 3), both with and without zeroing. Once the short term has died out, `Mp/Mt = a_l/a_lt`.

- **Heat loads on the two cooler stages** (the **Heat loads** tab). The first stage, at the shield temperature, takes the intercept, radiation from 300 K onto the shield, and any other loads you enter. The second stage, at the target temperature, takes the conduction into the target cell and cold clamps, an engine clamp that has to remove heat, radiation from the shield onto the cold surfaces, the direct line of sight down each tube bore, and its own other loads. Radiation is `σεA(T_hot⁴ − T_cold⁴)` with an effective emissivity for each pair of surfaces. Clamps that need heating are listed as heaters, not cooler loads. In loop mode, the heat the circulating gas carries is shown separately from conduction: the tubes are solved a second time with no flow, and the difference in each clamp's heat is the gas cooled on its way down (at the intercept and the cold clamps) and rewarmed on its way up (by the engine and the warm end).
- **Gas temperature against the wall (a diagnostic).** The model takes the gas at the wall temperature. The **Heat loads** tab also estimates how far the flowing gas lags the wall: it exchanges heat with the model's wall profile as fully developed laminar flow does (Nusselt number 3.66), catching up over a length `λ = ṁ c_p / (Nu π k_gas)` that is about 1 mm at 300 K and 1–1.5 cm at 5 K. The tab plots wall and gas along both tubes and gives the gas temperature at each clamp, where it reaches the target cell, and the heat it carries past the intercept; a warning appears when the lag is large. The estimate holds the wall profile fixed, so it is not fed back into the solution.
- **The cryocooler (PT420 capacity map).** The stage temperatures come from Cryomech's PT420-RM reference capacity chart, read off at its 30 points and interpolated bilinearly in the two loads. The first stage floats at whatever temperature the two loads give. A trim heater holds the second stage at the target set point; if the loads alone keep it warmer, the target runs at that temperature, and the whole model uses it. The tube solution depends on the target temperature, so it is re-solved with secant steps until the two agree. A performance setting derates the chart; the certified rating, 1.8 W at 4.2 K with 50 W at 45 K, is a little below the reference curve. With the map off, both temperatures are fixed and the map gives each stage's capacity there. The **Heat loads** tab plots the operating point on the chart; hover over it to read the loads that give any pair of stage temperatures.

- **Discharge schedule** (a separate page, `discharge_schedule.html`). A strong discharge pumps fast to a lower polarization and a weak one slowly to a higher one. Each discharge level is entered as its sealed-cell `Mss` and `Tb` (plus the discharge-on `T1p`, which splits `1/Tb` into `Γop + 1/T1p` and gives `Plim`), so its net production per pumping-cell atom is the straight line `f_d(P) = (Mss,d − P)/Tb,d`. Because the target holds ~120× the pumping cell's atoms and the system is MEOP-limited (`Np/(Tb G)` ≪ 1), the cells share nearly one polarization, and the time to a goal, `∫ dP/(dP/dt)`, is shortest when each `dP` uses the level with the larger `f_d(P)`. Target relaxation subtracts the same at a given P for every level, so the crossover `P* = (Mss,H/Tb,H − Mss,L/Tb,L)/(1/Tb,H − 1/Tb,L)` doesn't depend on the losses (in whole-system terms, `P* = (k_H P_H − k_L P_L)/(k_H − k_L)`). The page runs one level, a two-level switch (at `P*`, at a set pumping-cell polarization, or at a set time) or follows the upper envelope of `f_d(Pp)` with optional hysteresis and minimum dwell. Each constant-level stretch uses the exact two-cell solution, with switches located by bisection on it. It reports the time to a goal (absolute, or a share of the best single-level steady state) for each level alone and for the schedule, and warns when `Np/(Tb G)` exceeds 0.2 for a level, since the exchange then starts to limit the build-up. The default levels are placeholders to be replaced with measurements.

The **Formulas** tab lays these out as eleven steps in the order the model runs them. Each step gives its equations, what they apply to, how the result is reached, and the numbers they produce for the current settings.

**Not modeled:** the temperature drop between the second stage and the target cell, radiation exchanged with the tube walls (only the direct line of sight down each bore is counted), relaxation in the tubes, bends, 3-D flow in the cells, and quantum corrections to ³He transport below about 20 K. The borosilicate conductivity table and the default cell volumes and MEOP rates are approximate placeholders.

## Transfer tube joint stress

A separate page, **<https://jdmax.github.io/meop_target_sim/glass_joint/>**, models the thermal stress in the borosilicate transfer tube where it is bonded into the aluminum target cell, from 292 K down to 2 K. It compares a solid aluminum collar with the planned 1100-O feather seal and checks the glass against fracture. The model and how to rerun it are described in [glass_joint/README.md](glass_joint/README.md).

## Saving figures, settings and data

**Save / load** at the top of the page opens a panel with a title and notes field, a live preview, and four actions:

- **Download figure (PNG).** Choose any of the cell schematic, the two polarization charts, the PT420 map, a settings table and a results table. The figure is drawn on a white background at 2800 px wide, with the title, notes, date, model version and the page address in its header. With only the schematic and charts selected, it comes out close to 16:9 for a slide.
- **Save settings (JSON).** Every input, plus the key results and both sets of curves. **Load settings…** reads the file back into the page. After loading, the page says whether the current model reproduces the results recorded in the file, so a changed model is easy to spot. Files saved before the PT420 map load with the map off, so they keep their fixed temperatures.
- **Save curves (CSV).** Columns `phase,t_s,P_pump,P_target`, with the settings in `#` comment lines at the top. Read it with `pandas.read_csv(path, comment="#")`.

## References

- J. Maxwell and R. Milner, "A Concept for Polarized ³He Targets for High Luminosity Scattering Experiments in High Magnetic Field Environments," [arXiv:2106.06606](https://arxiv.org/abs/2106.06606).
- J.D. Maxwell, C.S. Epstein and R.G. Milner, "Diffusive Transfer of Polarized ³He Gas through Depolarizing Magnetic Gradients," [arXiv:1412.6167](https://arxiv.org/abs/1412.6167).

## Files

| File | Purpose |
|---|---|
| `index.html` | The page: controls, schematic, the Polarization, Heat loads and Formulas tabs |
| `discharge_schedule.html` | The discharge-schedule page. It shares the cell, tube and cooler settings with `index.html` |
| `he3_model.py` | The model. Both pages fetch this file and run it in Pyodide: `run()` for the main page, `run_schedule()` for the schedule page |
| `he3_model.js` | Line-for-line JavaScript port, used for instant previews while a slider is dragged and as a fallback if Python can't start |
| `tests/` | Checks of the schedule model: reduced-model hand calculations, the single-level trace against `run()`, and Python/JavaScript parity (`python tests/test_schedule.py`, `node tests/test_schedule.js`) |
| `glass_joint/` | Transfer-tube joint stress analysis: finite-element model, results and its page (`glass_joint/index.html`) |
| `.github/workflows/pages.yml` | Deploys both simulator pages, the model files and the joint-stress page to GitHub Pages on each push to `main` |

## Editing the model

Edit `he3_model.py` and push. The main page calls `run(params_json)`, the schedule page `run_schedule(params_json)`, and each expects the same JSON keys back. Both take the tube, cooler and exchange results from `system()`.

After each Python run, each page compares the result against the JavaScript port. While they agree, it uses the JavaScript port for instant slider previews. If they differ, the status line at the bottom of the Formulas tab says so and the page waits for Python instead. To restore instant previews, mirror the change in `he3_model.js`. If the change alters an equation, update its step in the `FX` list in `index.html`.

## Run locally

The page reads `he3_model.py` with `fetch`, which browsers block for `file://` pages. Serve the folder instead:

```
python3 -m http.server 8000
```

Then open http://localhost:8000, or http://localhost:8000/discharge_schedule.html for the schedule page.

To check the schedule model, and that the JavaScript port matches it (this needs Node.js), run from the repository root:

```
python3 tests/test_schedule.py
node tests/test_schedule.js
```

## Author
Written in 2026 by J. Maxwell (https://orcid.org/0000-0003-2710-4646).

## License

MIT License - see [LICENSE](LICENSE) file for details.
