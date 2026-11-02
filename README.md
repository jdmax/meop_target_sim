# ³He Two-Cell Exchange

Interactive model of a polarized ³He target with a warm MEOP pumping cell (~300 K) and a cold target cell (~5 K). The cells are joined by glass transfer tubes, filled to the operating pressure at temperature, and then valved off from the supply.

The physics runs in Python ([Pyodide](https://pyodide.org) 0.26.4, loaded from the jsDelivr CDN) directly in the browser. Nothing is computed on a server.

It is published at **<https://jdmax.github.io/meop_target_sim/>**.

## What it models

- **Tube temperature profiles.** Glass-wall and gas conduction along each leg, including the enthalpy carried by any net gas flow. Clamped sections (cold clamp, riser "engine", intercept) are isothermal; the heat each one removes or supplies is reported.
- **Gas exchange by diffusion.** Polarization diffuses at uniform pressure with conductance `G = A / ∫ ds/(n D)`. Because `n·D ∝ T^0.68` regardless of pressure, the cold end dominates.
- **Thermosiphon loop.** Buoyancy from the riser/downcomer density difference is balanced against laminar friction. The temperature profiles are re-solved for each trial flow.
- **Polarization dynamics.** An exact solution of the two-cell rate equations. The pumping cell is described by the steady polarization `Mss` and build-up time `Tb` measured in a sealed cell at the operating pressure and field, from a fit to `M(t) = Mss (1 − e^(−t/Tb))`:

  ```
  dPp/dt = (Mss − Pp)/Tb − (G/Np)(Pp − Pt)        laser and discharge on
  dPp/dt = −Pp/τp        − (G/Np)(Pp − Pt)        laser shuttered, discharge on
  dPt/dt = (G/Nt)(Pp − Pt) − Pt/τt
  ```

  In steady state the target reaches `Pt = Mss / (1 + a + b)`, where `a = τx/T1t` (exchange too slow, τx = Nt/G) and `b = (Nt/Np)(Tb/T1t)` (MEOP can't keep up with target relaxation).
- **Reading the target from the pumping cell.** After the build-up, the laser is shuttered and the pumping cell is zeroed. The discharge stays on so the optical polarimeter can watch polarized gas come back from the target. The page shows this in the terms of [Maxwell, Epstein & Milner, arXiv:1412.6167](https://arxiv.org/abs/1412.6167): `N = Np + Nt`, the exchange time `t_ex = Np·Nt/(N·G)`, the decay constants `τs` and `τl` (eq. 4), and the fit amplitudes `a_s`, `a_l`, `a_st` and `a_lt` (eq. 3), both with and without zeroing. Once the short term has died out, `Mp/Mt = a_l/a_lt`.

- **Heat loads on the two cooler stages** (the **Heat loads** tab). The first stage, at the shield temperature, takes the intercept, radiation from 300 K onto the shield, and any other loads you enter. The second stage, at the target temperature, takes the conduction into the target cell and cold clamps, an engine clamp that has to remove heat, radiation from the shield onto the cold surfaces, the direct line of sight down each tube bore, and its own other loads. Radiation is `σεA(T_hot⁴ − T_cold⁴)` with an effective emissivity for each pair of surfaces. Each stage is compared with the cooler capacity, which defaults to the PT420's 55 W at 45 K and 2.0 W at 4.2 K. Clamps that need heating are listed as heaters, not cooler loads.

**Not modeled:** radiation exchanged with the tube walls (only the direct line of sight down each bore is counted), relaxation in the tubes, bends, 3-D flow in the cells, and quantum corrections to ³He transport below about 20 K. The borosilicate conductivity table and the default cell volumes and MEOP rates are approximate placeholders.

## Saving figures, settings and data

**Save / load** at the top of the page opens a panel with a title and notes field, a live preview, and four actions:

- **Download figure (PNG).** Choose any of the cell schematic, the two charts, a settings table and a results table. The figure is drawn on a white background at 2800 px wide, with the title, notes, date, model version and the page address in its header. With only the schematic and charts selected, it comes out close to 16:9 for a slide.
- **Save settings (JSON).** Every input, plus the key results and both sets of curves. **Load settings…** reads the file back into the page. After loading, the page says whether the current model reproduces the results recorded in the file, so a changed model is easy to spot.
- **Save curves (CSV).** Columns `phase,t_s,P_pump,P_target`, with the settings in `#` comment lines at the top. Read it with `pandas.read_csv(path, comment="#")`.

## References

- J. Maxwell and R. Milner, "A Concept for Polarized ³He Targets for High Luminosity Scattering Experiments in High Magnetic Field Environments," [arXiv:2106.06606](https://arxiv.org/abs/2106.06606).
- J.D. Maxwell, C.S. Epstein and R.G. Milner, "Diffusive Transfer of Polarized ³He Gas through Depolarizing Magnetic Gradients," [arXiv:1412.6167](https://arxiv.org/abs/1412.6167).

## Files

| File | Purpose |
|---|---|
| `index.html` | The page: controls, schematic, the Polarization and Heat loads tabs, source listing |
| `he3_model.py` | The model. The page fetches this file and runs it in Pyodide; its listing is shown on the page |
| `he3_model.js` | Line-for-line JavaScript port, used for instant previews while a slider is dragged and as a fallback if Python can't start |
| `.github/workflows/pages.yml` | Deploys the three site files to GitHub Pages on each push to `main` |

## Editing the model

Edit `he3_model.py` and push. The page calls `run(params_json)` and expects the same JSON keys back.

After each Python run, the page compares the result against the JavaScript port. While they agree, it uses the JavaScript port for instant slider previews. If they differ, the status line under "Model source" says so and the page waits for Python instead. To restore instant previews, mirror the change in `he3_model.js`.

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
