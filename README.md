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
  dPp/dt = (Mss − Pp)/Tb − (G/Np)(Pp − Pt)        discharge on
  dPp/dt = −Pp/T1p,off   − (G/Np)(Pp − Pt)        discharge off
  dPt/dt = (G/Nt)(Pp − Pt) − Pt/T1t
  ```

  In steady state the target reaches `Pt = Mss / (1 + a + b)`, where `a = τx/T1t` (exchange too slow, τx = Nt/G) and `b = (Nt/Np)(Tb/T1t)` (MEOP can't keep up with target relaxation).

**Not modeled:** radiation down the bore, relaxation in the tubes, bends, 3-D flow in the cells, and quantum corrections to ³He transport below about 20 K. The borosilicate conductivity table and the default cell volumes and MEOP rates are approximate placeholders.

## Files

| File | Purpose |
|---|---|
| `index.html` | The page: controls, schematic, polarization trace, source listing |
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
