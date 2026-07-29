# SaLSPIA: Spectrally Adaptive LSPIA for B-Spline Fitting

MATLAB R2021b reference implementation accompanying:

> Xingxuan Peng and Yutong Li, “Spectral Adaptivity for Iterative
> Least-Squares B-Spline Curve and Surface Fitting.”

This release is synchronized with the revised manuscript. It contains the
SaLSPIA curve and tensor-product surface implementations, all comparison
methods used in the paper (LSPIA, ALSPIA, MLSPIA, NmLSPIA, and
LSPIA-Lin2018 where applicable), formal manuscript drivers, input data, and
generated CSV/LaTeX results.

## Requirements

- MATLAB R2021b.
- No external MATLAB toolbox is required.
- CPU time depends on system load. Compare the experiment configuration,
  iteration count (`IT`), and terminal relative error (`E_inf`) before
  comparing timings.

## Quick start

From the archive root:

```matlab
setup_paths
run(fullfile('experiments', 'run_table3_curve_comparison.m'))
```

Each timing driver performs five repeated solver runs and reports the mean and
sample standard deviation. The helper is `utils/run_with_avg_cpu.m`.

## Reproducibility map

### Tables

| Manuscript table | Content | Formal driver |
| --- | --- | --- |
| Table 1 | Ranks and effective condition numbers | `experiments/run_table1_condition_numbers.m` |
| Table 2 | Sensitivity to `M` (with sigma and delta checks) | `experiments/run_table2_parameter_sensitivity_curves.m`, `experiments/run_table2_parameter_sensitivity_surface.m` |
| Table 3 | Full-rank curve comparison, Ex. 4.1-4.3 | `experiments/run_table3_curve_comparison.m` |
| Table 4 | Theoretical-rate check, Ex. 4.3 | `experiments/run_table4_theoretical_rate.m` |
| Table 5 | Full-rank surface comparison, Ex. 4.4-4.6 | `experiments/run_table5_surface_comparison.m`, `experiments/run_table5_volcano_comparison.m` |
| Table 6 | Representative geometric residuals | `experiments/run_table6_geometric_residuals.m` |
| Table 7 | Missing-data curves and Franke surface | `experiments/run_table7_missing_data.m` |
| Table 8 | Spectral-weight ablation, Ex. 4.2, 4.3, 4.7 | `experiments/run_table8_ablation.m` |

`run_table7_missing_data.m` is the single formal entry for all four
rank-deficient cases. It runs the three missing-curve cases and the
`101 x 101` Franke case with `41 x 41` controls, then writes:

```text
results/table7_missing_data/table7_missing_data.csv
results/table7_missing_data/table7_missing_data.tex
```

For Ex. 4.8, the observed matrix has size `8627 x 1681`, MATLAB default
numerical rank `1639`, and 42 inactive columns.

### Figures

| Manuscript figure | Content | Source |
| --- | --- | --- |
| Figure 1 | Geometric interpretation of SaLSPIA | Manuscript schematic; not code-generated |
| Figure 2 | Reindeer fit, Ex. 4.3 | `run_table3_curve_comparison.m` |
| Figure 3 | Dini and Peaks fits, Ex. 4.4-4.5 | `run_table5_surface_comparison.m` |
| Figure 4 | Volcano iterations, Ex. 4.6 | `run_table5_volcano_comparison.m` |
| Figures 5-7 | Missing-data curves, Ex. 4.1, 4.2, 4.7 | `run_table7_missing_curves.m` |
| Figure 8 | Franke surface with four local holes, Ex. 4.8 | `run_table7_missing_surface.m` |
| Figure 9 | Ablation convergence, Ex. 4.7 | `generate_figure9_ablation.m` |

## SaLSPIA defaults and safeguards

All formal drivers use:

```matlab
salspia_params.c = 1e-4;
salspia_params.M = 10;
salspia_params.delta = 1e-8;
salspia_params.eps_saf = 1e-30;
```

Both `algorithms/salspia.m` and `algorithms/salspia_surf.m` implement:

- the initial `rho0 <= eps_saf` return;
- the iteration stop when `rho_k <= eps_saf`;
- the `b_k`/`c_k` safeguard;
- the stable positive-root formula separated by the signs of `Aq` and `Bq`.

## Directory layout

```text
algorithms/   SaLSPIA and comparison methods
utils/        B-spline, parameterization, residual, and timing helpers
experiments/  Formal drivers named by manuscript table/figure number
data/         Input data; provenance is documented in DATA_SOURCES.md
results/      Representative generated CSV, LaTeX, MAT, and figure outputs
```

Development-only scripts, unused surface-ablation code, obsolete
table-numbering drivers, and stale result files are intentionally excluded
from this release.

## Data

The prepared archive contains every input file needed by its formal drivers, so
no external non-bundled data are computationally required. See
`DATA_SOURCES.md` for provenance and redistribution notes. In particular,
redistribution permission must be confirmed before third-party curve files are
included in any public release.

## Citation

Public repository:
[https://github.com/Lily11898/SaLSPIA-Reproducibility](https://github.com/Lily11898/SaLSPIA-Reproducibility).
Machine-readable citation metadata are provided in `CITATION.cff`.

## License

The MATLAB code is released under the MIT License. Data files may have separate
terms; see `DATA_SOURCES.md`.
