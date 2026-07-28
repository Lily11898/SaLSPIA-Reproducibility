# Data sources and redistribution notes

## Reproducibility status

The prepared archive contains all input files used by the formal manuscript
drivers. Therefore:

> No external non-bundled data are required to run the archived experiments.

The source and redistribution status of every bundled file are listed below.

## Analytically generated data

The blob curve, helix curve, Dini surface, Peaks surface, and Franke surface are
generated directly inside the MATLAB drivers. They do not require data files or
external downloads.

## Maungawhau / Mt Eden LiDAR data

Bundled file:

- `data/maungawhau_hr.csv`

Used by:

- `experiments/run_table5_volcano_comparison.m`
- `experiments/run_table6_geometric_residuals.m`

Source:

- Toitu Te Whenua Land Information New Zealand (LINZ), Auckland LiDAR 1 m DEM
  (2013): https://data.linz.govt.nz/layer/53405-auckland-lidar-1m-dem-2013/

License:

- Creative Commons Attribution 4.0 International:
  https://creativecommons.org/licenses/by/4.0/

Attribution:

> Contains data sourced from Toitu Te Whenua Land Information New Zealand and
> licensed for reuse under CC BY 4.0.

The bundled CSV is a reformatted grid used in the manuscript experiments.

## Reindeer contour

Bundled file:

- `data/cur_data deer`

Used by:

- `experiments/run_table3_curve_comparison.m` (Example 4.3)
- `experiments/run_table4_theoretical_rate.m`
- `experiments/run_table6_geometric_residuals.m`
- parameter-sensitivity and ablation drivers

Upstream source:

- https://github.com/giacomoorsi/ProgressiveIterationApproximation

Redistribution note:

- The upstream repository page does not currently display an explicit
  redistribution license for this data file. The authors must confirm
  redistribution permission before publishing the prepared public archive.
- If permission is not obtained, remove this file from the public release and
  instruct users to obtain an authorized copy and place it at the exact path
  `data/cur_data deer`. In that case, the manuscript and README must not state
  that no external non-bundled data are required.

## G-shaped loop curve

Bundled file:

- `data/s_loop_curve_data.txt`

Used by:

- `experiments/run_table7_missing_curves.m`
- `experiments/run_table6_geometric_residuals.m` (Example 4.7)
- `experiments/run_table7_missing_data.m`
- parameter-sensitivity and ablation drivers

Upstream source:

- https://github.com/XuejiaoYuan/LSPIA

Redistribution note:

- The upstream repository page does not currently display an explicit
  redistribution license for this data file. The authors must confirm
  redistribution permission before publishing the prepared public archive.
- If permission is not obtained, remove this file from the public release and
  instruct users to obtain an authorized copy and place it at the exact path
  `data/s_loop_curve_data.txt`. In that case, the manuscript and README must
  identify the dependent experiments as requiring external data.
