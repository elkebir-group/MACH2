# MACH2 data

Clonal trees with observed location labelings of all instances used in the MACH2 paper, and the ground truth of the
simulated instances. The file formats are described in the [main README](../README.md#21-io-formats).
[`../analysis/run.ipynb`](../analysis/run.ipynb) generates the MACHINA and Metient inputs from these files and runs
MACH2, MACHINA and Metient with the settings of the paper.

```
data/
├── simulations/
│   ├── set1/m5/{mS,S,M,R}/<seed>.*                          first set:    40 instances
│   ├── set2/baseline/, set2/<parameter>=<value>/<seed>.*     second set:   70 instances
│   └── set3/number-of-locations=<n>/mutation-rate=<m>/driver-mutation-probability=<d>/
│            migration-rate=<g>/number-of-samples-per-location=<k>/migration-pattern=<p>/<seed>.*
│                                                             third set:  4860 instances
└── real/
    ├── lung/<patient>_<rank>.*     TRACERx, 126 patients, 266 instances
    ├── prostate/A<id>.*            Gundem et al. 2015, 10 patients
    ├── ovarian/<id>.*              McPherson et al. 2016, 7 patients (12 instances)
    └── breast/A<id>.*              Hoadley et al. 2016, 2 patients
```

## Files per instance

| File | Content | Simulations | Real |
|---|---|:-:|:-:|
| `<id>.mach2.tree` | clonal tree, one edge `parent child` per line | ✓ | ✓ |
| `<id>.mach2.labeling` | observed locations of each clone, `node loc1 loc2 …` | ✓ | ✓ |
| `<id>.gt.mach2.tree` | ground truth: refined tree (`X^L` = copy of clone `X` in location `L`) | ✓ | – |
| `<id>.gt.mach2.labeling` | ground truth: location of every refined-tree node | ✓ | – |
| `<id>.mutations.tsv` | number of mutations per clone, as given to Metient in the paper | – | ✓ |

## Number of instances

| Dataset | Instances | Breakdown |
|---|---:|---|
| First set of simulations (m5) | 40 | 4 migration patterns (mS, S, M, R) × 10 instances |
| Second set of simulations | 70 | 14 configurations (default baseline + 13 alternative settings of 6 parameters) × 5 replicates |
| Third set of simulations | 4860 | 972 configurations (3 × 3 × 3 × 3 × 3 × 4 parameter values) × 5 replicates |
| **Simulations, total** | **4970** | |
| Lung cancer (TRACERx) | 126 | highest-ranked clonal tree of 126 patients |
| Prostate cancer | 10 | 10 patients |
| Ovarian cancer | 12 | 7 patients; 5 with both ovaries as potential primary (5 × 2 + 2) |
| Breast cancer | 2 | 2 patients |
| **Real data, total** | **150** | 145 patients |
| Lung cancer, robustness analysis | +140 | trees ranked 2–5 of the 40 patients with ≥3 locations and multiple clonal trees |
