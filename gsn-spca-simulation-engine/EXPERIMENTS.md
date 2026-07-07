# Simulation Experiment Set Summary

This file briefly explains the 10 simulation parameter sets used by the package.

All experiments share the same planted structure idea:

- `PC1`: 15 target genes from the first three planted target clusters.
- `PC2`: 15 target genes from the last three planted target clusters.
- `NOISE-LOW`: 20 fixed noise genes inside the first 50 genes.
- `NOISE-HIGH`: added noise genes outside the first 50 genes, used only when `n_genes > 50`.

The four methods are compared on the same generated input for each parameter set:

- `ESPCA`
- `DM-ESPCA`
- `AWGE-ESPCA`
- `GSN-SPCA` / `ESPCA-LEE`

## Parameter Meanings

- `n_genes`: total number of genes in the simulated matrix.
- `d1`: signal strength/scale for PC1.
- `d2`: signal strength/scale for PC2.
- `c`: noise strength parameter.
- `p_target`: target-cluster edge density.
- `p_noise`: noise-edge density.
- `total_edges`: total generated edge count.
- `gsn_total_groups`: total GSN groups after standardization.

## 10 Parameter Sets

| config_id | family | Main purpose | n_genes | d1 | d2 | c | p_target | p_noise | total_edges | gsn_total_groups |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `baseline_current` | baseline | Baseline setting. Uses 50 genes, current signal/noise strength, and current target/noise edge density. | 50 | 10 | 5 | 5 | 0.9 | 0.7 | 195 | 86 |
| `density_60_40` | density | Lower-density network setting. Target density is 0.6 and noise density is 0.4. | 50 | 10 | 5 | 5 | 0.6 | 0.4 | 125 | 58 |
| `density_70_30` | density | Medium-low density setting. Target density is 0.7 and noise density is 0.3. | 50 | 10 | 5 | 5 | 0.7 | 0.3 | 114 | 42 |
| `density_80_20` | density | Medium-high target density with lower noise density. | 50 | 10 | 5 | 5 | 0.8 | 0.2 | 87 | 33 |
| `density_90_10` | density | High target density and very low noise density. Tests whether methods benefit from a cleaner network. | 50 | 10 | 5 | 5 | 0.9 | 0.1 | 79 | 24 |
| `noise_c10_current_density` | noise_strength | Same density as baseline, but noise strength `c` is increased from 5 to 10. | 50 | 10 | 5 | 10 | 0.9 | 0.7 | 195 | 86 |
| `noise_c20_current_density` | noise_strength | Same density as baseline, but noise strength `c` is increased from 5 to 20. | 50 | 10 | 5 | 20 | 0.9 | 0.7 | 195 | 86 |
| `genes_100_current_density` | gene_count | Increases total genes to 100 while keeping the same planted 50-gene structure and current density. Adds outside-50 noise genes. | 100 | 10 | 5 | 5 | 0.9 | 0.7 | 491 | 297 |
| `genes_1000_current_density` | gene_count | Scales total genes to 1000. Tests robustness when many outside-50 noise genes are present. | 1000 | 10 | 5 | 5 | 0.9 | 0.7 | 6535 | 4046 |
| `genes_10000_current_density` | gene_count | Large-scale gene-count stress test with 10000 genes. Tests whether methods still recover the planted PC1/PC2 structure under very large noise space. | 10000 | 10 | 5 | 5 | 0.9 | 0.7 | 66517 | 42044 |

## How To Read The Families

### Baseline

`baseline_current` is the reference experiment. Other experiments should be interpreted relative to this setting.

### Density Experiments

The `density_*` experiments keep `n_genes = 50`, `d1 = 10`, `d2 = 5`, and `c = 5`, but change `p_target` and `p_noise`.

They mainly test how each method behaves when the edge network becomes cleaner or sparser.

### Noise Strength Experiments

The `noise_c10_current_density` and `noise_c20_current_density` experiments keep the same edge density as the baseline, but increase the noise strength `c`.

They mainly test how robust each method is when noise magnitude becomes stronger.

### Gene Count Experiments

The `genes_100_*`, `genes_1000_*`, and `genes_10000_*` experiments keep the same planted PC1/PC2 target structure in the first 50 genes, but increase the total number of genes.

They mainly test scalability and false-positive behavior in the outside-50 noise region (`NOISE-HIGH`).
