# 1000 Genomes pairwise Hudson FST — reference values

## Source

Privé F, Aschard H, Carmi S, Folkersen L, Hoggart C, O'Reilly PF, Vilhjálmsson BJ.
**"Portability of 245 polygenic scores when derived from the UK Biobank and applied to 9 ancestry groups from the same cohort."**
*The American Journal of Human Genetics*, 2022;109(1):12–23.

- DOI: 10.1016/j.ajhg.2021.11.008
- PMC: [PMC8764121](https://pmc.ncbi.nlm.nih.gov/articles/PMC8764121/)

## File

`fst_prive_2022.csv` — 315 pairs across 26 populations from 1000 Genomes Phase 3, in long format:

```
pop1,pop2,fst_hudson
ACB,ASW,0.0020
...
```

## Coverage

- 26 populations: LWK, ESN, YRI, ACB, ASW, GWD, MSL, JPT, CHB, CHS, CDX, KHV, GIH, PJL, BEB, ITU, STU, PEL, MXL, CLM, PUR, FIN, CEU, GBR, IBS, TSI.
- 315 unique pairs (out of 325 possible). The 10 within-South-Asian pairs (GIH/PJL/BEB/ITU/STU × each other) are not included in the source matrix.
- Estimator: Hudson (`Fst_HUDSON`), as computed by the Privé et al. pipeline on dense 1KG variants.

## Use in this project

These values serve as the calibration target for our internal Hudson FST estimates: see `calibrate_and_fit_fst.py` and `pca/plots/fst_vs_reference_*.png`.
