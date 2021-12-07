# scrapR
R package for importing and processing [`scraps`](https://github.com/rnabioco/scraps) reference and output

### read reference (if converting to gene-level psi values)
```
saf <- scrapR::parse_saf_pf("polyadb32.hg38.saf.gz")
```

### read output
```
mat <- scrapR::scraps_to_matrix("R2_counts.tsv.gz", pf = saf)
```

### insert into Seurat object
```
so <- scrapR::scraps_to_seurat("R2_counts.tsv.gz", pf = saf, so)
```
