---
name: sc-dynamic-clustering
description: Run single-cell HVG selection, PCA/UMAP, and graph clustering with tunable n_pcs and resolution, then summarize cluster sizes. Use when the user asks for scRNA-seq clustering, T-cell subpopulation separation, Leiden/Louvain tuning, or iterative resolution search.
---

# sc-dynamic-clustering

Perform dimensionality reduction and graph clustering on cleaned single-cell data, with adjustable parameters for iterative search.

## Inputs

- `--input-h5ad`: cleaned dataset path (typically `qc_filtered.h5ad`)
- `--n-pcs`: principal components count
- `--resolution`: clustering resolution
- `--out-dir`: output folder

## Output Contract

Return JSON (stdout) with:

- `clustered_h5ad`
- `cluster_key` (`leiden` or fallback `louvain`)
- `n_clusters`
- `cluster_sizes` dictionary

Side effects in `--out-dir`:

- `clustered.h5ad` — clustered dataset with UMAP embeddings
- `clustering_summary.json` — same JSON persisted to disk

## Command

```bash
python scripts/dynamic_clustering.py \
  --input-h5ad /abs/path/to/qc_filtered.h5ad \
  --n-pcs 30 \
  --resolution 0.8 \
  --out-dir /abs/path/to/run_artifacts
```

## Notes

- Prefer Leiden; fallback to Louvain if Leiden dependencies are unavailable.
- This skill is intended for repeated calls during parameter search.
