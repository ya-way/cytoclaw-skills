---
name: sc-marker-extractor
description: Extract cluster-level differential expression markers from clustered scRNA-seq data and return structured JSON for each cluster. Use when the user asks for marker genes, differential expression by cluster, or top genes for cell type annotation.
---

# sc-marker-extractor

Compute differential markers per cluster and return a machine-friendly JSON object.

## Inputs

- `--input-h5ad`: clustered data path
- `--cluster-key`: cluster label column (`leiden` by default)
- `--top-n`: number of top markers per cluster
- `--out`: output JSON path

## Output Contract

Return JSON (stdout and written to `--out`):

- keys: cluster IDs
- values: ordered top marker gene lists

## Command

```bash
python scripts/extract_markers.py \
  --input-h5ad /abs/path/to/clustered.h5ad \
  --cluster-key leiden \
  --top-n 20 \
  --out /abs/path/to/markers.json
```
