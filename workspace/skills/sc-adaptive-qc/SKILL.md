---
name: sc-adaptive-qc
description: Perform adaptive single-cell QC filtering with user-selected mitochondrial and gene thresholds, then run doublet removal and save qc_filtered.h5ad. Use when the user requests adaptive QC, mito filtering, min genes filtering, or doublet detection for scRNA-seq.
---

# sc-adaptive-qc

Execute threshold-based QC filtering plus doublet removal, then persist a clean intermediate dataset.

## Inputs

- `--input`: raw 10x path (directory or .h5)
- `--max-mito-percent`: chosen mitochondrial cutoff (for example 5 or 10)
- `--min-genes`: chosen minimum genes per cell
- `--out-dir`: output folder

## Output Contract

Return JSON (stdout) with:

- `qc_filtered_h5ad`
- `cells_before`
- `cells_after_qc_filters`
- `cells_after_doublet_removal`
- `doublet_rate`

Side effects in `--out-dir`:

- `qc_filtered.h5ad` — filtered dataset
- `adaptive_qc_summary.json` — same JSON persisted to disk

## Command

```bash
python scripts/adaptive_qc.py \
  --input /abs/path/to/raw \
  --max-mito-percent 10 \
  --min-genes 200 \
  --out-dir /abs/path/to/run_artifacts
```

## Notes

- Prefer Scrublet when available.
- If Scrublet is unavailable, fallback to a conservative heuristic and report fallback usage in JSON.
