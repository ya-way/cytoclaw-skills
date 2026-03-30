#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

import numpy as np


def _load_adata(input_path: Path):
    import scanpy as sc

    if input_path.is_file() and input_path.suffix == ".h5":
        adata = sc.read_10x_h5(str(input_path))
    elif input_path.is_dir():
        adata = sc.read_10x_mtx(str(input_path), var_names="gene_symbols", cache=False)
    else:
        raise ValueError(f"Unsupported input path: {input_path}")
    adata.var_names_make_unique()
    return adata


def _to_dense_vector(x):
    if hasattr(x, "A1"):
        return x.A1
    if hasattr(x, "toarray"):
        return np.asarray(x.toarray()).ravel()
    return np.asarray(x).ravel()


def _compute_mito_pct(adata):
    mito_mask = adata.var_names.str.upper().str.startswith("MT-")
    if mito_mask.sum() == 0:
        return np.zeros(adata.n_obs, dtype=float)
    mito_counts = _to_dense_vector(adata.X[:, mito_mask].sum(axis=1))
    total_counts = _to_dense_vector(adata.X.sum(axis=1))
    return np.divide(
        mito_counts,
        np.maximum(total_counts, 1e-12),
        out=np.zeros_like(mito_counts, dtype=float),
        where=total_counts > 0,
    ) * 100.0


def _mark_doublets(adata):
    total_counts = _to_dense_vector(adata.X.sum(axis=1))
    n_genes = _to_dense_vector((adata.X > 0).sum(axis=1))
    scrublet_used = False

    try:
        import scrublet as scr

        scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.06)
        _scores, preds = scrub.scrub_doublets()
        is_doublet = np.asarray(preds, dtype=bool)
        scrublet_used = True
    except Exception:
        # Conservative fallback when Scrublet is not available.
        count_cut = np.quantile(total_counts, 0.995)
        gene_cut = np.quantile(n_genes, 0.995)
        is_doublet = (total_counts >= count_cut) & (n_genes >= gene_cut)

    return is_doublet, scrublet_used


def main():
    parser = argparse.ArgumentParser(description="Adaptive scRNA-seq QC filtering with doublet removal.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--max-mito-percent", required=True, type=float)
    parser.add_argument("--min-genes", required=True, type=int)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    input_path = Path(args.input).expanduser().resolve()
    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    adata = _load_adata(input_path)
    cells_before = int(adata.n_obs)

    mito_pct = _compute_mito_pct(adata)
    n_genes = _to_dense_vector((adata.X > 0).sum(axis=1))
    keep_qc = (mito_pct <= args.max_mito_percent) & (n_genes >= args.min_genes)

    adata_qc = adata[keep_qc].copy()
    cells_after_qc_filters = int(adata_qc.n_obs)

    is_doublet, scrublet_used = _mark_doublets(adata_qc)
    adata_clean = adata_qc[~is_doublet].copy()
    cells_after_doublet_removal = int(adata_clean.n_obs)
    doublet_rate = float(is_doublet.mean()) if len(is_doublet) else 0.0

    out_h5ad = out_dir / "qc_filtered.h5ad"
    adata_clean.write_h5ad(str(out_h5ad))

    result = {
        "input": str(input_path),
        "qc_filtered_h5ad": str(out_h5ad),
        "cells_before": cells_before,
        "cells_after_qc_filters": cells_after_qc_filters,
        "cells_after_doublet_removal": cells_after_doublet_removal,
        "doublet_rate": doublet_rate,
        "scrublet_used": scrublet_used,
        "params": {
            "max_mito_percent": args.max_mito_percent,
            "min_genes": args.min_genes,
        },
    }
    (out_dir / "adaptive_qc_summary.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(result, ensure_ascii=False))


if __name__ == "__main__":
    main()
