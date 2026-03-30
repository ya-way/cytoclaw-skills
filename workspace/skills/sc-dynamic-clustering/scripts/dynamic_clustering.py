#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

import numpy as np


def _cluster_sizes(adata, key):
    counts = adata.obs[key].value_counts().sort_index()
    return {str(k): int(v) for k, v in counts.items()}


def main():
    parser = argparse.ArgumentParser(description="Dynamic clustering for scRNA-seq.")
    parser.add_argument("--input-h5ad", required=True)
    parser.add_argument("--n-pcs", type=int, required=True)
    parser.add_argument("--resolution", type=float, required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    import scanpy as sc

    in_h5ad = Path(args.input_h5ad).expanduser().resolve()
    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(str(in_h5ad))

    # Standard preprocessing pipeline for clustering.
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=max(args.n_pcs, 2), svd_solver="arpack")
    sc.pp.neighbors(adata, n_pcs=args.n_pcs)
    sc.tl.umap(adata)

    cluster_key = "leiden"
    try:
        sc.tl.leiden(adata, resolution=args.resolution, key_added=cluster_key)
    except Exception:
        cluster_key = "louvain"
        sc.tl.louvain(adata, resolution=args.resolution, key_added=cluster_key)

    cluster_sizes = _cluster_sizes(adata, cluster_key)
    result = {
        "input_h5ad": str(in_h5ad),
        "clustered_h5ad": str(out_dir / "clustered.h5ad"),
        "cluster_key": cluster_key,
        "n_clusters": int(len(cluster_sizes)),
        "cluster_sizes": cluster_sizes,
        "params": {"n_pcs": args.n_pcs, "resolution": args.resolution},
    }

    adata.write_h5ad(result["clustered_h5ad"])
    (out_dir / "clustering_summary.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(result, ensure_ascii=False))


if __name__ == "__main__":
    main()
