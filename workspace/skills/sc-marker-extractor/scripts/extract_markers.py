#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Extract top marker genes per cluster.")
    parser.add_argument("--input-h5ad", required=True)
    parser.add_argument("--cluster-key", default="leiden")
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    import scanpy as sc

    in_h5ad = Path(args.input_h5ad).expanduser().resolve()
    out_json = Path(args.out).expanduser().resolve()
    out_json.parent.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(str(in_h5ad))
    if args.cluster_key not in adata.obs.columns:
        raise KeyError(f"Cluster key not found: {args.cluster_key}")

    sc.tl.rank_genes_groups(
        adata,
        groupby=args.cluster_key,
        method="wilcoxon",
        use_raw=False,
        pts=True,
    )

    cluster_series = adata.obs[args.cluster_key]
    if isinstance(cluster_series.dtype, pd.CategoricalDtype):
        groups = [str(x) for x in cluster_series.cat.categories]
    else:
        groups = sorted(map(str, np.unique(cluster_series.astype(str))))

    names = adata.uns["rank_genes_groups"]["names"]
    result = {}
    for g in groups:
        ranked = names[g][: args.top_n]
        result[str(g)] = [str(x) for x in ranked]

    out_json.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False))


if __name__ == "__main__":
    main()
