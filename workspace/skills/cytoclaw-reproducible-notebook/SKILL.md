---
name: cytoclaw-reproducible-notebook
description: |
  Generate a Jupyter notebook (.ipynb) containing the complete reproducible pipeline for the current analysis.
  Use after every completed analysis session. Always delivered alongside Analysis_Report.html.
  Triggered automatically at the end of the workflow, or when the user says "代码/pipeline/复现/notebook/ipynb".
---

# CytoClaw Reproducible Notebook Generator

## 触发条件

与 `Analysis_Report.html` 同时生成。每次完成分析后**必须**交付此文件。

## 生成方式

用 Python 的 `nbformat` 库创建 `.ipynb` 文件：

```python
import nbformat as nbf
import datetime

nb = nbf.v4.new_notebook()
nb.metadata.kernelspec = {
    "display_name": "Python 3",
    "language": "python",
    "name": "python3"
}
cells = []
```

若 `nbformat` 不可用，直接写 JSON 结构：

```python
import json

notebook = {
    "nbformat": 4, "nbformat_minor": 5,
    "metadata": {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.12"}
    },
    "cells": []
}

def md_cell(source):
    return {"cell_type": "markdown", "metadata": {}, "source": [source]}

def code_cell(source):
    return {"cell_type": "code", "metadata": {}, "source": [source],
            "execution_count": None, "outputs": []}
```

## Notebook 标准结构（按此顺序）

### Cell 1 — 标题（Markdown）

```markdown
# CytoClaw scRNA-seq 分析 Pipeline

**生成时间：** {datetime}
**数据文件：** {filename}
**环境：** scanpy {version} · harmonypy {version} · Python 3.12

> 本 notebook 完整复现了 CytoClaw 的自动化分析流程。
> 在同一 Python 环境下顺序执行所有 cell 即可复现全部结果。
```

### Cell 2 — 环境安装（Code）

```python
# 如果尚未安装，取消注释运行：
# !pip install scanpy==1.12 harmonypy==0.0.10 leidenalg scikit-learn matplotlib
```

### Cell 3 — Import（Code）

```python
import warnings; warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import scanpy as sc
from sklearn.metrics import silhouette_score
import harmonypy as hm
```

### Cell 4 — 加载数据（Code）

```python
adata = sc.read_h5ad('数据文件路径')
print(f"原始数据: {adata.n_obs} cells × {adata.n_vars} genes")
```

### Cell 5 — 质控过滤（Code）

包含 MT 标记、QC 指标计算、阈值过滤。使用与实际分析完全一致的参数。

### Cell 6 — QC 结果（Markdown）

```markdown
## 质控结果
- 原始细胞: {n_raw}
- 过滤后: {n_qc}（移除 {n_rm} 个低质量细胞）
- MT 超标率 > 10% 的细胞: {mt_hi}
```

### Cell 7 — 预处理 + PCA（Code）

Normalize → log1p → HVG → Scale → PCA。参数与分析一致。

### Cell 8 — UMAP + Leiden（Code）

```python
sc.pp.neighbors(adata, n_pcs=30)
sc.tl.umap(adata, random_state=42)
sc.tl.leiden(adata, resolution=0.5, random_state=42)
```

### Cell 9 — 批次效应检测（Code）

完整的 silhouette + 中心距计算代码。

### Cell 10 — 批次效应结论（Markdown）

用人话说明检测结果（有/无批次效应，关键数值）。

### Cell 11 — Harmony 矫正（Code，条件性）

若存在批次效应，包含 Harmony 代码。若无，本 cell 为空并注释说明。

### Cell 12 — 可视化（Code）

包含所有图表生成代码：
- Leiden UMAP（莫兰迪色板）
- Batch UMAP（Before/After Harmony）
- 差异表达热图
- 特征表达图

每张图后接 `plt.show()`（notebook 环境下直接显示）。

### Cell 13 — Marker 基因分析（Code）

```python
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False)
```

### Cell 14 — 保存结果（Code）

```python
adata.write('analysis_result.h5ad')
print("分析完成，结果已保存。")
```

### Cell 15 — 环境信息（Code）

```python
sc.logging.print_versions()
```

## 关键约束

1. **随机种子固定** — 所有 `random_state=42`，确保可复现
2. **参数与实际分析一致** — 从运行时变量回填，不用默认值
3. **每个 code cell 可独立运行** — 不依赖未定义的中间变量（除了 `adata`）
4. **Markdown cell 用中文** — 解释每步的目的和结论
5. **不包含飞书发送代码** — notebook 只关注分析逻辑
6. **图表直接 `plt.show()`** — 不保存文件，在 notebook 内联显示

## 输出路径

```
~/.openclaw/workspace/Reproducible_Pipeline.ipynb
```

## 飞书发送附言

交付时固定用此措辞：

```
另外，本次质控、Harmony 整合及差异分析的完整底层代码和运行日志已打包，您可以随时下载复现：
```

然后通过飞书文件 API 发送 `Analysis_Report.html` 和 `Reproducible_Pipeline.ipynb`。
