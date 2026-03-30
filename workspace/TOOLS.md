# TOOLS.md — 环境与工具

---

## Python 单细胞环境

**venv 路径：** `~/.openclaw/workspace/sc_venv`

```bash
# 激活后运行脚本
source ~/.openclaw/workspace/sc_venv/bin/activate && python your_script.py
# 或直接调用
~/.openclaw/workspace/sc_venv/bin/python your_script.py
```

**已安装关键包：** scanpy 1.12、anndata、harmonypy 0.0.10、leidenalg、matplotlib、numpy、scipy、scikit-learn

---

## 飞书消息发送（feishu-card skill）

工作目录为 `~/.openclaw/workspace`，所有 node 命令在此执行。

### 纯文本
```bash
node skills/feishu-card/send_safe.js \
  --target "TARGET_ID" \
  --text "消息内容"
```

### 带图片的卡片（主要用法）
```bash
# 1. 把卡片正文写入临时 md 文件（避免 shell 转义问题）
cat > /tmp/card_body.md << 'EOF'
卡片正文，支持 **加粗** 和 `等宽字体`。
EOF

# 2. 发送图片 + 卡片
node skills/feishu-card/send.js \
  --target "TARGET_ID" \
  --image-path "/absolute/path/to/figure.png" \
  --title "卡片标题" \
  --color orange \
  --text-file /tmp/card_body.md
```

**颜色选项：** `blue`（默认）/ `red` / `orange` / `green` / `purple` / `grey`

### 不含图片的带标题卡片
```bash
node skills/feishu-card/send_safe.js \
  --target "TARGET_ID" \
  --title "✅ 分析完成" \
  --color green \
  --text "两批样本中心距 3.50 → 1.30（↓63%）……"
```

### TARGET_ID 来源
- 用户 Open ID：`ou_` 开头（从对话上下文中获取当前发言者 ID）
- 群聊 Chat ID：`oc_` 开头

---

## scRNA-seq 标准工作流

### 必要 import
```python
import warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import scanpy as sc
from sklearn.metrics import silhouette_score
```

### QC 标准阈值（PBMC / 外周血）
```python
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
adata = adata[adata.obs.pct_counts_mt < 20].copy()   # 死细胞
adata = adata[adata.obs.n_genes_by_counts > 200].copy()
adata = adata[adata.obs.n_genes_by_counts < 6000].copy()
sc.pp.filter_genes(adata, min_cells=3)
```

### 批次效应定量检测
```python
umap = adata.obsm['X_umap']
# 自动找 batch 列（按优先级）
batch_cols = ['Sample_ID', 'sequencing_batch', 'batch', 'patient', 'donor']
batch_col = next((c for c in batch_cols if c in adata.obs.columns), None)

if batch_col and adata.obs[batch_col].nunique() >= 2:
    labels = adata.obs[batch_col].values
    sil = silhouette_score(umap, labels, sample_size=2000, random_state=42)
    unique = np.unique(labels)
    centers = {s: umap[labels==s].mean(0) for s in unique}
    k = list(centers.keys())
    dist = np.linalg.norm(centers[k[0]] - centers[k[1]])
    # 阈值判断
    batch_flag = sil > 0.1 or dist > 2.0
    print(f"batch_col={batch_col}  sil={sil:.4f}  dist={dist:.2f}  flag={batch_flag}")
```

### 默认图片配色

**莫兰迪聚类色板（9色）**
```python
MORANDI = ['#8FAADC','#C5B4A0','#D4A5A5','#9DB4A0','#B8A9C9',
           '#E6C6A0','#A0C4C8','#C9B1D0','#D1C6A8']
```

**批次对比**
```python
COL_BATCH = {'Patient_A': '#3A86FF', 'Patient_B': '#FF6B6B'}
# 若 batch 名称不同，用前两色即可
```

**特征表达**
```python
grey_red = LinearSegmentedColormap.from_list('gr', ['#EEEEEE','#CC0000'], N=256)
```

### 图片输出规范
```python
fig.savefig('/path/to/figure.png', dpi=150, bbox_inches='tight')
# 顶刊三件套用 dpi=200
```

---

## 文件位置约定

- 分析脚本：`~/.openclaw/workspace/`
- 生成图片：`~/.openclaw/workspace/`（绝对路径传给 feishu-card）
- 临时 md：`/tmp/`（feishu-card 卡片正文）
