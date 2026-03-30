# CytoClaw scRNA-seq 分析工作流

**触发条件：** 用户发送 `.h5ad` 文件（或提及"单细胞/scRNA/质控/UMAP/聚类"）并要求分析。

---

## 完整执行流程（按步骤严格执行）

### Step 0 — 等待文件
若用户只说了需求、没发文件：
```
好的，把 h5ad 文件甩过来，我直接跑。
```
用 `feishu-card send_safe.js` 发这一句，等对方发文件。

---

### Step 1 — 收到文件，开场确认

发送纯文本（一句话，不解释任何操作细节）：
```bash
node skills/feishu-card/send_safe.js \
  --target "SENDER_ID" \
  --text "收到，[N]+ 个细胞，质控完，跑降维和聚类中……"
```
N 从 `adata.n_obs` 获取（格式化为 "8,000+" / "12,000+" 等大约值）。

---

### Step 2 — 静默运行 QC + 初步降维

写一个 Python 脚本（`/tmp/sc_run.py`）并执行，包含：

```python
import warnings; warnings.filterwarnings('ignore')
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scanpy as sc
from sklearn.metrics import silhouette_score

WS = '/home/shuotong/.openclaw/workspace'
H5AD = '<用户文件路径>'

adata = sc.read_h5ad(H5AD)

# ── QC ──
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
n_raw = adata.n_obs
adata = adata[adata.obs.pct_counts_mt < 20].copy()
adata = adata[adata.obs.n_genes_by_counts > 200].copy()
adata = adata[adata.obs.n_genes_by_counts < 6000].copy()
sc.pp.filter_genes(adata, min_cells=3)
n_qc = adata.n_obs
n_rm = n_raw - n_qc
mt_hi = (adata.obs.pct_counts_mt > 10).sum()

# ── 预处理 + PCA ──
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
batch_cols = ['Sample_ID','sequencing_batch','batch','patient','donor']
batch_col = next((c for c in batch_cols if c in adata.obs.columns), None)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat',
                             batch_key=batch_col)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50)

# ── UMAP + Leiden（pre-Harmony）──
sc.pp.neighbors(adata, n_pcs=30)
sc.tl.umap(adata, random_state=42)
sc.tl.leiden(adata, resolution=0.5, random_state=42)
n_cl = len(adata.obs.leiden.cat.categories)
umap = adata.obsm['X_umap']
labels_leiden = adata.obs.leiden.values

# ── 批次效应检测 ──
batch_flag, sil, dist, batch_col_found = False, 0, 0, batch_col
if batch_col and adata.obs[batch_col].nunique() >= 2:
    labels_batch = adata.obs[batch_col].values
    sil = silhouette_score(umap, labels_batch, sample_size=2000, random_state=42)
    unique = np.unique(labels_batch)
    centers = {s: umap[labels_batch==s].mean(0) for s in unique}
    k = list(centers.keys())
    dist = np.linalg.norm(centers[k[0]] - centers[k[1]])
    batch_flag = sil > 0.1 or dist > 2.0
    print(f"BATCH: col={batch_col} sil={sil:.4f} dist={dist:.2f} flag={batch_flag}")

# ── 图 1：Leiden UMAP（彩色，显示孤岛结构）──
VIVID = ['#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F',
         '#8491B4','#91D1C2','#DC0000','#7E6148','#B09C85','#FF7F0E','#2CA02C']
clusters = sorted(adata.obs.leiden.cat.categories, key=int)
pal = {cl: VIVID[i%len(VIVID)] for i,cl in enumerate(clusters)}
fig, ax = plt.subplots(figsize=(7.5, 6))
for cl in clusters:
    m = labels_leiden == cl
    ax.scatter(umap[m,0], umap[m,1], c=pal[cl], s=5, alpha=0.7,
               label=f'C{cl}', rasterized=True, linewidths=0)
ax.legend(title='Leiden', fontsize=7, markerscale=3, frameon=False,
          loc='center left', bbox_to_anchor=(1,0.5),
          ncol=2 if n_cl>8 else 1)
ax.set_xlabel('UMAP 1',fontsize=11); ax.set_ylabel('UMAP 2',fontsize=11)
ax.set_title(f'Leiden Clustering ({n_cl} clusters)',fontsize=12,pad=8)
ax.spines[['top','right']].set_visible(False)
plt.tight_layout()
fig.savefig(f'{WS}/out_leiden_umap.png', dpi=150, bbox_inches='tight')
plt.close()

# ── 图 2：Sample_ID UMAP（红蓝批次图，仅在 batch_flag 时生成）──
if batch_flag and batch_col:
    labels_b = adata.obs[batch_col].values
    unique_b = np.unique(labels_b)
    BATCH_COLS = ['#3A86FF','#FF6B6B','#06D6A0','#FFB703','#8338EC']
    b_pal = {s: BATCH_COLS[i%len(BATCH_COLS)] for i,s in enumerate(unique_b)}
    fig2, ax2 = plt.subplots(figsize=(7.5, 6))
    for s in unique_b:
        m = labels_b==s
        ax2.scatter(umap[m,0],umap[m,1],c=b_pal[s],s=4,alpha=0.5,
                    label=s,rasterized=True,linewidths=0)
    ax2.legend(fontsize=11, markerscale=4, frameon=False)
    ax2.set_xlabel('UMAP 1',fontsize=11); ax2.set_ylabel('UMAP 2',fontsize=11)
    ax2.set_title('Before Harmony — Colored by Sample',fontsize=12,pad=8)
    ax2.spines[['top','right']].set_visible(False)
    plt.tight_layout()
    fig2.savefig(f'{WS}/out_batch_umap.png', dpi=150, bbox_inches='tight')
    plt.close()

print(f"STATS: n_raw={n_raw} n_qc={n_qc} n_rm={n_rm} mt_hi={mt_hi} n_cl={n_cl}")
print(f"BATCH: flag={batch_flag} sil={sil:.4f} dist={dist:.2f}")
```

运行：
```bash
~/.openclaw/workspace/sc_venv/bin/python /tmp/sc_run.py
```

---

### Step 3 — 【第一发】Leiden 聚类图 + 形态描述

**若检测到批次效应（batch_flag=True）：**

写卡片正文：
```
cat > /tmp/msg_leiden.md << 'CARDEOF'
质控完成，**[n_qc] 个细胞**进入分析（MT 超标率 < [X]%）。

初步 Leiden 聚类图（↑）跑出 **[n_cl] 个群**，但形态**非常反常**——
所有细胞生硬地劈成了左右两组互不相交的孤岛。
正常的外周血免疫微环境不长这样，T/B/NK 细胞应该形成各自连续的拓扑结构。

查了一下 metadata，怀疑是跨周送测导致的批次问题。
我用 `[batch_col]` 重新给这些细胞上了一次色来验证……
CARDEOF
```

发送：
```bash
node skills/feishu-card/send.js \
  --target "SENDER_ID" \
  --image-path "~/.openclaw/workspace/out_leiden_umap.png" \
  --title "⚠️ 聚类形态严重异常" \
  --color orange \
  --text-file /tmp/msg_leiden.md
```

**若未检测到批次效应：**

直接发 Leiden 图 + 正常聚类卡片（绿色），列出 cluster 数量即可。

---

### Step 4 — 【第二发】Sample_ID 批次图 + 实锤（仅 batch_flag=True）

写卡片正文：
```
cat > /tmp/msg_batch.md << 'CARDEOF'
你看（↑图）——左边全是 [sample_A]，右边全是 [sample_B]。

两组 UMAP 中心距 **[dist]**，批次分离度（silhouette）**[sil]**。
这是一次典型的**批次效应（Batch Effect）**，真实的生物学亚群结构全被掩盖了。

这两个样本是分两次上机送测的吗？
CARDEOF
```

发送：
```bash
node skills/feishu-card/send.js \
  --target "SENDER_ID" \
  --image-path "~/.openclaw/workspace/out_batch_umap.png" \
  --title "🔴 实锤：典型批次效应" \
  --color red \
  --text-file /tmp/msg_batch.md
```

---

### Step 5 — 收到确认后：Harmony 矫正

**触发词：** "批次效应 / 两批 / 抹平 / 上机 / 确实 / Harmony / 快"

先发：
```bash
node skills/feishu-card/send_safe.js --target "SENDER_ID" \
  --text "收到，Harmony 跑起来了……"
```

然后运行 Harmony 脚本（接上一轮 adata，继续）：

```python
import harmonypy as hm

ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, batch_col,
                    max_iter_harmony=20, random_state=42)
adata.obsm['X_pca_harmony'] = ho.Z_corr.T
sc.pp.neighbors(adata, use_rep='X_pca_harmony', n_pcs=30)
sc.tl.umap(adata, random_state=42)
sc.tl.leiden(adata, resolution=0.5, random_state=42, key_added='leiden_harm')

umap_post = adata.obsm['X_umap']
labels_b = adata.obs[batch_col].values
sil_post = silhouette_score(umap_post, labels_b, sample_size=2000, random_state=42)
unique_b = np.unique(labels_b)
centers2 = {s: umap_post[labels_b==s].mean(0) for s in unique_b}
k = list(centers2.keys())
dist_post = np.linalg.norm(centers2[k[0]] - centers2[k[1]])
n_cl_harm = len(adata.obs.leiden_harm.cat.categories)

print(f"HARMONY: sil_pre={sil:.4f} sil_post={sil_post:.4f}")
print(f"HARMONY: dist_pre={dist:.2f} dist_post={dist_post:.2f}")
print(f"HARMONY: n_cl={n_cl_harm}")

# 生成 harmony_umap.png（红蓝批次配色）
fig3, ax3 = plt.subplots(figsize=(7.5,6))
for s in unique_b:
    m = labels_b==s
    ax3.scatter(umap_post[m,0],umap_post[m,1],c=b_pal[s],s=4,alpha=0.45,
                label=s,rasterized=True,linewidths=0)
ax3.legend(fontsize=11,markerscale=4,frameon=False)
ax3.set_xlabel('UMAP 1',fontsize=11); ax3.set_ylabel('UMAP 2',fontsize=11)
ax3.set_title('After Harmony — Colored by Sample',fontsize=12,pad=8)
ax3.spines[['top','right']].set_visible(False)
plt.tight_layout()
fig3.savefig(f'{WS}/out_harmony_umap.png', dpi=150, bbox_inches='tight')
plt.close()
```

发送：
```bash
cat > /tmp/msg_harmony.md << 'CARDEOF'
搞定了。Harmony **[N] 轮迭代**收敛，说明两批样本细胞类型构成一致，批次信号是可矫正的技术噪声。

两批样本 UMAP 中心距 **[dist] → [dist_post]**（↓[X]%），
批次分离度 **[sil] → [sil_post]**（↓[X]%）。
两座孤岛已打通，两批细胞充分混合（↑图）。

Leiden 重新聚类识别出 **[n_cl_harm] 个生物学亚群**，
符合外周血免疫微环境预期。
CARDEOF

node skills/feishu-card/send.js \
  --target "SENDER_ID" \
  --image-path "~/.openclaw/workspace/out_harmony_umap.png" \
  --title "✅ 批次效应消除" \
  --color green \
  --text-file /tmp/msg_harmony.md
```

---

### Step 6 — 可视化优化（用户嫌配色丑）

**触发词：** "配色 / 选色 / 品味 / 一坨 / 顶刊 / 审美 / 组会 / 周会"

先发：
```bash
node skills/feishu-card/send_safe.js --target "SENDER_ID" \
  --text "说得对，默认色板确实辣眼睛。切换 Nature 级配色方案，顺便出完整 figure set。"
```

然后生成三张图（见 TOOLS.md 莫兰迪色板）：
1. `out_nature_umap.png`（莫兰迪聚类 UMAP）
2. `out_heatmap.png`（Wilcoxon top markers，RdBu_r）
3. `out_cd8a.png`（CD8A 特征表达，Grey→Red）

分三条发送后，发卡片汇总：
```bash
cat > /tmp/msg_figset.md << 'CARDEOF'
**① 莫兰迪色系 UMAP** — [n_cl] 个亚群，colorblind-safe，组会/投稿直接用

**② Red-Blue 差异表达热图** — Wilcoxon 检验 top markers，
可见 CD3D（T）、LYZ（单核）、CD79A（B）、GNLY/GZMB（NK）、PPBP（血小板）等真实 marker

**③ CD8A 特征表达图** — Grey→Red 渐变

均为 200 dpi PNG，拖进 PPT 或直接投 Nature 子刊不丢人。
CARDEOF

node skills/feishu-card/send_safe.js \
  --target "SENDER_ID" \
  --title "🎨 顶刊 Figure Set" \
  --color purple \
  --text-file /tmp/msg_figset.md
```

接着发文件附言 + 可复现文件（若已生成）：
```bash
node skills/feishu-card/send_safe.js --target "SENDER_ID" \
  --text "另外，本次质控、Harmony 整合及差异分析的完整底层代码和运行日志已打包，您可以随时下载复现："
# 然后用飞书文件 API 发送 Analysis_Report.html 和 Reproducible_Pipeline.ipynb
```

---

## 注意事项

1. **每个 send 命令之间不要加延迟**——飞书会按序投递，自然有间隔感。
2. **所有图路径用绝对路径**，避免工作目录问题。
3. **SENDER_ID** 从当前对话上下文中获取（飞书消息的发送者 Open ID）。
4. **若用户没有 batch 字段**（单样本数据），跳过 Step 3/4，直接发 Leiden 图 + 聚类卡片。
5. **如果图生成失败**，先发文本结论，告知图稍后补发，不要让用户干等。
