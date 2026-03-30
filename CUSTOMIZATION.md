# CytoClaw 定制指南

本文档说明如何通过修改 `SOUL.md`、`IDENTITY.md`、`TOOLS.md` 来调整 CytoClaw 的行为、人格和分析参数，**无需修改任何代码**。

---

## 1. 人格调整（IDENTITY.md + SOUL.md）

### 1.1 改变说话风格

打开 `workspace/IDENTITY.md`，修改 `Vibe` 字段：

```markdown
- **Vibe:** 像实验室里最靠谱的那个生信博后——专业但不端着，说人话不说论文摘要……
```

**调整方向举例：**

| 想要的效果 | 修改为 |
|---|---|
| 更严肃正式 | "像资深 PI 做学术汇报——严谨精确，结论有据可查" |
| 更活泼轻松 | "像隔壁工位的生信小伙伴——会用表情包，偶尔吐槽" |
| 更简洁冷静 | "像 Unix 命令行——只给结论和关键数值，不废话" |

### 1.2 改变语言

默认简体中文。如需切换：

```markdown
- **语言:** 永远用英文，除非用户主动切中文。
```

### 1.3 控制主动性

在 `SOUL.md` 的 `⑧ 专业主动性` 部分：

- **更主动**：保持现有设定（自动检测批次效应、自动出报告）
- **更被动**：改为 "等用户明确要求后再执行下一步，不要自动推进流程"
- **关闭审美干预**：删除 `⑥ 可视化标准` 中的"自动切换 Nature 级配色"段落

---

## 2. 分析参数调整（SOUL.md + TOOLS.md）

### 2.1 QC 过滤阈值

在 `SOUL.md` 的 `③ scRNA-seq 分析协议` 和 `TOOLS.md` 的 `QC 标准阈值` 部分：

```python
adata = adata[adata.obs.pct_counts_mt < 20].copy()    # ← MT 阈值
adata = adata[adata.obs.n_genes_by_counts > 200].copy()  # ← 最低基因数
adata = adata[adata.obs.n_genes_by_counts < 6000].copy() # ← 最高基因数
```

| 参数 | 默认值 | 调松（保留更多细胞） | 调紧（更严格过滤） |
|---|---|---|---|
| MT% 上限 | 20% | 25-30% | 10-15% |
| 最低基因数 | 200 | 100 | 500 |
| 最高基因数 | 6000 | 8000 | 4000 |

### 2.2 批次效应检测阈值

在 `SOUL.md` 第三步：

```python
batch_flag = sil > 0.1 or dist > 2.0
```

- **更敏感**（容易触发预警）：`sil > 0.05 or dist > 1.5`
- **更宽松**（减少误报）：`sil > 0.2 or dist > 3.0`

### 2.3 聚类分辨率

在 `cytoclaw-scrna-workflow/SKILL.md` 的 Step 2：

```python
sc.tl.leiden(adata, resolution=0.5, random_state=42)
```

- `0.3`：粗聚类（适合样本量小或追求大群分类）
- `0.5`：默认（适合大多数 PBMC 数据）
- `0.8-1.0`：细聚类（适合发现稀有亚群）

### 2.4 Batch 列名自动检测

在 `TOOLS.md` 和 skill 中，batch 列按优先级依次匹配：

```python
batch_cols = ['Sample_ID', 'sequencing_batch', 'batch', 'patient', 'donor']
```

如果你的数据用了不同的列名，把它加到列表前面即可。

---

## 3. 配色方案（TOOLS.md）

### 3.1 聚类 UMAP 色板

默认莫兰迪色系（9 色）：

```python
MORANDI = ['#8FAADC','#C5B4A0','#D4A5A5','#9DB4A0','#B8A9C9',
           '#E6C6A0','#A0C4C8','#C9B1D0','#D1C6A8']
```

替换为其他方案：

| 方案 | 色板 | 适用场景 |
|---|---|---|
| Nature 经典 | `['#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F','#8491B4','#91D1C2']` | 需要高对比度 |
| Cell 柔和 | `['#5B8FA8','#D4A574','#7BA68C','#C17C74','#8B7EB8','#A8C4D4']` | 偏学术优雅 |
| Viridis | 使用 `plt.cm.viridis` | colorblind-safe 首选 |

### 3.2 批次对比色

```python
COL_BATCH = {'Patient_A': '#3A86FF', 'Patient_B': '#FF6B6B'}
```

改为你实际的样本名和对应颜色。

### 3.3 特征表达渐变

```python
grey_red = LinearSegmentedColormap.from_list('gr', ['#EEEEEE','#CC0000'], N=256)
```

可替换为 Grey→Purple、Grey→Blue 等方案。

---

## 4. 飞书卡片风格（cytoclaw-feishu-interaction SKILL）

### 4.1 颜色语义

在 `skills/cytoclaw-feishu-interaction/SKILL.md` 中定义了卡片颜色与语义的对应关系。如需调整：

- 把 `orange` 预警改为 `red`（更强烈）
- 把 `purple` 成果展示改为 `green`（更克制）

### 4.2 消息节奏

默认是"分步推送"（先图后文、分多条发）。如需改为一次性推送所有结果，修改 `SOUL.md` 第四步为：

```
将所有图和结论合并为一条卡片消息发送。
```

---

## 5. 交付物控制

### 5.1 关闭 HTML 报告

如果不需要 `Analysis_Report.html`，在 `SOUL.md` 的 `⑥ 交付物` 部分删除相关条目。

### 5.2 关闭 Notebook 交付

同上，删除 `Reproducible_Pipeline.ipynb` 相关条目。

### 5.3 修改附言措辞

在 `SOUL.md` 和 `cytoclaw-reproducible-notebook/SKILL.md` 中修改固定附言文案。

---

## 6. 模型选择

在 `openclaw.json` 中修改 `agents.defaults.model.primary`：

| 模型 | 特点 | 推荐场景 |
|---|---|---|
| Seed 2.0 Lite | 便宜、快、支持 vision | 日常使用（默认） |
| DeepSeek-V3 | 中文理解强、推理稳 | 需要更好的中文对话 |
| GPT-4o Mini | 均衡、vision 强 | 需要读图分析 |

---

## 快速验证

修改完配置后，用以下流程验证效果：

1. 重启 OpenClaw：`openclaw restart`
2. 在飞书中发送一个 `.h5ad` 文件
3. 观察：
   - [ ] 开场确认消息是否符合预期语气
   - [ ] 聚类图是否使用了你设定的配色
   - [ ] 批次效应检测是否按你的阈值判定
   - [ ] 卡片颜色是否符合语义
   - [ ] HTML 报告和 ipynb 是否自动交付
