# CytoClaw — 飞书单细胞生物学 AI 分析平台

CytoClaw（细胞夹）是一个嵌入飞书工作流的单细胞生物信息学 AI 分析平台，基于 [OpenClaw](https://github.com/nicepkg/openclaw) 网关构建。覆盖 scRNA-seq、scATAC-seq、空间转录组、TCR/BCR 免疫组库、多组学整合、基因调控网络等细胞生物学核心分析场景，内置 **98 个专业 Skills**。

用户在飞书中发送数据文件，CytoClaw 自动完成分析并以卡片消息推送结果、HTML 报告和可复现 Notebook。

## 能力矩阵

| 模块 | Skills | 覆盖领域 |
|---|---|---|
| **核心分析流程** | 9 | scRNA-seq 全流程、HTML 报告生成、ipynb 交付、飞书交互人格 |
| **单细胞分析** | 20 | 预处理/QC/聚类/注释/轨迹推断/doublet/批次整合/scATAC/多模态/Perturb-seq |
| **空间转录组** | 11 | 数据加载/预处理/空间域识别/反卷积/细胞通讯/多组学/可视化 |
| **TCR/BCR 免疫组库** | 5 | MiXCR/VDJtools/scirpy/ImmCantation/组库可视化 |
| **基因调控网络** | 5 | SCENIC/共表达网络/差异网络/多组学 GRN/扰动模拟 |
| **组学工作流** | 13 | scRNA/ATAC/spatial/multiome/RNA-seq DE/通路/GRN/ChIP/甲基化/剪接/CRISPR |
| **通路与差异分析** | 11 | GO/GSEA/KEGG/Reactome/WikiPathways/DESeq2/edgeR/批次矫正 |
| **数据可视化** | 11 | matplotlib/seaborn/plotly/热图/火山图/多面板/网络图/配色 |
| **文献与数据库** | 13 | PubMed/bioRxiv/GEO/CellxGene Census/基因数据库/文献综述/科学写作 |
| **飞书集成** | 3 | 卡片消息/多维表格/在线表格 |
| **基础设施** | 2 | FastQC 测序质控、RiboClaw UI 可视化 |

## Demo 流程

```
用户（飞书）                          CytoClaw
    │                                    │
    ├─ 发送 h5ad 文件 ─────────────────→│
    │                                    ├─ "收到，8000+ 细胞，跑着了"
    │                                    ├─ 静默运行 QC + 降维 + 聚类
    │                                    ├─ 定量检测批次效应
    │                                    │
    │←── ⚠️ Leiden UMAP + 形态预警卡片 ──┤
    │←── 🔴 Sample_ID 批次图 + 实锤卡片 ─┤
    │                                    │
    ├─ "确实分两次上机的" ──────────────→│
    │                                    ├─ 运行 Harmony 矫正
    │←── ✅ 矫正后 UMAP + 绿色卡片 ─────┤
    │                                    │
    ├─ "配色有点丑" ───────────────────→│
    │                                    ├─ 切换莫兰迪色系
    │←── 🎨 Nature 级 Figure Set ────────┤
    │←── 📄 Report.html + Pipeline.ipynb ┤
```

## 快速部署

### 前置条件

- Linux / macOS / WSL2
- Python 3.12+
- Node.js 18+
- [OpenClaw](https://github.com/nicepkg/openclaw) 已安装并可运行

### Step 1：克隆本仓库

```bash
git clone <repo-url> cytoclaw-skills
cd cytoclaw-skills
```

### Step 2：创建 Python 单细胞环境

```bash
python3 -m venv ~/.openclaw/workspace/sc_venv
source ~/.openclaw/workspace/sc_venv/bin/activate
pip install -r requirements.txt
```

验证：
```bash
python -c "import scanpy; import harmonypy; print('scanpy', scanpy.__version__, '| harmonypy OK')"
```

### Step 3：部署 Skills 到 OpenClaw workspace

```bash
# 复制系统人格文件
cp workspace/IDENTITY.md ~/.openclaw/workspace/
cp workspace/SOUL.md ~/.openclaw/workspace/
cp workspace/TOOLS.md ~/.openclaw/workspace/

# 复制所有 skills
cp -r workspace/skills/* ~/.openclaw/workspace/skills/

# 安装 feishu-card 和 feishu-bitable 的 Node 依赖
cd ~/.openclaw/workspace/skills/feishu-card && npm install
cd ~/.openclaw/workspace/skills/feishu-bitable && npm install
```

### Step 4：配置 OpenClaw

1. 复制模板：
```bash
cp openclaw.example.json ~/.openclaw/openclaw.json
```

2. 填入你的密钥（编辑 `~/.openclaw/openclaw.json`）：
   - `<YOUR_OPENROUTER_API_KEY>` — [OpenRouter](https://openrouter.ai/) API key
   - `<YOUR_DEEPSEEK_API_KEY>` — [DeepSeek](https://platform.deepseek.com/) API key
   - `<YOUR_BRAVE_SEARCH_API_KEY>` — [Brave Search](https://brave.com/search/api/) API key
   - `<YOUR_FEISHU_APP_ID>` / `<YOUR_FEISHU_APP_SECRET>` — 飞书开放平台应用凭证
   - `<GENERATE_A_RANDOM_TOKEN>` — 随机生成网关 token（`openssl rand -hex 24`）

### Step 5：配置飞书应用

1. 前往 [飞书开放平台](https://open.feishu.cn/) 创建企业自建应用
2. 开通权限：`im:message`、`im:message:send_as_bot`、`im:resource`、`sheets:spreadsheet`、`bitable:app`
3. 配置事件订阅 URL：`http://<你的服务器>:18789/feishu/event`
4. 启用机器人能力

### Step 6：启动

```bash
openclaw start
```

在飞书中 @ CytoClaw 或私聊发送 `.h5ad` 文件即可开始分析。

## 目录结构

```
cytoclaw-skills/
├── README.md                              ← 本文件
├── CUSTOMIZATION.md                       ← 人格/参数/配色定制指南
├── requirements.txt                       ← Python 依赖
├── openclaw.example.json                  ← 脱敏配置模板
└── workspace/
    ├── IDENTITY.md                        ← CytoClaw 身份定义
    ├── SOUL.md                            ← 行为准则与分析协议
    ├── TOOLS.md                           ← 环境路径与配色规范
    └── skills/                            ← 98 个专业 Skills
        ├── cytoclaw-scrna-workflow/        ← 核心 scRNA-seq 全流程编排
        ├── cytoclaw-analysis-report/       ← HTML 报告生成模板
        ├── cytoclaw-reproducible-notebook/ ← ipynb 可复现 Pipeline 交付
        ├── cytoclaw-feishu-interaction/    ← 飞书交互人格与卡片美学
        ├── feishu-card/                   ← 飞书卡片消息（JS）
        ├── feishu-bitable/                ← 飞书多维表格
        ├── feishu-sheets-skill/           ← 飞书表格
        ├── openclaw-fastqc/               ← FastQC 质控
        ├── riboclaw-visualization/        ← UI 可视化输出
        ├── bio-single-cell-*/             ← 14 个单细胞分析模块
        ├── bio-spatial-transcriptomics-*/ ← 11 个空间转录组模块
        ├── bio-tcr-bcr-analysis-*/        ← 5 个免疫组库模块
        ├── bio-gene-regulatory-networks-*/← 5 个基因调控网络模块
        ├── bio-wf-*/                      ← 13 个组学工作流
        ├── bio-pathway-analysis-*/        ← 6 个通路分析模块
        ├── bio-differential-expression-*/ ← 5 个差异表达模块
        ├── bio-data-visualization-*/      ← 6 个生物数据可视化模块
        ├── labclaw-viz-*/                 ← 5 个可视化工具
        ├── labclaw-lit-*/                 ← 8 个文献与数据库
        └── ...                            ← 更多专业模块
```

## 核心架构

```
飞书用户消息
    ↓
OpenClaw Gateway (port 18789)
    ↓
Feishu Plugin (事件订阅 → 消息解析)
    ↓
LLM Agent (Seed 2.0 Lite / DeepSeek-V3)
    ├── 读取 IDENTITY.md + SOUL.md → 建立人格
    ├── 读取 TOOLS.md → 知道 Python 环境路径
    ├── 匹配 cytoclaw-scrna-workflow SKILL → 编排全流程
    ├── 调用 sc_venv/bin/python 运行 scanpy 脚本
    └── 通过 feishu-card skill 发送结果卡片
```

## 定制与调参

参见 [CUSTOMIZATION.md](CUSTOMIZATION.md)，涵盖：

- 人格风格调整（博后口吻 ↔ PI 正式 ↔ 命令行简洁）
- QC 阈值、批次效应检测灵敏度、聚类分辨率
- 配色方案切换（莫兰迪 / Nature 经典 / Viridis）
- 飞书卡片颜色语义、消息节奏
- 交付物开关（HTML 报告、ipynb notebook）
- 模型选择

## 注意事项

- **Python 环境路径**：TOOLS.md 和 SKILL.md 中硬编码了 `~/.openclaw/workspace/sc_venv`，请确保 venv 安装在此路径
- **图片发送**：feishu-card 需要飞书应用具有 `im:resource` 权限才能上传图片
- **模型选择**：默认使用 Seed 2.0 Lite（便宜且支持 vision），可替换为其他 OpenAI 兼容模型
- **批次效应检测阈值**：silhouette > 0.1 或 UMAP 中心距 > 2.0 触发预警（见 SOUL.md）
- **Skill 架构**：98 个 skills 分层协作——核心流程 skill 编排分析逻辑，领域 skill 提供专业知识，飞书 skill 处理消息投递
