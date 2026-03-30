# CytoClaw — 飞书单细胞分析 AI 助手

CytoClaw（细胞夹）是一个嵌入飞书工作流的 scRNA-seq AI 分析师，基于 [OpenClaw](https://github.com/nicepkg/openclaw) 网关构建。用户在飞书中发送 `.h5ad` 文件，CytoClaw 自动完成质控、降维聚类、批次效应检测与矫正、顶刊级可视化，并将结果以卡片消息推送回飞书。

## 能力概览

| Skill | 功能 |
|---|---|
| **cytoclaw-scrna-workflow** | 核心 scRNA-seq 全流程：QC → PCA → UMAP → Leiden → 批次效应检测 → Harmony 矫正 → 顶刊出图 |
| **feishu-card** | 飞书富文本卡片消息（支持图片、Markdown、颜色标题、按钮） |
| **feishu-bitable** | 飞书多维表格读写（记录分析任务） |
| **feishu-sheets-skill** | 飞书在线表格操作（导出 QC 统计表） |
| **openclaw-fastqc** | 上游 FASTQ 测序质控（FastQC + MultiQC） |
| **riboclaw-visualization** | 前端 UI 可视化链路 JSON 输出 |

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
├── requirements.txt                       ← Python 依赖
├── openclaw.example.json                  ← 脱敏配置模板
└── workspace/
    ├── IDENTITY.md                        ← CytoClaw 身份定义
    ├── SOUL.md                            ← 行为准则与分析协议
    ├── TOOLS.md                           ← 环境路径与配色规范
    └── skills/
        ├── cytoclaw-scrna-workflow/        ← scRNA-seq 全流程编排
        ├── feishu-card/                   ← 飞书卡片消息（JS）
        ├── feishu-bitable/                ← 飞书多维表格
        ├── feishu-sheets-skill/           ← 飞书表格
        ├── openclaw-fastqc/               ← FastQC 质控
        └── riboclaw-visualization/        ← UI 可视化输出
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

## 注意事项

- **Python 环境路径**：TOOLS.md 和 SKILL.md 中硬编码了 `~/.openclaw/workspace/sc_venv`，请确保 venv 安装在此路径
- **图片发送**：feishu-card 需要飞书应用具有 `im:resource` 权限才能上传图片
- **模型选择**：默认使用 Seed 2.0 Lite（便宜且支持 vision），可替换为其他 OpenAI 兼容模型
- **批次效应检测阈值**：silhouette > 0.1 或 UMAP 中心距 > 2.0 触发预警（见 SOUL.md）
