# EasyMultiProfiler V3.0

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R Package](https://img.shields.io/badge/R-%E2%89%A54.0-blue)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-3.8%2B-blue)](https://www.python.org/)
[![Docker](https://img.shields.io/badge/Docker-Ready-blue)](https://www.docker.com/)

> **统一多组学数据分析平台 (Unified Multi-Omics Analysis Platform)**
> 
> R包 + Web界面 + Docker部署，一站式解决多组学数据分析需求

---

## 📋 目录

- [快速开始](#-快速开始)
- [功能特性](#-功能特性)
- [安装指南](#-安装指南)
- [使用教程](#-使用教程)
- [架构设计](#-架构设计)
- [API文档](#-api文档)
- [相关链接](#-相关链接)

---

## 🚀 快速开始

### 方式1: Docker一键部署（推荐）

```bash
# 克隆仓库
git clone https://github.com/xielab2017/EasyMultiProfiler-V3.git
cd EasyMultiProfiler-V3

# Docker启动
docker-compose up -d

# 访问 http://localhost:8080
```

### 方式2: 一键安装脚本

```bash
# 自动安装R包和Web环境
bash scripts/install.sh

# 启动服务
bash scripts/start.sh
```

### 方式3: R包单独使用

```r
# 安装R包
devtools::install_github("xielab2017/EasyMultiProfiler-V3", subdir = "r-package")

# 加载使用
library(EasyMultiProfiler)
```

---

## ✨ 功能特性

### 支持的组学类型

| 组学类型 | R函数 | Web界面 | 核心技术 |
|----------|-------|---------|----------|
| **RNA-seq** | `EMP_rnaseq_analysis()` | ✅ | DESeq2/edgeR/limma |
| **单细胞RNA-seq** | `EMP_scrnaseq_analysis()` | ✅ | Seurat, SingleR |
| **蛋白质组学** | `EMP_proteomics_analysis()` | ✅ | limma, DEP |
| **ChIP-seq** | `EMP_chipseq_analysis()` | ✅ | ChIPseeker, MACS2 |
| **CUT&Tag** | `EMP_cutntag_analysis()` | ✅ | ChIPseeker, SEACR |
| **CUT&RUN** | `EMP_cutnrun_analysis()` | ✅ | ChIPseeker, SEACR |
| **微生物组** | `EMP_microbiome_analysis()` | ✅ | phyloseq, vegan |
| **代谢组学** | `EMP_metabolome_analysis()` | ✅ | MetaboAnalystR |
| **多组学整合** | `EMP_multiomics_integration()` | ✅ | MOFA2, mixOmics |

### 分析流程

```
数据导入 → 质控 → 标准化 → 差异分析 → 富集分析 → 可视化 → 报告生成
```

### 可视化类型

- 📊 火山图、MA图、热图
- 🗺️ PCA、UMAP、t-SNE降维
- 📈 箱线图、小提琴图、散点图
- 🧬 基因组覆盖图、Peak注释图
- 🔗 网络图、桑基图

---

## 📦 安装指南

### 系统要求

| 组件 | 最低要求 | 推荐配置 |
|------|----------|----------|
| R | 4.0.0 | 4.3.0+ |
| Python | 3.8 | 3.11+ |
| Node.js | 16 | 18+ |
| 内存 | 8GB | 16GB+ |
| 磁盘 | 10GB | 50GB+ |

### Docker安装

```bash
# 1. 确保已安装Docker和docker-compose
# https://docs.docker.com/get-docker/

# 2. 克隆仓库
git clone https://github.com/xielab2017/EasyMultiProfiler-V3.git
cd EasyMultiProfiler-V3

# 3. 启动服务
docker-compose up -d

# 4. 查看日志
docker-compose logs -f

# 5. 停止服务
docker-compose down
```

### 本地安装

#### 步骤1: 安装R包

```r
# 安装devtools
install.packages("devtools")

# 安装EasyMultiProfiler V3
devtools::install_github(
    "xielab2017/EasyMultiProfiler-V3",
    subdir = "r-package",
    dependencies = TRUE
)
```

#### 步骤2: 安装Web后端

```bash
cd EasyMultiProfiler-V3/web/backend

# 创建虚拟环境
python3 -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# 安装依赖
pip install -r requirements.txt
```

#### 步骤3: 安装Web前端

```bash
cd ../frontend
npm install
npm run build
```

#### 步骤4: 启动服务

```bash
# 终端1: 启动后端
cd web/backend
python app.py

# 终端2: 启动前端 (开发模式)
cd web/frontend
npm start

# 或生产模式
# 前端构建后，Flask会自动服务静态文件
```

---

## 📖 使用教程

### R包使用示例

#### 单细胞RNA-seq分析

```r
library(EasyMultiProfiler)

# 读取数据
counts <- read.csv("counts.csv", row.names = 1)
metadata <- read.csv("metadata.csv", row.names = 1)

# 执行分析
result <- EMP_scrnaseq_analysis(
    counts = counts,
    metadata = metadata,
    params = list(
        qc = list(min_genes = 200, max_mt_percent = 5),
        clustering = list(resolution = 0.8),
        annotation = list(enable = TRUE)
    ),
    output_dir = "./scRNA_results"
)

# 查看结果
print(result)
DimPlot(result$seurat_object, reduction = "umap")
```

#### ChIP-seq分析

```r
# Peak分析
result <- EMP_chipseq_analysis(
    peak_file = "peaks.narrowPeak",
    params = list(
        annotation = list(genome = "hg38"),
        enrichment = list(go = TRUE, kegg = TRUE)
    ),
    output_dir = "./chipseq_results"
)

# 查看注释结果
head(result$peak_annotation)
```

#### 多组学整合

```r
# ChIP-seq + RNA-seq联合分析
result <- EMP_integrate_chipseq_rnaseq(
    chipseq_data = chip_counts,
    rnaseq_data = rna_counts,
    output_dir = "./integration_results"
)

# 查看因子关联
plot_factor_correlation(result$model)
```

### Web界面使用

1. **上传数据**
   - 点击"上传数据"按钮
   - 选择CSV/TSV/Excel格式的数据文件
   - 可选上传metadata文件

2. **选择分析模块**
   - 从9个模块中选择合适的分析类型
   - 查看模块功能介绍

3. **配置参数**
   - 根据数据特点调整质控参数
   - 选择分析方法（如DESeq2/edgeR）
   - 设置可视化选项

4. **运行分析**
   - 点击"开始分析"
   - 实时查看分析进度

5. **查看结果**
   - 交互式图表查看
   - 下载PDF报告
   - 导出数据表格

---

## 🏗️ 架构设计

### 系统架构图

```
┌─────────────────────────────────────────────────────────────────┐
│                        用户层                                    │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐          │
│  │   Web界面    │  │   R控制台    │  │   API调用    │          │
│  └──────┬───────┘  └──────┬───────┘  └──────┬───────┘          │
└─────────┼─────────────────┼─────────────────┼──────────────────┘
          │                 │                 │
          └─────────────────┴─────────────────┘
                            │
┌───────────────────────────▼───────────────────────────────────┐
│                      接口层 (API)                              │
│              Flask REST API / R Functions                      │
└───────────────────────────┬───────────────────────────────────┘
                            │
┌───────────────────────────▼───────────────────────────────────┐
│                      分析层                                    │
│  ┌────────────┐ ┌────────────┐ ┌────────────┐                │
│  │  Seurat    │ │ ChIPseeker │ │   MOFA2    │                │
│  │  (单细胞)  │ │(ChIP-seq)  │ │(多组学)    │                │
│  └────────────┘ └────────────┘ └────────────┘                │
│  ┌────────────┐ ┌────────────┐ ┌────────────┐                │
│  │  DESeq2    │ │  limma     │ │ clusterPro-│                │
│  │  (差异)    │ │ (蛋白质)   │ │filer(富集) │                │
│  └────────────┘ └────────────┘ └────────────┘                │
└───────────────────────────────────────────────────────────────┘
```

### 核心组件

| 组件 | 技术栈 | 说明 |
|------|--------|------|
| **前端** | React + Ant Design | 用户交互界面 |
| **后端** | Flask + Python | API服务和任务调度 |
| **R包** | R + Bioconductor | 核心分析算法 |
| **数据库** | SQLite/PostgreSQL | 任务和结果存储 |
| **容器** | Docker + Docker Compose | 部署和分发 |

---

## 🔌 API文档

### REST API端点

#### 健康检查
```http
GET /api/health

Response:
{
    "status": "ok",
    "version": "3.0.0",
    "service": "EasyMultiProfiler"
}
```

#### 获取分析模块列表
```http
GET /api/modules

Response:
{
    "success": true,
    "modules": [
        {
            "id": "rnaseq",
            "name": "RNA-seq分析",
            "icon": "📊",
            "features": ["差异表达", "火山图", "热图"]
        }
    ]
}
```

#### 上传数据
```http
POST /api/upload
Content-Type: multipart/form-data

file: <数据文件>

Response:
{
    "success": true,
    "file_id": "uuid",
    "samples": 100,
    "features": 20000
}
```

#### 提交分析任务
```http
POST /api/analyze
Content-Type: application/json

{
    "file_id": "uuid",
    "module": "rnaseq",
    "params": {
        "de": {"method": "DESeq2"},
        "enrichment": {"database": "GO_KEGG"}
    }
}

Response:
{
    "success": true,
    "task_id": "task-uuid",
    "status": "queued"
}
```

#### 查询任务状态
```http
GET /api/status/{task_id}

Response:
{
    "task_id": "task-uuid",
    "status": "running",
    "progress": 75,
    "message": "正在运行差异分析..."
}
```

#### 获取结果
```http
GET /api/results/{task_id}

Response:
{
    "success": true,
    "plots": [...],
    "tables": [...],
    "report_url": "/results/task-uuid/report.pdf"
}
```

---

## 🔗 相关链接

### 项目仓库

| 项目 | 链接 | 说明 |
|------|------|------|
| **EasyMultiProfiler-V3** | https://github.com/xielab2017/EasyMultiProfiler-V3 | 本仓库（统一版） |
| **EasyMultiProfiler** (R包) | https://github.com/xielab2017/EasyMultiProfiler | R包源代码 |
| **EasyMultiProfiler-V2** (Web) | https://github.com/xielab2017/EasyMultiProfiler-V2 | Web版源代码 |

### 文档资源

- 📖 [完整文档](https://easymultiprofiler.xielab.net/docs)
- 🎓 [使用教程](https://easymultiprofiler.xielab.net/tutorial)
- 📊 [示例 gallery](https://easymultiprofiler.xielab.net/gallery)
- 💬 [讨论区](https://github.com/xielab2017/EasyMultiProfiler-V3/discussions)

### 依赖项目

- [Seurat](https://satijalab.org/seurat/) - 单细胞分析
- [ChIPseeker](https://guangchuangyu.github.io/software/ChIPseeker/) - ChIP-seq分析
- [MOFA2](https://biofam.github.io/MOFA2/) - 多组学整合
- [Flask](https://flask.palletsprojects.com/) - Web后端
- [React](https://reactjs.org/) - Web前端

---

## 📄 引用

如果您使用了 EasyMultiProfiler，请引用：

```
Li X, et al. EasyMultiProfiler: An Efficient and Scalable Multi-Omics Analysis Platform.
Science China Life Sciences (2025), DOI: 10.1007/s11427-025-3035-0
```

---

## 📞 联系我们

- 📧 邮箱: contact@xielab.net
- 🏠 主页: https://xielab.net
- 🐛 问题反馈: https://github.com/xielab2017/EasyMultiProfiler-V3/issues
- 💡 功能建议: https://github.com/xielab2017/EasyMultiProfiler-V3/discussions

---

## 📜 许可证

本项目采用 MIT 许可证 - 详见 [LICENSE](LICENSE) 文件

---

<div align="center">

**Made with ❤️ by XieLab**

</div>
