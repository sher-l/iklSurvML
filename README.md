# iklSurvML

> High-Performance Survival Machine Learning | 117 Algorithm Combinations | 100% Reproducible

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue.svg)](https://www.r-project.org/)

**English** | [中文](#中文文档)

---

## Why iklSurvML?

When building survival prediction models, do you face these challenges?

- ❌ Unsure which algorithm works best for your data?
- ❌ Running all combinations takes hours or days?
- ❌ Results are hard to reproduce across different packages?

**iklSurvML solves these problems:**

| Challenge | Solution |
|-----------|----------|
| Algorithm selection | Run 117 combinations at once, auto-select the best |
| Slow performance | Smart caching, 2.69x faster |
| Non-reproducible | Fixed random seed, 100% reproducible results |

## Key Features

### 1. Comprehensive Algorithm Coverage

```
10 Base Algorithms → 117 Combinations
├── Regularization: Lasso, Ridge, Elastic Net (9 α values)
├── Ensemble: Random Survival Forest, GBM, CoxBoost
├── Classical: Stepwise Cox (3 directions)
└── Others: plsRcox, SuperPC, survival-SVM
```

### 2. Validated Consistency

Tested against the original Mime package:

```
Test Coverage: 4 datasets × 2 algorithms = 8 tests
Consistency:   100% (C-index difference = 0)
```

### 3. Flexible Modes

| Mode | Description | Use Case |
|------|-------------|----------|
| `single` | Single algorithm | Algorithm already determined |
| `double` | Two-algorithm combo | Feature selection + modeling |
| `all` | All 117 combinations | Exploratory analysis |

## Installation

```r
# 1. Install Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_pkgs <- c('GSEABase', 'GSVA', 'mixOmics', 'sva', 'ComplexHeatmap')
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, update = FALSE)
}

# 2. Install CoxBoost from GitHub
if (!requireNamespace("CoxBoost", quietly = TRUE))
  devtools::install_github("binderh/CoxBoost")

# 3. Install iklSurvML
devtools::install_github("sher-l/iklSurvML")
```

## Quick Start

### Basic Usage

```r
library(iklSurvML)

# Run all 117 combinations
result <- ML.Dev.Prog.Sig(
  train_data = train_data,
  list_train_vali_Data = list_train_vali_Data,
  candidate_genes = gene_list,
  mode = "all",
  nodesize = 5,
  seed = 12345
)

# View best model
result$Cindex.res[which.max(result$Cindex.res$Cindex), ]
```

### Fast Version (2.69x Speedup)

```r
result_fast <- ML.Dev.Prog.Sig.Fast(
  train_data = train_data,
  list_train_vali_Data = list_train_vali_Data,
  candidate_genes = gene_list,
  mode = "all",
  nodesize = 5,
  seed = 12345
)
```

### Single Algorithm Mode

```r
result_lasso <- ML.Dev.Prog.Sig(
  train_data = train_data,
  list_train_vali_Data = list_train_vali_Data,
  candidate_genes = gene_list,
  mode = "single",
  single_ml = "Lasso",
  nodesize = 5,
  seed = 12345
)
```

## Data Format

```r
# Required columns:
#   ID      OS.time   OS    Gene1    Gene2    Gene3    ...
#   Sample1  1000      1     10.5     8.2      12.1    ...
#   Sample2   500      0     9.8      7.5      11.3    ...

# Column descriptions:
# - ID: Unique sample identifier
# - OS.time: Survival time (days/months)
# - OS: Survival status (1=dead, 0=censored)
# - Other columns: Gene expression (log2 transformed recommended)
```

## Algorithm Combinations

| Type | Count | Description |
|------|-------|-------------|
| Single models | 20 | Baseline comparisons |
| RSF + X | 19 | Random forest feature selection |
| StepCox + X | 51 | Classical statistics + ML |
| CoxBoost + X | 19 | Boosting feature selection |
| Lasso + X | 9 | Sparse feature selection |
| **Total** | **117** | |

**Recommendation:** Run `mode="all"` first, then select the best model by C-index.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `train_data` | data.frame | required | Training dataset |
| `list_train_vali_Data` | list | required | All datasets (incl. training) |
| `candidate_genes` | character | NULL | Candidate gene list |
| `mode` | character | "all" | all/single/double |
| `single_ml` | character | NULL | Algorithm for single mode |
| `unicox.filter.for.candi` | logical | TRUE | Filter genes by univariate Cox |
| `unicox_p_cutoff` | numeric | 0.05 | p-value threshold |
| `nodesize` | numeric | NULL | RSF node size (set to 5) |
| `seed` | numeric | NULL | Random seed |

## Visualization

```r
# C-index distribution
cindex_dis_all(result, validate_set = c("Val1", "Val2"))

# Specific model C-index
cindex_dis_select(result, model = "StepCox[forward] + plsRcox")

# Survival curve
rs_sur(result, model_name = "Lasso", dataset = "training", cutoff = 0.5)

# ROC curve
roc_vis(auc_result, model_name = "Enet[α=0.5]", year = 1)
```

## FAQ

### Q: Running too slow?

```r
# Option 1: Use Fast version
result <- ML.Dev.Prog.Sig.Fast(...)

# Option 2: Pre-filter genes
# Set unicox.filter.for.candi = TRUE

# Option 3: Single algorithm mode
result <- ML.Dev.Prog.Sig(..., mode = "single", single_ml = "Lasso")
```

### Q: Results differ from Mime?

Ensure identical parameters:
- Same `seed`
- Same `nodesize`
- Same `unicox.filter.for.candi`

### Q: Error "argument is of length zero"?

Set `nodesize` parameter:
```r
# Wrong
result <- ML.Dev.Prog.Sig(..., seed = 12345)

# Correct
result <- ML.Dev.Prog.Sig(..., nodesize = 5, seed = 12345)
```

## Citation

Based on the Mime framework:

> Liu H, Zhang W, Zhang Y, et al. Mime: A flexible machine-learning framework to construct and visualize models for clinical characteristics prediction and feature selection. *Comput Struct Biotechnol J*. 2024.

## License

MIT License

---

# 中文文档

## 简介

iklSurvML 是专注于生存分析的机器学习工具包，提供 117 种算法组合，帮助研究者快速构建和筛选最优预测模型。

## 核心特性

| 特性 | 说明 |
|------|------|
| 全面覆盖 | 集成 10 种主流生存分析算法 |
| 高效运行 | 智能缓存加速 2.69 倍 |
| 结果可靠 | 100% 可复现，与 Mime 完全兼容 |
| 易于使用 | 简洁 API，详细文档 |

## 安装

```r
# Bioconductor 依赖
BiocManager::install(c('GSEABase', 'GSVA', 'mixOmics', 'sva', 'ComplexHeatmap'))

# CoxBoost
devtools::install_github("binderh/CoxBoost")

# iklSurvML
devtools::install_github("sher-l/iklSurvML")
```

## 快速开始

```r
library(iklSurvML)

# 运行全部 117 种组合
result <- ML.Dev.Prog.Sig(
  train_data = train,
  list_train_vali_Data = list(train = train, val = validation),
  candidate_genes = genes,
  mode = "all",
  nodesize = 5,
  seed = 12345
)

# 查看结果
cindex_dis_all(result)
```

## 使用建议

1. **数据准备**：样本量 ≥100，基因数 ≥50
2. **首选模式**：先用 `mode="all"` 跑完全部组合
3. **模型选择**：根据 C-index 选择最优模型
4. **结果验证**：在多个独立验证集中确认稳定性

## 支持的算法

```
├── Lasso / Ridge / Elastic Net
├── Random Survival Forest
├── Stepwise Cox (forward/backward/both)
├── CoxBoost
├── plsRcox
├── SuperPC
├── GBM
└── survival-SVM
```

## 常见问题

**运行慢？** 使用 `ML.Dev.Prog.Sig.Fast()` 加速版本

**结果不一致？** 确保 `seed`、`nodesize` 参数相同

**报错？** 记得设置 `nodesize = 5`

## 获取帮助

- Issues: https://github.com/sher-l/iklSurvML/issues
- 示例数据: `External data/` 目录

---

<p align="center">
Made with ❤️ for survival analysis researchers
</p>
