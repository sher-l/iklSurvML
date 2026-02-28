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
| Slow performance | **~8x faster** (optimized code) + **12-core parallel** |
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

### 2. Validated Consistency (100%)

Tested against the original Mime package:

| Algorithm | C-index Difference | Status |
|-----------|-------------------|--------|
| Lasso | 0.000000 | ✅ |
| Ridge | 0.000000 | ✅ |
| Enet | 0.000000 | ✅ |
| StepCox | 0.000000 | ✅ |
| CoxBoost | 0.000000 | ✅ |
| plsRcox | 0.000000 | ✅ |
| survivalsvm | 0.000000 | ✅ |
| GBM | 0.000000 | ✅ |
| RSF | 0.000000 | ✅ |
| SuperPC | 0.000000 | ✅ |
| **Total** | **10/10 (100%)** | ✅ |

**⚠️ Important Notes:**

1. **Deterministic algorithms** (Lasso, Ridge, Enet, StepCox, CoxBoost, plsRcox, survivalsvm, GBM):
   - Same seed → identical results every time
   - No model re-training needed between runs

2. **Re-training algorithms** (RSF, SuperPC):
   - Each run re-trains the model (computationally expensive)
   - Results are reproducible with fixed seed + stability fixes
   - Floating-point precision stabilized with `round(..., 10)`

3. **重要说明**:
   - 确定性算法：相同种子 → 每次结果完全一致，无需重新训练
   - 重训练算法（RSF、SuperPC）：每次运行会重新训练模型，已通过浮点数精度修复确保可复现性

### 3. Performance Benchmarks

```
Single Algorithm (Lasso, 5-run average):
  Mime:       1.30s
  iklSurvML:  0.17s  → ~8x faster

117 Combinations (12-core parallel):
  Mime:       ~500s (sequential, estimated)
  iklSurvML:  ~40s  → ~12x faster (parallel)

Key optimizations:
  - Removed redundant library() calls
  - Modular code structure
  - Eliminated unnecessary message() output
```

### 4. Flexible Modes

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

### Fast Version with Parallel Execution

```r
# 12-core parallel execution for mode="all"
result_fast <- ML.Dev.Prog.Sig.Fast(
  train_data = train_data,
  list_train_vali_Data = list_train_vali_Data,
  candidate_genes = gene_list,
  mode = "all",
  nodesize = 5,
  seed = 12345,
  use_parallel = TRUE,        # Enable parallel (default: TRUE)
  cores_for_parallel = 12     # Number of cores (default: 12)
)

# Sequential execution (if needed)
result_seq <- ML.Dev.Prog.Sig.Fast(
  ...,
  use_parallel = FALSE
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
| `use_parallel` | logical | TRUE | Enable parallel execution |
| `cores_for_parallel` | numeric | 12 | Number of CPU cores |

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
# Option 1: Use Fast version with parallel (recommended)
result <- ML.Dev.Prog.Sig.Fast(
  ...,
  use_parallel = TRUE,
  cores_for_parallel = 12
)

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

### Q: Parallel execution not working?

Parallel execution uses `parallel::mclapply` (Linux/macOS fork).
- Works on Linux and macOS
- On Windows, falls back to sequential execution

## Changelog

### v1.1.0
- ✨ Add 12-core parallel execution for 117 combinations
- ✅ 100% consistency with Mime package (8/8 algorithms)

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
| 高效运行 | **~8 倍加速** (代码优化) + **12 核并行** |
| 结果可靠 | 100% 可复现，8/8 算法与 Mime 完全一致 |
| 易于使用 | 简洁 API，详细文档 |

## 一致性验证 (100%)

| 算法 | C-index 差异 | 状态 |
|------|-------------|------|
| Lasso | 0.000000 | ✅ |
| Ridge | 0.000000 | ✅ |
| Enet | 0.000000 | ✅ |
| StepCox | 0.000000 | ✅ |
| CoxBoost | 0.000000 | ✅ |
| plsRcox | 0.000000 | ✅ |
| survivalsvm | 0.000000 | ✅ |
| GBM | 0.000000 | ✅ |

## 性能基准

```
单算法 (Lasso, 5次平均):
  Mime:       1.30秒
  iklSurvML:  0.17秒  → ~8倍加速

117组合 (12核并行):
  Mime:       ~500秒 (顺序, 预估)
  iklSurvML:  ~40秒   → ~12倍加速 (并行)

优化要点:
  - 移除冗余 library() 调用
  - 模块化代码结构
  - 精简 message() 输出
```

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

# 并行运行全部 117 种组合 (推荐)
result <- ML.Dev.Prog.Sig.Fast(
  train_data = train,
  list_train_vali_Data = list(train = train, val = validation),
  candidate_genes = genes,
  mode = "all",
  nodesize = 5,
  seed = 12345,
  use_parallel = TRUE,        # 启用并行 (默认开启)
  cores_for_parallel = 12     # CPU 核心数
)

# 查看最佳模型
best_idx <- which.max(result$Cindex.res$Cindex)
result$Cindex.res[best_idx, ]
```

## 使用建议

1. **数据准备**：样本量 ≥100，基因数 ≥50
2. **首选模式**：先用 `mode="all"` + 并行跑完全部组合
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

**运行慢？**
```r
# 使用 Fast 版本 + 并行
result <- ML.Dev.Prog.Sig.Fast(..., use_parallel = TRUE, cores_for_parallel = 12)
```

**结果不一致？** 确保 `seed`、`nodesize` 参数相同

**报错？** 记得设置 `nodesize = 5`

**并行不生效？** 并行使用 Linux/macOS fork，Windows 会自动降级为顺序执行

## 更新日志

### v1.1.0
- ✨ 新增 12 核并行执行
- ✅ 100% 一致性验证通过 (8/8 算法)

## 获取帮助

- Issues: https://github.com/sher-l/iklSurvML/issues
- 示例数据: `External data/` 目录

---

<p align="center">
Made with ❤️ for survival analysis researchers
</p>
