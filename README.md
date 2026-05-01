# iklSurvML

> High-Performance Survival Machine Learning | Fixed 117 Algorithm Combinations | 100% Reproducible

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
| Algorithm selection | Run candidate combinations, then evaluate on a pre-specified validation strategy |
| Slow performance | **~8x faster** (optimized code) + **12-core parallel** |
| Non-reproducible | Fixed random seed, 100% reproducible results |

## Key Features

### 1. Comprehensive Algorithm Coverage

```
10 Base Algorithms → Fixed 117 All-Mode Combinations
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

**Reproducibility Guarantee:**

All 10 algorithms are 100% reproducible with proper seed handling:
- Same `seed` → identical C-index every time
- All random operations (CV splits, bootstrap, etc.) are seeded
- RSF: `var.select()` seeded for variable selection stability
- CoxBoost: `parallel=FALSE` + seeded for deterministic boosting
- SuperPC: `superpc.cv()` seeded for cross-validation stability

**Test Evidence**: See `test/` directory for validation scripts

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
| `double` | Two distinct algorithms | Feature selection + modeling; self-combinations such as `Lasso + Lasso` or `StepCox + StepCox` are not supported |
| `all` | Fixed 117 Mime-style combinations | Exploratory analysis |

`mode = "all"` now requires the complete fixed 117-model Mime-style grid.
`model_grid` is retained only for compatibility and must be `"117"` if set.
If a weak dataset causes a first-stage selector or backend model to fail, the
run stops with a clear message. Set `allow_partial = TRUE` only when you
intentionally want to compare the successfully fitted subset.

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

### Optional classification engines

Survival modeling installs the required survival backends. Some classification
helpers in `ML.Dev.Pred.Category.Sig()` use optional engines and are enabled
only when their packages are installed:

| Method | Extra package(s) |
|--------|------------------|
| `nb` | `klaR` |
| `adaboost` | `fastAdaboost` |
| `cancerclass` | `cancerclass`, `Biobase`, `pROC` |

If you request one of these methods without its extra packages, iklSurvML fails
early with a method-specific dependency message. Leaving `methods = NULL`
automatically skips unavailable optional engines.

## Quick Start

### Basic Usage

```r
library(iklSurvML)

# Run the fixed 117 all-mode combinations
result <- ML.Dev.Prog.Sig(
  train_data = train_data,
  list_train_vali_Data = list_train_vali_Data,
  candidate_genes = gene_list,
  mode = "all",
  nodesize = 5,
  seed = 12345
)


# Inspect validation performance after pre-specifying the model-selection cohort.
# Avoid choosing the final model from training rows or from repeated external-test looks.
validation_cindex <- subset(result$Cindex.res, ID != "train")
validation_cindex[order(validation_cindex$Cindex, decreasing = TRUE), ]
```

### Fast Version with Optional Parallel Execution

```r
# Sequential all-mode is the fork-safe default
result_fast <- ML.Dev.Prog.Sig.Fast(
  train_data = train_data,
  list_train_vali_Data = list_train_vali_Data,
  candidate_genes = gene_list,
  mode = "all",
  nodesize = 5,
  seed = 12345,
  use_parallel = FALSE,       # Default: FALSE (fork-safe)
  cores_for_parallel = 12     # Number of cores (default: 12)
)

# Opt in to parallel execution on Linux/macOS when desired.
# GBM tasks stay in the parent process to avoid fork-related crashes.
result_parallel <- ML.Dev.Prog.Sig.Fast(
  ...,
  use_parallel = TRUE
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
| StepCox + X | 51 | Uses both/backward/forward as first-stage selectors |
| CoxBoost + X | 18 | Keeps the Mime-style omission of CoxBoost + RSF; RSF + CoxBoost remains |
| Lasso + X | 9 | Sparse feature selection; omits Enet/Ridge after Lasso selection |
| **Total** | **117** | Fixed all-mode grid |

**Recommendation:** Pre-specify the validation/cohort split used for model comparison before running `mode="all"`; do not repeatedly pick the highest C-index on the final external test set.

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
| `model_grid` | character | "117" | Compatibility parameter; only "117" is supported for all-mode |
| `use_parallel` | logical | FALSE | Enable forked parallel execution for non-GBM tasks |
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
# Option 1: Use Fast version with safe sequential all-mode (default)
result <- ML.Dev.Prog.Sig.Fast(
  ...,
  use_parallel = FALSE,
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
- GBM model fits are run in the parent process for fork safety

## Changelog

### v1.3.0
- 🐛 Fix 24 bugs across two comprehensive audit rounds
  - **CRITICAL**: predict type `"link"` for Cox models, C-index sign correction, factor→numeric removal, closure variable capture
  - **HIGH**: StepCox direction parameter, CoxBoost combination coverage, tryCatch error handling
  - **MEDIUM**: input validation, SuperPC null checks, parallel Windows warning, etc.
  - **LOW**: message formatting, parameter documentation, default values
- 🔧 Keep a fixed `model_grid = "117"` all-mode grid aligned with Mime-style 117 combinations
- 🔧 CoxBoost combinations now use non-zero CoxBoost coefficients for first-stage feature selection
- 🔧 C-index now evaluates risk-score ordering directly instead of refitting a Cox model on each validation cohort
- 🧹 Remove the dead duplicate `R/IMPRES` source tree; keep only the packaged runtime copies under `inst/extdata` and remove their hard-coded legacy RNG seed
- 🔢 Total all-mode algorithm combinations: fixed 117
- 🧹 Remove duplicate validation blocks, clean up code structure

### v1.2.0
- 🐛 Fix seed handling for RSF `var.select()`, survivalsvm, CoxBoost, SuperPC
- ✅ All 10 algorithms now 100% reproducible with same seed
- 🏷️ StepCox model names now include direction (e.g., `StepCox[both]`)

### v1.1.0
- ✨ Add 12-core parallel execution for all-mode combinations
- ✅ 100% consistency with Mime package (10/10 algorithms)

## Citation

Based on the Mime framework:

> Liu H, Zhang W, Zhang Y, et al. Mime: A flexible machine-learning framework to construct and visualize models for clinical characteristics prediction and feature selection. *Comput Struct Biotechnol J*. 2024.

## License

MIT License

---

# 中文文档

## 简介

iklSurvML 是专注于生存分析的机器学习工具包，固定提供 Mime 风格 117 种 all-mode 算法组合，帮助研究者快速构建和筛选最优预测模型。

## 核心特性

| 特性 | 说明 |
|------|------|
| 全面覆盖 | 集成 10 种主流生存分析算法 |
| 高效运行 | **~8 倍加速** (代码优化) + **12 核并行** |
| 结果可靠 | 100% 可复现，10/10 算法与 Mime 完全一致 |
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
| RSF | 0.000000 | ✅ |
| SuperPC | 0.000000 | ✅ |
| **总计** | **10/10 (100%)** | ✅ |

**可复现性保证：**

所有 10 种算法都通过正确的种子处理实现 100% 可复现：
- 相同 `seed` → 每次运行 C-index 完全一致
- 所有随机操作（CV 划分、bootstrap 等）都设置了种子
- RSF：`var.select()` 设置种子确保变量选择稳定性
- CoxBoost：禁用并行 + 设置种子确保确定性提升
- SuperPC：`superpc.cv()` 设置种子确保交叉验证稳定性

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

# 安全顺序运行固定 117 种组合
result <- ML.Dev.Prog.Sig.Fast(
  train_data = train,
  list_train_vali_Data = list(train = train, val = validation),
  candidate_genes = genes,
  mode = "all",
  nodesize = 5,
  seed = 12345,
  use_parallel = FALSE,       # 默认关闭并行，更稳妥
  cores_for_parallel = 12     # CPU 核心数
)


# 在预先指定的验证队列中查看表现；不要用训练集行或反复查看外部测试集来选最终模型
validation_cindex <- subset(result$Cindex.res, ID != "train")
validation_cindex[order(validation_cindex$Cindex, decreasing = TRUE), ]
```

## 使用建议

1. **数据准备**：样本量 ≥100，基因数 ≥50
2. **首选模式**：先用 `mode="all"` 跑完全部组合；需要加速时再显式开启并行
3. **模型选择**：预先指定内部验证/调参策略；避免在训练集或最终外部测试集上反复挑选最高 C-index
4. **结果验证**：锁定模型后，在多个独立验证集中只评估一次稳定性

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
# 使用 Fast 版本；默认顺序执行更稳妥
result <- ML.Dev.Prog.Sig.Fast(..., use_parallel = FALSE)
```

**结果不一致？** 确保 `seed`、`nodesize` 参数相同

**报错？** 记得设置 `nodesize = 5`

**并行不生效？** 并行使用 Linux/macOS fork，Windows 会自动降级为顺序执行；GBM 会留在主进程执行以避免 fork 崩溃

## 更新日志

### v1.3.0
- 🐛 两轮全面审计共修复 24 个 Bug
  - **严重**: Cox 模型 predict type 修正为 `"link"`、C-index 符号修正、factor→numeric 移除、闭包变量捕获
  - **高危**: StepCox 方向参数、CoxBoost 组合覆盖、tryCatch 错误处理
  - **中等**: 输入验证、SuperPC 空值检查、并行 Windows 警告 等
  - **低危**: 消息格式化、参数文档、默认值
- 🔧 固定使用 `model_grid = "117"` 的 Mime 风格 all-mode 组合
- 🔧 CoxBoost 组合作为一阶段时改为使用非零 CoxBoost 系数筛选特征
- 🔧 C-index 改为直接评估风险分数排序，不再在每个验证集重新拟合 Cox 模型
- 🧹 删除无效重复的 `R/IMPRES` 源码树；仅保留 `inst/extdata` 中当前运行会加载的副本，并移除其中写死的历史 RNG 种子
- 🔢 all-mode 算法组合总数: 固定 117
- 🧹 清除重复验证代码块，优化代码结构

### v1.2.0
- 🐛 修复 RSF `var.select()`、survivalsvm、CoxBoost、SuperPC 的种子处理
- ✅ 所有 10 种算法使用相同种子 100% 可复现
- 🏷️ StepCox 模型名称现在包含方向 (如 `StepCox[both]`)

### v1.1.0
- ✨ 新增 12 核并行执行
- ✅ 100% 一致性验证通过 (10/10 算法)

## 获取帮助

- Issues: https://github.com/sher-l/iklSurvML/issues
- 示例数据: `External_data/` 目录

---

<p align="center">
Made with ❤️ for survival analysis researchers
</p>
