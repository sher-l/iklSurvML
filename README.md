# iklSurvML

<!-- badges: start -->
<!-- badges: end -->

Optimized Survival Machine Learning with 117 Algorithm Combinations.

**English** | [中文](#中文文档)

## Features

- **117 ML Combinations**: 10 algorithms (RSF, Enet, StepCox, CoxBoost, plsRcox, SuperPC, GBM, survival-SVM, Ridge, Lasso) with feature selection combinations
- **2.69x Faster**: Feature selection caching optimization
- **100% Compatible**: Results match original Mime package
- **Modular Design**: Clean, maintainable code structure

## Installation

```r
# Install dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_deps <- c('GSEABase', 'GSVA', 'mixOmics', 'sva', 'ComplexHeatmap')
for (pkg in bioc_deps) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, update = FALSE)
}

# Install CoxBoost
if (!requireNamespace("CoxBoost", quietly = TRUE))
  devtools::install_github("binderh/CoxBoost")

# Install iklSurvML
devtools::install_github("sher-l/iklSurvML")
```

## Quick Start

```r
library(iklSurvML)

# Load example data
load("./External data/Example.cohort.Rdata")
load("./External data/genelist.Rdata")

# Run all 117 combinations
result <- ML.Dev.Prog.Sig(
  train_data = list_train_vali_Data$Dataset1,
  list_train_vali_Data = list_train_vali_Data,
  candidate_genes = genelist,
  mode = "all",
  seed = 5201314
)

# Or use optimized version (2.69x faster)
result_fast <- ML.Dev.Prog.Sig.Fast(
  train_data = list_train_vali_Data$Dataset1,
  list_train_vali_Data = list_train_vali_Data,
  candidate_genes = genelist,
  mode = "all",
  seed = 5201314
)
```

## 117 Combinations

| Category | Count | Algorithms |
|----------|-------|------------|
| Single models | 20 | RSF, Enet(9), StepCox(3), CoxBoost, plsRcox, SuperPC, GBM, survival-SVM, Ridge, Lasso |
| RSF + X | 19 | RSF + (CoxBoost, Enet(9), GBM, Lasso, plsRcox, Ridge, SuperPC, survival-SVM, StepCox(3)) |
| StepCox + X | 51 | 3 directions × 17 combinations |
| CoxBoost + X | 19 | CoxBoost + (Enet(9), GBM, Lasso, plsRcox, Ridge, StepCox(3), SuperPC, survival-SVM) |
| Lasso + X | 9 | Lasso + (CoxBoost, GBM, plsRcox, RSF, StepCox(3), SuperPC, survival-SVM) |
| **Total** | **117** | |

## API

### Main Functions

- `ML.Dev.Prog.Sig()` - Original implementation (backward compatible)
- `ML.Dev.Prog.Sig.Fast()` - Optimized version with caching

### Visualization

- `cindex_dis_all()` - Plot C-index of all models
- `cindex_dis_select()` - Plot C-index of specific model
- `rs_sur()` - Survival curve by risk score
- `auc_dis_all()` - Plot AUC of all models
- `roc_vis()` - ROC curve visualization

---

# 中文文档

基于 117 种机器学习算法组合的优化生存分析包。

## 特点

- **117 种 ML 组合**：10 种算法（RSF, Enet, StepCox, CoxBoost, plsRcox, SuperPC, GBM, survival-SVM, Ridge, Lasso）的特征选择组合
- **2.69 倍加速**：特征选择缓存优化
- **100% 兼容**：结果与原始 Mime 包完全一致
- **模块化设计**：清晰、可维护的代码结构

## 安装

```r
# 安装依赖
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_deps <- c('GSEABase', 'GSVA', 'mixOmics', 'sva', 'ComplexHeatmap')
for (pkg in bioc_deps) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, update = FALSE)
}

# 安装 CoxBoost
if (!requireNamespace("CoxBoost", quietly = TRUE))
  devtools::install_github("binderh/CoxBoost")

# 安装 iklSurvML
devtools::install_github("sher-l/iklSurvML")
```

## 快速开始

```r
library(iklSurvML)

# 加载示例数据
load("./External data/Example.cohort.Rdata")
load("./External data/genelist.Rdata")

# 运行全部 117 种组合
result <- ML.Dev.Prog.Sig(
  train_data = list_train_vali_Data$Dataset1,
  list_train_vali_Data = list_train_vali_Data,
  candidate_genes = genelist,
  mode = "all",
  seed = 5201314
)

# 或使用优化版本（2.69 倍加速）
result_fast <- ML.Dev.Prog.Sig.Fast(
  train_data = list_train_vali_Data$Dataset1,
  list_train_vali_Data = list_train_vali_Data,
  candidate_genes = genelist,
  mode = "all",
  seed = 5201314
)
```

## 117 种组合详解

| 类别 | 数量 | 算法 |
|------|------|------|
| 单模型 | 20 | RSF, Enet(9), StepCox(3), CoxBoost, plsRcox, SuperPC, GBM, survival-SVM, Ridge, Lasso |
| RSF + X | 19 | RSF + (CoxBoost, Enet(9), GBM, Lasso, plsRcox, Ridge, SuperPC, survival-SVM, StepCox(3)) |
| StepCox + X | 51 | 3 个方向 × 17 种组合 |
| CoxBoost + X | 19 | CoxBoost + (Enet(9), GBM, Lasso, plsRcox, Ridge, StepCox(3), SuperPC, survival-SVM) |
| Lasso + X | 9 | Lasso + (CoxBoost, GBM, plsRcox, RSF, StepCox(3), SuperPC, survival-SVM) |
| **总计** | **117** | |

## 主要函数

### 模型构建

- `ML.Dev.Prog.Sig()` - 原始实现（向后兼容）
- `ML.Dev.Prog.Sig.Fast()` - 带缓存的优化版本

### 可视化

- `cindex_dis_all()` - 绘制所有模型的 C-index
- `cindex_dis_select()` - 绘制特定模型的 C-index
- `rs_sur()` - 按风险评分绘制生存曲线
- `auc_dis_all()` - 绘制所有模型的 AUC
- `roc_vis()` - ROC 曲线可视化

## 引用

如果您使用本包，请引用：

- Liu H, Zhang W, Zhang Y, et al. Mime: A flexible machine-learning framework to construct and visualize models for clinical characteristics prediction and feature selection. Comput Struct Biotechnol J. 2024.

## 联系方式

如有问题请在 Issues 中提出。

## License

MIT
