# iklSurvML

> 高性能生存分析机器学习框架 | 117种算法组合 | 100%结果可复现

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue.svg)](https://www.r-project.org/)

---

**English** | [中文](#中文文档)

## 为什么选择 iklSurvML？

构建生存预测模型时，你是否面临这些问题：

- ❌ 不知道哪种算法最适合你的数据？
- ❌ 跑完全部算法组合需要数小时甚至数天？
- ❌ 不同包的结果难以复现？

**iklSurvML 解决这些问题：**

| 痛点 | 解决方案 |
|------|----------|
| 算法选择困难 | 一键运行 117 种组合，自动筛选最优模型 |
| 运行速度慢 | 智能缓存机制，加速 2.69 倍 |
| 结果不可复现 | 固定随机种子，结果 100% 可复现 |

## 核心优势

### 1. 全面覆盖主流算法

```
10 种基础算法 → 117 种组合
├── 正则化方法: Lasso, Ridge, Elastic Net (9种α值)
├── 集成学习: Random Survival Forest, GBM, CoxBoost
├── 经典统计: Stepwise Cox (3种方向)
└── 其他: plsRcox, SuperPC, survival-SVM
```

### 2. 经过验证的一致性

与原始 Mime 包对比测试结果：

```
测试数据: 4组数据 × 2种算法 = 8次测试
一致率:   100% (C-index差异 = 0)
```

### 3. 灵活的使用模式

| 模式 | 说明 | 适用场景 |
|------|------|----------|
| `single` | 单一算法 | 已确定算法，快速验证 |
| `double` | 两个算法组合 | 特征选择 + 模型构建 |
| `all` | 全部117种组合 | 探索性分析，筛选最优 |

## 安装

```r
# 1. 安装 Bioconductor 依赖
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_pkgs <- c('GSEABase', 'GSVA', 'mixOmics', 'sva', 'ComplexHeatmap')
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, update = FALSE)
}

# 2. 安装 CoxBoost (GitHub)
if (!requireNamespace("CoxBoost", quietly = TRUE))
  devtools::install_github("binderh/CoxBoost")

# 3. 安装 iklSurvML
devtools::install_github("sher-l/iklSurvML")
```

## 快速开始

### 基础用法

```r
library(iklSurvML)

# 准备数据
# train_data: 训练集，包含 ID, OS.time, OS 和基因表达列
# list_train_vali_Data: 训练集和验证集的列表

# 运行全部 117 种组合
result <- ML.Dev.Prog.Sig(
  train_data = train_data,
  list_train_vali_Data = list_train_vali_Data,
  candidate_genes = gene_list,
  mode = "all",
  nodesize = 5,
  seed = 12345  # 固定种子确保可复现
)

# 查看最优模型
result$Cindex.res[which.max(result$Cindex.res$Cindex), ]
```

### 使用加速版本

```r
# 大数据集推荐使用 Fast 版本
result_fast <- ML.Dev.Prog.Sig.Fast(
  train_data = train_data,
  list_train_vali_Data = list_train_vali_Data,
  candidate_genes = gene_list,
  mode = "all",
  nodesize = 5,
  seed = 12345
)
```

### 单一算法模式

```r
# 只运行 Lasso
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

## 数据格式要求

```r
# 训练数据格式
#   ID      OS.time   OS    Gene1    Gene2    Gene3    ...
#   Sample1  1000      1     10.5     8.2      12.1    ...
#   Sample2   500      0     9.8      7.5      11.3    ...

# 列说明:
# - ID: 样本唯一标识符
# - OS.time: 生存时间（天/月）
# - OS: 生存状态（1=死亡，0=删失）
# - 其他列: 基因表达量（建议 log2 转换）
```

## 可视化函数

```r
# C-index 分布图
cindex_dis_all(result, validate_set = c("Validation1", "Validation2"))

# 特定模型的 C-index
cindex_dis_select(result, model = "StepCox[forward] + plsRcox")

# 生存曲线
rs_sur(result, model_name = "Lasso", dataset = "training", cutoff = 0.5)

# ROC 曲线
roc_vis(auc_result, model_name = "Enet[α=0.5]", year = 1)
```

## 算法组合速查表

| 组合类型 | 数量 | 典型应用 |
|----------|------|----------|
| 单模型 | 20 | 基线对比 |
| RSF + X | 19 | 结合随机森林的特征选择 |
| StepCox + X | 51 | 经典统计方法 + ML |
| CoxBoost + X | 19 | Boosting 特征选择 |
| Lasso + X | 9 | 稀疏特征选择 |

**推荐策略：** 先用 `mode="all"` 跑完全部组合，根据 C-index 选择最优模型。

## 参数说明

### ML.Dev.Prog.Sig()

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `train_data` | data.frame | 必需 | 训练集 |
| `list_train_vali_Data` | list | 必需 | 所有数据集（含训练集） |
| `candidate_genes` | character | NULL | 候选基因列表 |
| `mode` | character | "all" | 运行模式：all/single/double |
| `single_ml` | character | NULL | single 模式指定的算法 |
| `unicox.filter.for.candi` | logical | TRUE | 是否用单因素 Cox 筛选基因 |
| `unicox_p_cutoff` | numeric | 0.05 | 单因素 Cox p 值阈值 |
| `nodesize` | numeric | NULL | RSF 节点大小（建议设为 5） |
| `seed` | numeric | NULL | 随机种子 |

## 常见问题

### Q: 运行很慢怎么办？

```r
# 方案1: 使用 Fast 版本
result <- ML.Dev.Prog.Sig.Fast(...)

# 方案2: 先筛选基因
# 使用 unicox.filter.for.candi = TRUE 减少候选基因数量

# 方案3: 只运行特定算法
result <- ML.Dev.Prog.Sig(..., mode = "single", single_ml = "Lasso")
```

### Q: 结果与 Mime 不一致？

确保使用相同的参数：
- `seed` 必须相同
- `nodesize` 必须相同
- `unicox.filter.for.candi` 必须相同

### Q: 报错 "argument is of length zero"？

通常是因为 `nodesize` 参数未设置：
```r
# 错误写法
result <- ML.Dev.Prog.Sig(..., seed = 12345)

# 正确写法
result <- ML.Dev.Prog.Sig(..., nodesize = 5, seed = 12345)
```

## 引用

本包基于 Mime 框架开发，使用时请引用：

> Liu H, Zhang W, Zhang Y, et al. Mime: A flexible machine-learning framework to construct and visualize models for clinical characteristics prediction and feature selection. *Comput Struct Biotechnol J*. 2024.

## 许可证

MIT License

---

# 中文文档

## 简介

iklSurvML 是一个专注于生存分析的机器学习工具包，提供 117 种算法组合，帮助研究者快速构建和筛选最优预测模型。

## 主要特性

- **一站式解决方案**：集成 10 种主流生存分析算法
- **高效运行**：智能缓存加速，支持大规模数据
- **结果可靠**：100% 可复现，与原 Mime 包完全兼容
- **易于使用**：简洁的 API 设计，详细文档

## 安装

```r
# Bioconductor 依赖
BiocManager::install(c('GSEABase', 'GSVA', 'mixOmics', 'sva', 'ComplexHeatmap'))

# CoxBoost
devtools::install_github("binderh/CoxBoost")

# iklSurvML
devtools::install_github("sher-l/iklSurvML")
```

## 30秒上手

```r
library(iklSurvML)

# 加载数据
load("your_data.Rdata")

# 运行分析
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

## 获取帮助

- 提交 Issue: https://github.com/sher-l/iklSurvML/issues
- 查看示例: `External data/` 目录

---

<p align="center">
  Made with ❤️ for survival analysis researchers
</p>
