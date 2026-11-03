---
title: 习题 3.38：由平方和估计误差与组间变异
short_title: 3.38 由平方和估计误差与组间变异
description: 第9版习题3.38的完整中文题目、数据与答案解析
---

(exercise-3-38)=
# 习题 3.38：由平方和估计误差与组间变异

[习题与答案](../../index.md) / [第 3 章　单因子试验与方差分析](../chapter03.md)

第 9 版 · 第 3 章 · 原书章末习题见印刷页 125～134。

## 题目

4 个水平，每水平 3 次重复，总平方和 330.56，处理平方和 250.65。

**（a）** 估计误差方差；

**（b）** 处理解释多少响应变异？

## 答案与解析

### （a）误差方差估计

**模型与模型假设（固定效应）**：$y_{ij}=\mu_i+\epsilon_{ij}=\mu+\tau_i+\epsilon_{ij}$，$i=1,\ldots,a$，$j=1,\ldots,n_i$；$\epsilon_{ij}\overset{\mathrm{iid}}\sim N(0,\sigma^2)$。误差独立、正态且同方差，处理分配及抽样符合题意；不平衡设计使用实际 $n_i$。取 $\sum_i n_i\tau_i=0$ 识别参数；平衡时等价于 $\sum_i\tau_i=0$。

**原假设和备择假设**：$H_0:\mu_1=\cdots=\mu_a$；$H_1$：至少一个总体均值不同。

**平方和公式与检验统计量**：令 $N=\sum_i n_i$，$\bar y_i=\sum_jy_{ij}/n_i$，$\bar y_{\cdot\cdot}=\sum_i\sum_jy_{ij}/N$。

$$SS_A=\sum_i n_i(\bar y_i-\bar y_{\cdot\cdot})^2,\qquad SS_E=\sum_i\sum_j(y_{ij}-\bar y_i)^2,$$
$$SS_T=\sum_i\sum_j(y_{ij}-\bar y_{\cdot\cdot})^2=SS_A+SS_E,$$
$$MS_A=\frac{SS_A}{a-1},\quad MS_E=\frac{SS_E}{N-a},\quad F_0=\frac{MS_A}{MS_E}\overset{H_0}{\sim}F_{a-1,N-a}.$$

当 $F_0>F_{1-\alpha,a-1,N-a}$ 时拒绝 $H_0$，$P=P(F_{a-1,N-a}\ge F_0)$。本文分位数均为**下侧累积概率分位数**。

**点估计量**：$\hat\mu_i=\bar y_i$，$\hat\sigma^2=MS_E$；各组样本方差为 $s_i^2=\sum_j(y_{ij}-\bar y_i)^2/(n_i-1)$。以下为 ANOVA Table 与计算结果。

| 来源 | 自由度 | 平方和 | 均方 | $F$ | $P$ 值 |
| --- | --- | --- | --- | --- | --- |
| 处理 | 3 | 250.65 | 83.55 | 8.36441 | 0.00754584 |
| 误差 | 8 | 79.91 | 9.98875 | — | — |
| 总计 | 11 | 330.56 | — | — | — |

$\hat\sigma^2=MS_E=(330.56-250.65)/8=\mathbf{9.98875}$。

### （b）解释变异的比例

固定效应描述性估计为 $\hat R^2=SS_A/SS_T=250.65/330.56=\mathbf{75.83\%}$。

若按后续题的随机效应解释，$\hat\sigma_\tau^2=(MS_A-MS_E)/n=24.52042$，占比估计量为 $\hat\rho=\hat\sigma_\tau^2/(\hat\sigma_\tau^2+\hat\sigma^2)=\mathbf{71.06\%}$。两种比例定义不同，应随模型说明。

---

**答案来源**：参考：`ch03.pdf`，原题 3.33；本页按第 9 版编号。

← [上一题 3.37](./exercise-3-37.md) · [本章题目目录](../chapter03.md) · [下一题 3.39](./exercise-3-39.md) →
