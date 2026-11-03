---
title: 习题 3.2：六水平试验的 P 值
short_title: 3.2 六水平试验的 P 值
description: 第9版习题3.2的完整中文题目、数据与答案解析
---

(exercise-3-2)=
# 习题 3.2：六水平试验的 P 值

[习题与答案](../../index.md) / [第 3 章　单因子试验与方差分析](../chapter03.md)

第 9 版 · 第 3 章 · 原书章末习题见印刷页 125～134。

## 题目

一位试验者实施了一个单因子试验，该因子有 6 个水平，每个水平重复 3 次。计算得到 F 统计量为 $F_0=5.81$。请给出 $P$ 值的上下界。

## 答案与解析

**模型与模型假设（固定效应）**：$y_{ij}=\mu_i+\epsilon_{ij}=\mu+\tau_i+\epsilon_{ij}$，$i=1,\ldots,a$，$j=1,\ldots,n_i$；$\epsilon_{ij}\overset{\mathrm{iid}}\sim N(0,\sigma^2)$。误差独立、正态且同方差，处理分配及抽样符合题意；不平衡设计使用实际 $n_i$。取 $\sum_i n_i\tau_i=0$ 识别参数；平衡时等价于 $\sum_i\tau_i=0$。

**原假设和备择假设**：$H_0:\mu_1=\cdots=\mu_a$；$H_1$：至少一个总体均值不同。

**平方和公式与检验统计量**：令 $N=\sum_i n_i$，$\bar y_i=\sum_jy_{ij}/n_i$，$\bar y_{\cdot\cdot}=\sum_i\sum_jy_{ij}/N$。

$$SS_A=\sum_i n_i(\bar y_i-\bar y_{\cdot\cdot})^2,\qquad SS_E=\sum_i\sum_j(y_{ij}-\bar y_i)^2,$$
$$SS_T=\sum_i\sum_j(y_{ij}-\bar y_{\cdot\cdot})^2=SS_A+SS_E,$$
$$MS_A=\frac{SS_A}{a-1},\quad MS_E=\frac{SS_E}{N-a},\quad F_0=\frac{MS_A}{MS_E}\overset{H_0}{\sim}F_{a-1,N-a}.$$

当 $F_0>F_{1-\alpha,a-1,N-a}$ 时拒绝 $H_0$，$P=P(F_{a-1,N-a}\ge F_0)$。本文分位数均为**下侧累积概率分位数**。

**点估计量**：$\hat\mu_i=\bar y_i$，$\hat\sigma^2=MS_E$；各组样本方差为 $s_i^2=\sum_j(y_{ij}-\bar y_i)^2/(n_i-1)$。以下为 ANOVA Table 与计算结果。

原题未给原始观测或平方和，无法唯一恢复平方和的绝对数值。令 $m=MS_E>0$，下表保留其尺度，不虚构平方和。

| 来源 | 自由度 | 平方和 | 均方 | $F$ | $P$ 值 |
| --- | --- | --- | --- | --- | --- |
| 处理 | 5 | 29.05m | 5.81m | 5.81 | 0.00594169 |
| 误差 | 12 | 12m | $m$ | — | — |
| 总计 | 17 | 41.05m | — | — | — |

$P=0.00594169$，故 $0.005<P<0.01$。处理与误差自由度为 $(5,12)$。

---

**答案来源**：参考：`ch03.pdf`，原题 3.2；本页按第 9 版编号。

← [上一题 3.1](./exercise-3-1.md) · [本章题目目录](../chapter03.md) · [下一题 3.3](./exercise-3-3.md) →
