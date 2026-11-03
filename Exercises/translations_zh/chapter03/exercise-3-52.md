---
title: 习题 3.52：已知四个总体均值的功效规划
short_title: 3.52 已知四个总体均值的功效规划
description: 第9版习题3.52的完整中文题目、数据与答案解析
---

(exercise-3-52)=
# 习题 3.52：已知四个总体均值的功效规划

[习题与答案](../../index.md) / [第 3 章　单因子试验与方差分析](../chapter03.md)

第 9 版 · 第 3 章 · 原书章末习题见印刷页 125～134。

## 题目

假设四个正态总体的均值分别为 $\mu_1=50$、$\mu_2=60$、$\mu_3=50$、$\mu_4=60$。希望拒绝“各总体均值相等”这一原假设的概率至少为 0.90，每个总体应抽取多少个观测？假定 $\alpha=0.05$，且误差方差的合理估计为 $\sigma^2=25$。

## 答案与解析

**规划用模型、假设与检验结构**

**模型与模型假设（固定效应）**：$y_{ij}=\mu_i+\epsilon_{ij}=\mu+\tau_i+\epsilon_{ij}$，$i=1,\ldots,a$，$j=1,\ldots,n_i$；$\epsilon_{ij}\overset{\mathrm{iid}}\sim N(0,\sigma^2)$。误差独立、正态且同方差，处理分配及抽样符合题意；不平衡设计使用实际 $n_i$。取 $\sum_i n_i\tau_i=0$ 识别参数；平衡时等价于 $\sum_i\tau_i=0$。

**原假设和备择假设**：$H_0:\mu_1=\cdots=\mu_a$；$H_1$：至少一个总体均值不同。

**平方和公式与检验统计量**：令 $N=\sum_i n_i$，$\bar y_i=\sum_jy_{ij}/n_i$，$\bar y_{\cdot\cdot}=\sum_i\sum_jy_{ij}/N$。

$$SS_A=\sum_i n_i(\bar y_i-\bar y_{\cdot\cdot})^2,\qquad SS_E=\sum_i\sum_j(y_{ij}-\bar y_i)^2,$$
$$SS_T=\sum_i\sum_j(y_{ij}-\bar y_{\cdot\cdot})^2=SS_A+SS_E,$$
$$MS_A=\frac{SS_A}{a-1},\quad MS_E=\frac{SS_E}{N-a},\quad F_0=\frac{MS_A}{MS_E}\overset{H_0}{\sim}F_{a-1,N-a}.$$

当 $F_0>F_{1-\alpha,a-1,N-a}$ 时拒绝 $H_0$，$P=P(F_{a-1,N-a}\ge F_0)$。本文分位数均为**下侧累积概率分位数**。

**点估计量**：$\hat\mu_i=\bar y_i$，$\hat\sigma^2=MS_E$；各组样本方差为 $s_i^2=\sum_j(y_{ij}-\bar y_i)^2/(n_i-1)$。以下为 ANOVA Table 与计算结果。

样本量规划尚未产生观测值，实际平方和未知，以下 ANOVA Table 仅列结构；不能把非中心参数当作已观测平方和。

| 来源 | 自由度 | 平方和 | 均方 |
| --- | --- | --- | --- |
| 处理 | 3 | $SS_A$ | $SS_A/3$ |
| 误差 | $4(n-1)$ | $SS_E$ | $SS_E/[4(n-1)]$ |
| 总计 | $4n-1$ | $SS_T$ | — |

总均值 55，$\sum_i(\mu_i-55)^2=100$，非中心参数 $\lambda=100n/25=4n$。计算

$$P\{F_{3,4(n-1)}(\lambda)>F_{0.95,3,4(n-1)}\}\ge0.90$$

得到每总体至少 **5 个**，总数 20。原答案给出 9 个，对“已指定四个均值”的本题偏保守；它接近仅给最大均值差时采用最不利均值配置的结果，二者不是同一个规划问题。

---

**答案来源**：参考：`ch03.pdf`，原题 3.45；本页按第 9 版编号。

← [上一题 3.51](./exercise-3-51.md) · [本章题目目录](../chapter03.md) · [下一题 3.53](./exercise-3-53.md) →
