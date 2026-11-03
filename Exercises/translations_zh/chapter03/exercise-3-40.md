---
title: 习题 3.40：四水平试验的组内相关区间
short_title: 3.40 四水平试验的组内相关区间
description: 第9版习题3.40的完整中文题目、数据与答案解析
---

(exercise-3-40)=
# 习题 3.40：四水平试验的组内相关区间

[习题与答案](../../index.md) / [第 3 章　单因子试验与方差分析](../chapter03.md)

第 9 版 · 第 3 章 · 原书章末习题见印刷页 125～134。

## 题目

求 3.38 中组内相关系数的 95% 区间。

## 答案与解析

**共同模型、假设与方差估计依据**

**模型与模型假设（随机效应）**：$y_{ij}=\mu+u_i+\epsilon_{ij}$，$u_i\overset{\mathrm{iid}}\sim N(0,\sigma_\tau^2)$，$\epsilon_{ij}\overset{\mathrm{iid}}\sim N(0,\sigma^2)$；随机效应与误差相互独立，抽取的水平代表目标水平总体。不同组独立，同组观测的边际协方差为 $\sigma_\tau^2$，不能把同组观测的边际独立性当作假设。

**原假设和备择假设**：$H_0:\sigma_\tau^2=0$；$H_1:\sigma_\tau^2>0$。

**平方和公式与检验统计量**：令 $N=\sum_i n_i$，$\bar y_i=\sum_jy_{ij}/n_i$，$\bar y_{\cdot\cdot}=\sum_i\sum_jy_{ij}/N$。

$$SS_A=\sum_i n_i(\bar y_i-\bar y_{\cdot\cdot})^2,\qquad SS_E=\sum_i\sum_j(y_{ij}-\bar y_i)^2,$$
$$SS_T=\sum_i\sum_j(y_{ij}-\bar y_{\cdot\cdot})^2=SS_A+SS_E,$$
$$MS_A=\frac{SS_A}{a-1},\quad MS_E=\frac{SS_E}{N-a},\quad F_0=\frac{MS_A}{MS_E}\overset{H_0}{\sim}F_{a-1,N-a}.$$

当 $F_0>F_{1-\alpha,a-1,N-a}$ 时拒绝 $H_0$，$P=P(F_{a-1,N-a}\ge F_0)$。本文分位数均为**下侧累积概率分位数**。

**方差分量估计量**：$\hat\sigma^2=MS_E$，$\hat\sigma_\tau^2=(MS_A-MS_E)/n_0$，其中 $n_0=[N-\sum_i n_i^2/N]/(a-1)$，平衡时 $n_0=n$。表中组样本均值仅作描述，不是随机效应的 BLUP；各组样本方差为 $s_i^2=\sum_j(y_{ij}-\bar y_i)^2/(n_i-1)$。以下为 ANOVA Table 与计算结果。

| 来源 | 自由度 | 平方和 | 均方 | $F$ | $P$ 值 |
| --- | --- | --- | --- | --- | --- |
| 处理 | 3 | 250.65 | 83.55 | 8.36441 | 0.00754584 |
| 误差 | 8 | 79.91 | 9.98875 | — | — |
| 总计 | 11 | 330.56 | — | — | — |

采用平衡单因子随机效应模型 $y_{ij}=\mu+u_i+\epsilon_{ij}$，$u_i\sim N(0,\sigma_\tau^2)$、$\epsilon_{ij}\sim N(0,\sigma^2)$，两者独立。$a=4,n=3$，$r=MS_A/MS_E=(250.65/3)/(79.91/8)$。

令 $\theta=\sigma_\tau^2/\sigma^2$，枢轴量为 $r/(1+n\theta)\sim F_{3,8}$。先求

$$L=(r/F_{0.975,3,8}-1)/3,\quad U=(r/F_{0.025,3,8}-1)/3,$$

再用 $\rho=\theta/(1+\theta)$，把 $[L/(1+L),U/(1+U)]$ 与 $[0,1]$ 相交。**95% 区间结果为 $(0.15359,0.97573)$**。分子是组间方差分量，不是误差方差。

---

**答案来源**：参考：`ch03.pdf`，原题 3.35；本页按第 9 版编号。

**关联题目**：[3.38](./exercise-3-38.md)

← [上一题 3.39](./exercise-3-39.md) · [本章题目目录](../chapter03.md) · [下一题 3.41](./exercise-3-41.md) →
