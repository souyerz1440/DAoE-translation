---
title: 习题 3.41：六水平试验的组内相关区间
short_title: 3.41 六水平试验的组内相关区间
description: 第9版习题3.41的完整中文题目、数据与答案解析
---

(exercise-3-41)=
# 习题 3.41：六水平试验的组内相关区间

[习题与答案](../../index.md) / [第 3 章　单因子试验与方差分析](../chapter03.md)

第 9 版 · 第 3 章 · 原书章末习题见印刷页 125～134。

## 题目

求 3.39 中组内相关系数的 95% 区间。

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
| 处理 | 5 | 750.5 | 150.1 | 24.0561 | 1.26758e-08 |
| 误差 | 24 | 149.75 | 6.23958 | — | — |
| 总计 | 29 | 900.25 | — | — | — |

使用平衡随机效应模型 $y_{ij}=\mu+u_i+\epsilon_{ij}$，各随机效应和误差相互独立且分别正态，同组误差方差相同。$a=6,n=5$，$r=(750.50/5)/(149.75/24)$。

估计目标为 $\rho=\sigma_\tau^2/(\sigma_\tau^2+\sigma^2)$。利用 $r/(1+5\theta)\sim F_{5,24}$，$\theta=\sigma_\tau^2/\sigma^2$，构造

$$L=(r/F_{0.975,5,24}-1)/5,\quad U=(r/F_{0.025,5,24}-1)/5,$$

再变换为 $[L/(1+L),U/(1+U)]\cap[0,1]$。**95% 区间结果为 $(0.56990,0.96775)$**。原答案一处误用每组 3 次重复，本题必须使用 5。

---

**答案来源**：参考：`ch03.pdf`，原题 3.36；本页按第 9 版编号。

**关联题目**：[3.39](./exercise-3-39.md)

← [上一题 3.40](./exercise-3-40.md) · [本章题目目录](../chapter03.md) · [下一题 3.42](./exercise-3-42.md) →
