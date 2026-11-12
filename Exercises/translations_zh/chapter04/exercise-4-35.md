---
title: 习题 4.35：拉丁方固定效应的最小二乘估计
short_title: 4.35 拉丁方固定效应的最小二乘估计
description: 第9版习题4.35的中文题目、数据与答案解析
---

(exercise-4-35)=
# 习题 4.35：拉丁方固定效应的最小二乘估计

[习题与答案](../../index.md) / [第 4 章　随机区组、拉丁方与相关设计](../chapter04.md)

第 9 版 · 第 4 章 · 原书章末习题见印刷页 171～178。题目依据第 9 版，参考答案按题意匹配。

## 题目

对 $p\times p$ 拉丁方，行效应 $\alpha_i$、列效应 $\beta_k$ 和处理效应 $\tau_j$ 均固定，求这些模型参数的最小二乘估计。

## 答案与解析

**模型与假设**：$y_{ij}=\mu+\tau_{L(i,j)}+\rho_i+\kappa_j+\epsilon_{ij}$，各组效应和为零；各误差独立同分布为 $N(0,\sigma^2)$，所有因子效应可加，没有交互作用。$L(i,j)$ 是该单元的处理标签。行、列及处理均为固定效应。

**原假设与备择假设**：$H_0:\tau_1=\cdots=\tau_p=0$；$H_1$：至少一个处理效应不为零。设 $C=G^2/p^2$，各因子水平总和为 $U_l$，

$$SS_T=\sum y^2-C,\qquad SS_{\text{因子}}=\frac1p\sum_lU_l^2-C,$$
$$SS_E=SS_T-SS_A-SS_R-SS_C,\qquad F_0=\frac{SS_A/(p-1)}{SS_E/[(p-1)(p-2)]}.$$

显著性水平取 $0.05$；处理、行、列自由度各为 $p-1$，总自由度 $p^2-1$。


在 $\sum_i\alpha_i=\sum_k\beta_k=\sum_j\tau_j=0$ 下，因每个处理在各行各列各出现一次，正规方程为

$$p^2\hat\mu=G,\quad p\hat\mu+p\hat\alpha_i=R_i,\quad p\hat\mu+p\hat\beta_k=C_k,\quad p\hat\mu+p\hat\tau_j=T_j.$$

故 $\hat\mu=G/p^2$，$\hat\alpha_i=R_i/p-\hat\mu$，$\hat\beta_k=C_k/p-\hat\mu$，$\hat\tau_j=T_j/p-\hat\mu$。这是点估计公式；题目无数值数据，不虚构 ANOVA 数值。

---

**答案来源**：按题意对应 `Exercises/ch04.pdf` 的 4.32；公式与数值由编者复核，补充完整推导与诊断。

← [上一题 4.34](./exercise-4-34.md) · [本章题目目录](../chapter04.md) · [下一题 4.36](./exercise-4-36.md) →
