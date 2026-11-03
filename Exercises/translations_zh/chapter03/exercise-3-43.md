---
title: 习题 3.43：误差方差占比的区间
short_title: 3.43 误差方差占比的区间
description: 第9版习题3.43的完整中文题目、数据与答案解析
---

(exercise-3-43)=
# 习题 3.43：误差方差占比的区间

[习题与答案](../../index.md) / [第 3 章　单因子试验与方差分析](../chapter03.md)

第 9 版 · 第 3 章 · 原书章末习题见印刷页 125～134。

## 题目

对平衡单因子随机效应模型，推导 $\sigma^2/(\sigma_\tau^2+\sigma^2)$ 的 $100(1-\alpha)\%$ 区间。

## 答案与解析

在平衡随机效应模型中，$y_{ij}=\mu+u_i+\epsilon_{ij}$，$u_i\sim N(0,\sigma_\tau^2)$、$\epsilon_{ij}\sim N(0,\sigma^2)$，相互独立。目标为 $\eta=\sigma^2/(\sigma_\tau^2+\sigma^2)=1-\rho$。

令 $r=MS_A/MS_E$，$\theta=\sigma_\tau^2/\sigma^2$，由 $r/(1+n\theta)\sim F_{a-1,N-a}$，先求

$$L=(r/F_{1-\alpha/2,a-1,N-a}-1)/n,\qquad U=(r/F_{\alpha/2,a-1,N-a}-1)/n.$$

对递减函数 $\eta=1/(1+\theta)$ 变换，所求 **$100(1-\alpha)\%$ 区间结果**为

$$\boxed{\left[\frac1{1+U},\frac1{1+L}\right]\cap[0,1]}.$$

也可先得到 $\rho$ 的区间 $[L_\rho,U_\rho]$，再取 $[1-U_\rho,1-L_\rho]$。本题只求一般式，没有具体数据，结果为上述符号区间。

---

**答案来源**：参考答案文件无直接对应题；解析由编者补充。

**关联题目**：[3.40](./exercise-3-40.md)

← [上一题 3.42](./exercise-3-42.md) · [本章题目目录](../chapter03.md) · [下一题 3.44](./exercise-3-44.md) →
