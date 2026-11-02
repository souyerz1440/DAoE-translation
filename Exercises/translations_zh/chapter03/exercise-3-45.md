---
title: 习题 3.45：线性组合的方差
short_title: 3.45 线性组合的方差
description: 第9版习题3.45的完整中文题目、数据与答案解析
---

(exercise-3-45)=
# 习题 3.45：线性组合的方差

[习题与答案](../../index.md) / [第 3 章　单因子试验与方差分析](../chapter03.md)

第 9 版 · 第 3 章 · 原书章末习题见印刷页 125～134。

## 题目

证明处理观测总和的线性组合 $\sum_i c_i y_{i\cdot}$ 的方差为 $\sigma^2\sum_i n_i c_i^2$。

## 答案与解析

在固定效应、独立同方差误差下，$y_{i\cdot}=\sum_j y_{ij}$，所以 $\operatorname{Var}(y_{i\cdot})=n_i\sigma^2$，不同组协方差为零。故

$$\operatorname{Var}\left(\sum_i c_i y_{i\cdot}\right)=\sum_i c_i^2 n_i\sigma^2.$$

注意本题是**组总和**；若使用组均值 $\bar y_i$，相应公式是 $\sigma^2\sum_i c_i^2/n_i$。随机效应模型的无条件组间方差还包含 $\sigma_\tau^2$，不能直接套用。

---

**答案来源**：参考：`ch03.pdf`，原题 3.38S；本页按第 9 版编号。

← [上一题 3.44](./exercise-3-44.md) · [本章题目目录](../chapter03.md) · [下一题 3.46](./exercise-3-46.md) →
