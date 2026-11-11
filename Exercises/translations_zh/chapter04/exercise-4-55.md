---
title: 习题 4.55：BIBD 组内效应估计量的方差
short_title: 4.55 BIBD 组内效应估计量的方差
description: 第9版习题4.55的中文题目、数据与答案解析
---

(exercise-4-55)=
# 习题 4.55：BIBD 组内效应估计量的方差

[习题与答案](../../index.md) / [第 4 章　随机区组、拉丁方与相关设计](../chapter04.md)

第 9 版 · 第 4 章 · 原书章末习题见印刷页 171～178。题目依据第 9 版，参考答案按题意匹配。

## 题目

证明组内处理效应估计量的方差为 $\operatorname{Var}(\hat\tau_i)=k(a-1)\sigma^2/(\lambda a^2)$。

## 答案与解析

由区组消元的正规方程，$C=(\lambda a/k)(I-J/a)$，$\hat\tau=C^+Q$，且 $\operatorname{Var}(Q)=\sigma^2C$，所以

$$\operatorname{Cov}(\hat\tau)=\sigma^2C^+=\frac{k\sigma^2}{\lambda a}(I-J/a).$$

对角元为 $k(a-1)\sigma^2/(\lambda a^2)$；非对角元为 $-k\sigma^2/(\lambda a^2)$。因此 $\operatorname{Var}(\hat\tau_i-\hat\tau_j)=2k\sigma^2/(\lambda a)$，各处理效应估计并非相互独立。

---

**答案来源**：按题意对应 `Exercises/ch04.pdf` 的 4.52S；公式与数值由编者复核，补充完整推导与诊断。

← [上一题 4.54](./exercise-4-54.md) · [本章题目目录](../chapter04.md) · [下一题 4.56](./exercise-4-56.md) →
