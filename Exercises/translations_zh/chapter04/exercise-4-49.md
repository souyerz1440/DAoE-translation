---
title: 习题 4.49：BIBD 调整处理平方和的证明
short_title: 4.49 BIBD 调整处理平方和的证明
description: 第9版习题4.49的中文题目、数据与答案解析
---

(exercise-4-49)=
# 习题 4.49：BIBD 调整处理平方和的证明

[习题与答案](../../index.md) / [第 4 章　随机区组、拉丁方与相关设计](../chapter04.md)

第 9 版 · 第 4 章 · 原书章末习题见印刷页 171～178。题目依据第 9 版，参考答案按题意匹配。

## 题目

证明 BIBD 的调整处理平方和为 $k\sum_iQ_i^2/(\lambda a)$。

## 答案与解析

设处理—区组关联矩阵为 $N$，$T$、$B$ 是处理、区组总和向量。消去区组后的正规方程为 $C\hat\tau=Q$，其中

$$C=rI-\frac1kNN^T,\quad Q=T-\frac1kNB.$$

BIBD 中 $NN^T=(r-\lambda)I+\lambda J$，并有 $r(k-1)=\lambda(a-1)$，故 $C=(\lambda a/k)(I-J/a)$。又 $\mathbf1^TQ=0$，在 $\sum_i\tau_i=0$ 下 $\hat\tau=kQ/(\lambda a)$。因此

$$SS_{A\mid B}=\hat\tau^TQ=\frac{k}{\lambda a}\sum_iQ_i^2.$$

该量等于加入处理后减少的残差平方和，自由度为 $a-1$。

---

**答案来源**：按题意对应 `Exercises/ch04.pdf` 的 4.46S；公式与数值由编者复核，补充完整推导与诊断。

← [上一题 4.48](./exercise-4-48.md) · [本章题目目录](../chapter04.md) · [下一题 4.50](./exercise-4-50.md) →
