---
title: 习题 3.7：补全 ANOVA 表（一）
short_title: 3.7 补全 ANOVA 表（一）
description: 第9版习题3.7的完整中文题目、数据与答案解析
---

(exercise-3-7)=
# 习题 3.7：补全 ANOVA 表（一）

[习题与答案](../../index.md) / [第 3 章　单因子试验与方差分析](../chapter03.md)

第 9 版 · 第 3 章 · 原书章末习题见印刷页 125～134。

## 题目

下面给出了部分单因子方差分析软件输出。请补全所有问号处的数值；对于 $P$ 值，可以给出上下界。

| 来源 | 自由度 | 平方和 | 均方 | $F$ | $P$ |
| --- | --- | --- | --- | --- | --- |
| 因子 | 3 | 36.15 | ? | ? | ? |
| 误差 | ? | ? | ? | — | — |
| 总计 | 19 | 196.04 | — | — | — |

## 答案与解析

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
| 处理 | 3 | 36.15 | 12.05 | 1.20583 | 0.339523 |
| 误差 | 16 | 159.89 | 9.99312 | — | — |
| 总计 | 19 | 196.04 | — | — | — |

补全依据为 $df_E=df_T-df_A=16$，$SS_E=196.04-36.15=159.89$。$P=0.339523$，在 0.05 水平不拒绝各均值相等。这里 $N=20,a=4$；汇总表本身不能证明每组都恰为 5 次重复。

---

**答案来源**：参考答案文件无直接对应题；解析由编者补充。

← [上一题 3.6](./exercise-3-6.md) · [本章题目目录](../chapter03.md) · [下一题 3.8](./exercise-3-8.md) →
