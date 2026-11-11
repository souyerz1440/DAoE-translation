---
title: 习题 4.37：多个拉丁方的联合分析
short_title: 4.37 多个拉丁方的联合分析
description: 第9版习题4.37的中文题目、数据与答案解析
---

(exercise-4-37)=
# 习题 4.37：多个拉丁方的联合分析

[习题与答案](../../index.md) / [第 4 章　随机区组、拉丁方与相关设计](../chapter04.md)

第 9 版 · 第 4 章 · 原书章末习题见印刷页 171～178。题目依据第 9 版，参考答案按题意匹配。

## 题目

为了增加处理重复，运行 $n$ 个 $p\times p$ 拉丁方。模型为

$$y_{ijkh}=\mu+\rho_h+\alpha_{i(h)}+\tau_j+\beta_{k(h)}+(\tau\rho)_{jh}+\epsilon_{ijkh},$$

其中 $\rho_h$ 为方效应，行列嵌套于方。采用 $\sum_h\rho_h=0$，每方行、列效应分别和为零，$\sum_j\tau_j=0$，交互作用按处理及方分别求和均为零。

**（a）** 写出正规方程并求参数估计。

**（b）** 写出方差分析表。


## 答案与解析

### （a）模型、正规方程与估计

误差独立同分布为 $N(0,\sigma^2)$；除显式加入的处理×方交互作用外，不含其他交互作用。所有效应暂按固定效应解释。令 $G_h,R_{ih},C_{kh},T_j,T_{jh}$ 分别为方、方内行、方内列、处理、处理×方总和。正规方程可写为

$$np^2\hat\mu=G,\quad p^2(\hat\mu+\hat\rho_h)=G_h,$$
$$p(\hat\mu+\hat\rho_h+\hat\alpha_{i(h)})=R_{ih},\quad p(\hat\mu+\hat\rho_h+\hat\beta_{k(h)})=C_{kh},$$
$$np(\hat\mu+\hat\tau_j)=T_j,\quad p[\hat\mu+\hat\rho_h+\hat\tau_j+\widehat{(\tau\rho)}_{jh}]=T_{jh}.$$

因此 $\hat\mu=\bar y$，$\hat\rho_h=\bar y_h-\bar y$，$\hat\tau_j=\bar y_j-\bar y$，$\hat\alpha_{i(h)}=\bar y_{ih}-\bar y_h$，$\hat\beta_{k(h)}=\bar y_{kh}-\bar y_h$，$\widehat{(\tau\rho)}_{jh}=\bar y_{jh}-\bar y_j-\bar y_h+\bar y$。


### （b）符号 ANOVA 表

设 $C=G^2/(np^2)$，$SS_T=\sum y^2-C$，各平方和及自由度如下。均方等于平方和除以自由度。


| 来源 | 平方和 | 自由度 | 固定效应 F |
| --- | --- | --- | --- |
| 处理 | $\sum_jT_j^2/(np)-C$ | $p-1$ | $MS_A/MS_E$ |
| 方 | $\sum_hG_h^2/p^2-C$ | $n-1$ | $MS_S/MS_E$ |
| 行（方内） | $\sum_{ih}R_{ih}^2/p-\sum_hG_h^2/p^2$ | $n(p-1)$ | $MS_R/MS_E$ |
| 列（方内） | $\sum_{kh}C_{kh}^2/p-\sum_hG_h^2/p^2$ | $n(p-1)$ | $MS_C/MS_E$ |
| 处理×方 | $\sum_{jh}T_{jh}^2/p-C-SS_A-SS_S$ | $(p-1)(n-1)$ | $MS_{AS}/MS_E$ |
| 误差 | 由总平方和减去以上各项 | $n(p-1)(p-2)$ | — |
| 总计 | $\sum y^2-C$ | $np^2-1$ | — |

处理检验 $H_0:\tau_j=0$ 对至少一个非零；交互作用检验 $H_0:(\tau\rho)_{jh}=0$ 对至少一个非零。若方为随机抽取并把处理×方也作为随机项，处理检验分母应改为 $MS_{AS}$，自由度为 $(p-1,(p-1)(n-1))$。有交互作用时还应解释各方中的处理表现。

---

**答案来源**：按题意对应 `Exercises/ch04.pdf` 的 4.34；公式与数值由编者复核，补充完整推导与诊断。

← [上一题 4.36](./exercise-4-36.md) · [本章题目目录](../chapter04.md) · [下一题 4.38](./exercise-4-38.md) →
