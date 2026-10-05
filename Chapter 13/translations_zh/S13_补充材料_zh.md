# 第 13 章补充材料

## S13.1 随机模型的期望均方

考虑正文式（13.1）的两因子随机效应平衡方差分析模型：

$$
y_{ijk}=\mu+\tau_i+\beta_j+(\tau\beta)_{ij}+\epsilon_{ijk},\qquad i=1,\ldots,a,\quad j=1,\ldots,b,\quad k=1,\ldots,n.
$$

正文式（13.3）列出了期望均方，但没有正式推导。直接应用期望算子即可求得。例如，

$$
E(MS_A)=E\left(\frac{SS_A}{a-1}\right)=\frac{E(SS_A)}{a-1}.
$$

$SS_A$ 是行因子的平方和。模型中的 $\tau_i,\beta_j,(\tau\beta)_{ij}$ 独立服从均值为零、方差分别为 $\sigma_\tau^2,\sigma_\beta^2,\sigma_{\tau\beta}^2$ 的正态分布。平方和及其期望为

$$
SS_A=\frac1{bn}\sum_{i=1}^ay_{i\cdot\cdot}^2-\frac{y_{\cdots}^2}{abn},\qquad
E(SS_A)=\frac1{bn}E\left(\sum_{i=1}^ay_{i\cdot\cdot}^2\right)-\frac1{abn}E(y_{\cdots}^2).
$$

由于

$$
y_{i\cdot\cdot}=bn\mu+bn\tau_i+n\beta_\cdot+n(\tau\beta)_{i\cdot}+\epsilon_{i\cdot\cdot},
$$

利用随机分量均值为零以及相互独立的性质，有

$$
\begin{aligned}
\frac1{bn}E\left(\sum_i y_{i\cdot\cdot}^2\right)
&=\frac1{bn}\left[a(bn\mu)^2+a(bn)^2\sigma_\tau^2+abn^2\sigma_\beta^2+abn^2\sigma_{\tau\beta}^2+abn\sigma^2\right]\\
&=abn\mu^2+abn\sigma_\tau^2+an\sigma_\beta^2+an\sigma_{\tau\beta}^2+a\sigma^2.
\end{aligned}
$$

另一方面，

$$
y_{\cdots}=abn\mu+bn\tau_\cdot+an\beta_\cdot+n(\tau\beta)_{\cdot\cdot}+\epsilon_{\cdots},
$$

所以

$$
\begin{aligned}
\frac1{abn}E(y_{\cdots}^2)
&=\frac1{abn}\left[(abn\mu)^2+a(bn)^2\sigma_\tau^2+b(an)^2\sigma_\beta^2+abn^2\sigma_{\tau\beta}^2+abn\sigma^2\right]\\
&=abn\mu^2+bn\sigma_\tau^2+an\sigma_\beta^2+n\sigma_{\tau\beta}^2+\sigma^2.
\end{aligned}
$$

收集各项并除以 $a-1$，得到

$$
\begin{aligned}
E(MS_A)&=\frac{(a-1)\sigma^2+n(a-1)\sigma_{\tau\beta}^2+bn(a-1)\sigma_\tau^2}{a-1}\\
&=\sigma^2+n\sigma_{\tau\beta}^2+bn\sigma_\tau^2,
\end{aligned}
$$

与正文式（13.3）的第一个结果相同。

## S13.2 混合模型的期望均方

正文第 13.3 节介绍了若干混合模型，其期望均方依赖于模型假定。这里考虑受限模型，下一节讨论非受限模型。受限模型对固定因子 A 的约束为

$$
\tau_\cdot=0,\qquad(\tau\beta)_{\cdot j}=0,\qquad
V[(\tau\beta)_{ij}]=\frac{a-1}{a}\sigma_{\tau\beta}^2.
$$

求随机因子 B 的期望均方：

$$
E(MS_B)=\frac{E(SS_B)}{b-1},\qquad
E(SS_B)=\frac1{an}E\left(\sum_{j=1}^by_{\cdot j\cdot}^2\right)-\frac1{abn}E(y_{\cdots}^2).
$$

利用约束，

$$
y_{\cdot j\cdot}=an\mu+an\beta_j+\epsilon_{\cdot j\cdot},
$$

因而

$$
\frac1{an}E\left(\sum_jy_{\cdot j\cdot}^2\right)
=\frac{b(an\mu)^2+b(an)^2\sigma_\beta^2+abn\sigma^2}{an}
=abn\mu^2+abn\sigma_\beta^2+b\sigma^2.
$$

又因 $y_{\cdots}=abn\mu+an\beta_\cdot+\epsilon_{\cdots}$，有

$$
\frac1{abn}E(y_{\cdots}^2)
=\frac{(abn\mu)^2+b(an)^2\sigma_\beta^2+abn\sigma^2}{abn}
=abn\mu^2+an\sigma_\beta^2+\sigma^2.
$$

所以

$$
E(MS_B)=\frac{(b-1)\sigma^2+an(b-1)\sigma_\beta^2}{b-1}
=\sigma^2+an\sigma_\beta^2.
$$

其他期望均方可类似推导。

## S13.3 受限与非受限混合模型

考虑非受限模型

$$
y_{ijk}=\mu+\alpha_i+\gamma_j+(\alpha\gamma)_{ij}+\epsilon_{ijk},\qquad
i=1,\ldots,a,\quad j=1,\ldots,b,\quad k=1,\ldots,n,
$$

假定 $\alpha_\cdot=0$、$V[(\alpha\gamma)_{ij}]=\sigma_{\alpha\gamma}^2$，且所有随机效应互不相关。与受限模型不同，不要求交互作用对固定因子水平求和为零。受限模型允许更一般的相关结构，但现代软件有的提供两种选择，有的只采用非受限模型，因此两种形式都值得关注。

下面推导随机因子 B 的期望均方。它与受限模型不同，而差异的关键正是交互作用假定。仍有

$$
E(MS_B)=\frac1{b-1}\left[\frac1{an}E\left(\sum_jy_{\cdot j\cdot}^2\right)-\frac1{abn}E(y_{\cdots}^2)\right].
$$

先考虑

$$
\begin{aligned}
y_{\cdot j\cdot}&=an\mu+n\alpha_\cdot+an\gamma_j+n(\alpha\gamma)_{\cdot j}+\epsilon_{\cdot j\cdot}\\
&=an\mu+an\gamma_j+n(\alpha\gamma)_{\cdot j}+\epsilon_{\cdot j\cdot}.
\end{aligned}
$$

$\alpha_\cdot=0$，但交互作用和不为零。于是

$$
\frac1{an}E\left(\sum_jy_{\cdot j\cdot}^2\right)
=abn\mu^2+abn\sigma_\gamma^2+bn\sigma_{\alpha\gamma}^2+b\sigma^2.
$$

由于

$$
y_{\cdots}=abn\mu+an\gamma_\cdot+n(\alpha\gamma)_{\cdot\cdot}+\epsilon_{\cdots},
$$

有

$$
\frac1{abn}E(y_{\cdots}^2)=abn\mu^2+an\sigma_\gamma^2+n\sigma_{\alpha\gamma}^2+\sigma^2.
$$

相减并除以 $b-1$，得到

$$
E(MS_B)=\frac{(b-1)\sigma^2+n(b-1)\sigma_{\alpha\gamma}^2+an(b-1)\sigma_\gamma^2}{b-1}
=\sigma^2+n\sigma_{\alpha\gamma}^2+an\sigma_\gamma^2,
$$

与正文式（13.12）一致。直接使用期望算子推导较繁琐，正文的规则可节省很多工作。还有适用于不平衡设计的其他规则和算法，参见 Milliken 和 Johnson（1984）。

## S13.4 样本量不等的随机模型和混合模型

若各单元观测数不同，方差分析通常更复杂。第 15 章简要讨论两因子固定效应不平衡设计，这里对随机或混合模型提供一些建议。

不平衡随机或混合模型通常没有平衡设计中那样的精确 F 检验，正文基于独立均方的 Satterthwaite 合成检验公式也不能直接套用。较方便的方法是极大似然或 REML，正文第 13.2 节介绍的方法及 SAS、JMP 均可用于不平衡设计。其局限是方差分量推断常依赖大样本近似，而设计试验的次数通常不多。一般讨论见 Searle（1987）。

## S13.5 修正大样本方法的背景

正文第 13.6.2 节讨论了可写为均方线性组合的方差分量的置信区间。大样本理论表明，当 $\min(f_1,\ldots,f_Q)\to\infty$ 时，

$$
Z=\frac{\hat\sigma_0^2-\sigma_0^2}{\sqrt{V(\hat\sigma_0^2)}}\ \xrightarrow{d}\ N(0,1),\qquad
V(\hat\sigma_0^2)=2\sum_{i=1}^Q\frac{c_i^2\theta_i^2}{f_i}.
$$

$\theta_i=E(MS_i)$ 是第 i 个均方所估计的方差分量线性组合，$f_i$ 是其自由度。因此，双侧 $100(1-\alpha)\%$ 大样本置信区间为

$$
\hat\sigma_0^2-z_{\alpha/2}\sqrt{V(\hat\sigma_0^2)}
\le\sigma_0^2\le
\hat\sigma_0^2+z_{\alpha/2}\sqrt{V(\hat\sigma_0^2)}.
$$

实际计算用 $MS_i$ 替换 $\theta_i$。这也是正文中 JMP 构造近似区间的依据。自由度大时效果较好，自由度小时可能不可靠。Welch（1956）提出了改善方法；Graybill 和 Wang（1980）的修正使某些特殊情形的区间成为精确区间，在其他情形中也常有较好的近似表现，其结果见正文式（13.25）。

## S13.6 用修正大样本方法构造方差分量比值的置信区间

实际常关注方差分量的比值。例如，正文例 13.1 的测量系统能力研究中，量具总方差为 $\sigma_\beta^2+\sigma_{\tau\beta}^2+\sigma^2$，产品方差为 $\sigma_\tau^2$。为将量具变异表示为产品变异的百分比，可考察

$$
\frac{\sigma_\beta^2+\sigma_{\tau\beta}^2+\sigma^2}{\sigma_\tau^2}.
$$

设所关注比值为 $\sigma_1^2/\sigma_2^2$，两方差可分别用均方线性组合估计：

$$
\frac{\hat\sigma_1^2}{\hat\sigma_2^2}
=\frac{\sum_{i=1}^Pc_iMS_i}{\sum_{j=P+1}^Qc_jMS_j}.
$$

其 $100(1-\alpha)\%$ 单侧置信下限为

$$
L=\frac{\hat\sigma_1^2}{\hat\sigma_2^2}
\left[\frac{2+k_4/(k_1k_2)-\sqrt{V_L}}{2(1-k_5/k_2^2)}\right],
$$

其中

$$
\begin{aligned}
V_L&=\left(2+\frac{k_4}{k_1k_2}\right)^2-4\left(1-\frac{k_5}{k_2^2}\right)\left(1-\frac{k_3}{k_1^2}\right),\\
k_1&=\sum_{i=1}^Pc_iMS_i,\qquad k_2=\sum_{j=P+1}^Qc_jMS_j,\\
k_3&=\sum_{i=1}^PG_i^2c_i^2MS_i^2+\sum_{i=1}^{P-1}\sum_{t>i}^PG_{it}^*c_ic_tMS_iMS_t,\\
k_4&=\sum_{i=1}^P\sum_{j=P+1}^QG_{ij}c_ic_jMS_iMS_j,\\
k_5&=\sum_{j=P+1}^QH_j^2c_j^2MS_j^2.
\end{aligned}
$$

$G_i,H_j,G_{ij},G_{it}^*$ 的定义沿用正文修正大样本方法中的量。更多细节见 Burdick 和 Graybill（1992）。原补充材料中若干公式及交叉引用沿用旧版编号，这里按当前正文调整；比值公式中的分母也按上下文恢复为 $\hat\sigma_2^2$。[^1]

## 补充参考文献

Graybill, F. A. and C. M. Wang（1980），“Confidence Intervals on Nonnegative Linear Combinations of Variances”，*Journal of the American Statistical Association*，75，869–873。

Welch, B. L.（1956），“On Linear Combinations of Several Variances”，*Journal of the American Statistical Association*，51，132–148。

[^1]: 译者注：原补充材料的受限模型交互作用方差比例印为 $a/(a-1)$，应与正文一致为 $(a-1)/a$；随机模型的 $E(SS_A)$ 化简中漏了 $\sigma_\tau^2$ 项的 $(a-1)$ 因子；S13.6 置信下限的前置比值分母重复印为 $\hat\sigma_1^2$，应为 $\hat\sigma_2^2$。这里均按推导修正。
