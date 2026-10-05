# 第 10 章补充材料

## S10.1 回归系数的协方差矩阵

教材第 10.3 节已说明，线性回归模型 $\mathbf y=\mathbf X\boldsymbol\beta+\boldsymbol\varepsilon$ 的最小二乘估计量

$$
\hat{\boldsymbol\beta}=(\mathbf X'\mathbf X)^{-1}\mathbf X'\mathbf y
$$

是无偏的，并给出其协方差矩阵为 $\sigma^2(\mathbf X'\mathbf X)^{-1}$〔式（10.19）〕。这一结果可直接推导。由于 $(\mathbf X'\mathbf X)^{-1}\mathbf X'$ 是常数矩阵，而 $\mathbf y$ 是随机向量，利用标量随机变量乘以常数时方差乘以常数平方的矩阵形式，有

$$
\begin{aligned}
\operatorname{Var}(\hat{\boldsymbol\beta})
&=(\mathbf X'\mathbf X)^{-1}\mathbf X'\operatorname{Var}(\mathbf y)
  \big[(\mathbf X'\mathbf X)^{-1}\mathbf X'\big]'\\
&=\sigma^2(\mathbf X'\mathbf X)^{-1}\mathbf X'\mathbf X
  (\mathbf X'\mathbf X)^{-1}\\
&=\sigma^2(\mathbf X'\mathbf X)^{-1}.
\end{aligned}
$$

其中 $\operatorname{Var}(\mathbf y)=\sigma^2\mathbf I$；推导还用到了矩阵乘积的转置要按相反顺序排列各因子，以及 $(\mathbf X'\mathbf X)^{-1}$ 的对称性。


## S10.2 回归模型与设计试验

教材例 10.2～10.5 展示了回归方法分析设计试验数据的几种用法。例 10.2 对有四次中心点试验的 $2^3$ 因子设计拟合主效应模型。由于设计正交，$(\mathbf X'\mathbf X)^{-1}$ 为对角矩阵，各回归系数估计量之间的协方差均为零，其方差为

$$
\operatorname{Var}(\hat\beta_0)=\frac{\sigma^2}{12}=0.0833\sigma^2,
\qquad
\operatorname{Var}(\hat\beta_i)=\frac{\sigma^2}{8}=0.125\sigma^2,
\quad i=1,2,3.
$$

例 10.3 考察同一问题，但假设原来的 12 个观测中缺失一个。用其余 11 个观测拟合一阶模型后，系数估计值变化不大；不过，$(\mathbf X'\mathbf X)^{-1}$ 表明缺失观测使各系数估计量的方差增大，并使不同系数估计量之间出现协方差。

例 10.4 研究因子水平不准确的影响，其系数估计量也不再互不相关。不过，偏离正交设计并不意味着各系数的方差必然增大，还须考虑实际因子水平及各列的平方和。按该例的编码水平复算，三个斜率系数的相对方差约为 0.12289、0.11753、0.10845，均小于例 10.2 的 0.125；截距的相对方差则略有增大。这两个例子中的协方差并不算很大，通常不会妨碍结果解释。

## S10.3 调整后的 $R^2$

教材多次指出，调整后的 $R^2$ 比普通 $R^2$ 更适合比较不同大小的模型，因为它不会随模型变量数增加而必然不减。由式（10.27），

$$
R_{\mathrm{adj}}^2
=1-\frac{SS_E/df_E}{SS_T/df_T}
=1-\frac{MS_E}{SS_T/df_T}.
$$

对同一组响应数据，分母 $SS_T/df_T$ 固定；增加或删除自变量时，误差均方 $MS_E$ 会改变。只有增加变量使 $MS_E$ 降低，调整后的 $R^2$ 才会增大。由于增加变量会使误差自由度减少 1，新变量必须让残差平方和的减少量超过原模型的误差均方[^1]，才能使新模型的调整后 $R^2$ 增大。


## S10.4 逐步回归及其他变量选择方法

教材讨论回归时着重于拟合完整模型。分析设计试验数据时，试验者通常已通过方差分析或效应估计值的正态概率图，对模型形式有较明确的判断。

另一些回归应用来自非预先设计的研究。数据可能是对某个过程日常采集的观测，也可能是从历史记录或资料库取得的存档数据。这类研究往往有中等数量甚至大量候选自变量，分析者希望从中选出适合建立回归模型的子集。离群值、自变量之间的强相关等特点会使问题更加复杂。

选择变量的方法大致包括逐步型方法和全子集回归。逐步型方法每一步向当前模型加入一个变量，或从中删除一个变量。**前向选择**从不含候选变量的模型出发，逐个加入变量，直至得到最终方程；**后向剔除**从包含全部候选变量的模型出发，逐个删除变量。通常所说的**逐步回归**结合前向与后向操作；这些基本方法还有多种变体。

对于 $K$ 个候选变量，全子集回归考察可能的 $2^K$ 个变量子集，并找出有望成为实用模型的方案。$K$ 稍大时，候选模型数便迅速增加；因此已有算法通过隐式搜索，避免逐一显式拟合全部方程。关于变量选择方法的进一步讨论，参见 Montgomery、Peck 和 Vining（2006）或 Myers（1990）的回归教材。

## S10.5 预测响应的方差

教材第 10.5.2 节给出在目标点 $\mathbf x_0'=[1,x_{01},x_{02},\ldots,x_{0k}]$ 处，预测平均响应的方差〔式（10.40）〕：

$$
\operatorname{Var}[\hat y(\mathbf x_0)]
=\sigma^2\mathbf x_0'(\mathbf X'\mathbf X)^{-1}\mathbf x_0.
$$

由式（10.39）的 $\hat y(\mathbf x_0)=\mathbf x_0'\hat{\boldsymbol\beta}$，容易得到

$$
\begin{aligned}
\operatorname{Var}[\hat y(\mathbf x_0)]
&=\operatorname{Var}(\mathbf x_0'\hat{\boldsymbol\beta})\\
&=\mathbf x_0'\operatorname{Var}(\hat{\boldsymbol\beta})\mathbf x_0\\
&=\sigma^2\mathbf x_0'(\mathbf X'\mathbf X)^{-1}\mathbf x_0.
\end{aligned}
$$

Design-Expert 使用教材式（10.41）计算并显示目标点平均响应的置信区间；该结果位于优化菜单的点预测功能中。程序绘制预测标准误的等高线图时，也使用式（10.40）。

## S10.6 预测误差的方差

教材第 10.6 节给出目标点 $\mathbf x_0$ 处未来单次观测的预测区间，其依据是预测误差的方差。未来观测 $y_0$ 的点预测与预测误差分别为

$$
\hat y(\mathbf x_0)=\mathbf x_0'\hat{\boldsymbol\beta},
\qquad e_p=y_0-\hat y(\mathbf x_0).
$$

由于未来观测独立于根据已有数据得到的点预测，

$$
\begin{aligned}
\operatorname{Var}(e_p)
&=\operatorname{Var}(y_0)+\operatorname{Var}[\hat y(\mathbf x_0)]\\
&=\sigma^2+\sigma^2\mathbf x_0'(\mathbf X'\mathbf X)^{-1}\mathbf x_0\\
&=\sigma^2\big[1+\mathbf x_0'(\mathbf X'\mathbf X)^{-1}\mathbf x_0\big].
\end{aligned}
$$

用 $\hat\sigma^2=MS_E$ 代替未知的 $\sigma^2$，再取平方根，就得到教材式（10.42）预测区间所用的标准误。

## S10.7 回归模型中的杠杆值

教材第 10.7.2 节正式定义了设计试验（以及一般回归数据）中每个观测的杠杆值。它是“帽子矩阵”

$$
\mathbf H=\mathbf X(\mathbf X'\mathbf X)^{-1}\mathbf X'
$$

的第 $i$ 个对角元素，即

$$
h_{ii}=\mathbf x_i'(\mathbf X'\mathbf X)^{-1}\mathbf x_i,
$$

其中 $\mathbf x_i'$ 是 $\mathbf X$ 的第 $i$ 行。

杠杆值有两种解释。其一，$h_{ii}$ 衡量设计点离设计空间中心有多远。例如，在编码单位下，$2^k$ 因子设计的所有立方体顶点距中心都为 $\sqrt{k}$[^2]；如果每个点的重复次数相同，它们的杠杆值也相同。其二，杠杆值表示某个观测对拟合模型可能产生的最大影响。接近饱和的设计中，许多甚至所有设计点都可能具有很高的杠杆值。

单个观测的杠杆值最大为 $h_{ii}=1$；若某个设计点有 $n$ 次相同条件下的重复观测，这些观测各自的最大杠杆值为 $1/n$。高杠杆值并不理想：若某个观测的杠杆值等于 1，模型将在该点精确通过它，因此对该点的离群或异常观测十分敏感。重复该设计点可以降低单个观测的杠杆值。

[^1]: 译者注：原文末句称不满足上述条件时“新模型的调整后 $R^2$ 更大”，与前文和公式相反，应为“不会更大”；译文按公式改正。

[^2]: 译者注：原文称 $2^k$ 立方体顶点到中心的距离为 $k$。各编码坐标均为 $\pm1$，欧氏距离应为 $\sqrt{k}$；若指距离的平方，才等于 $k$。译文按欧氏距离修正。
