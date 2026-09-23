# 第 3 章补充材料

## S3.1 因子效应的定义

正如 3.2 节和 3.3 节所指出的，单因子试验的模型有两种写法：**均值模型**（means model）和**效应模型**（effects model）。我们一般使用效应模型

$$
y_{ij} = \mu + \tau_{i} + \varepsilon_{ij} \left\{ \begin{array}{l} i = 1, 2, \dots, a \\ j = 1, 2, \dots, n \end{array} \right.
$$

其中为简单起见，我们处理的是平衡情形（所有因子水平或处理都重复相同的次数）。回忆一下，在写出这个模型时，第 $i$ 个因子水平均值 $\mu_{i}$ 被分解为两个分量，即 $\mu_{i} = \mu + \tau_{i}$，其中 $\tau_{i}$ 是第 $i$ 个**处理效应**（treatment effect），$\mu$ 是总均值。我们通常定义 $\mu = \frac{\displaystyle \sum_{i=1}^{a} \mu_{i}}{a}$，这蕴含着 $\sum_{i=1}^{a} \tau_{i} = 0$。

这实际上是一个任意的定义，定义总“均值”还有其他方式。例如，我们可以定义

$$
\mu = \sum_{i=1}^{a} w_{i} \mu_{i} \quad \text{其中} \quad \sum_{i=1}^{a} w_{i} = 1
$$

这样得到的处理效应满足

$$
\sum_{i=1}^{a} w_{i} \tau_{i} = 0
$$

这里总均值是各个处理均值的加权平均。当各处理中的观测个数不相等时，权重 $w_{i}$ 可以取为该处理样本量的比例 $n_{i}/N$。

## S3.2 期望均方

在 3.3.1 节中，我们推导了单因子**方差分析**（analysis of variance）中误差均方的期望值。我们给出了处理均方期望值的结果，但省略了推导过程。该推导相当直接。

考虑

$$
E\left(MS_{\text{处理}}\right) = E\left(\frac{SS_{\text{处理}}}{a-1}\right)
$$

对于平衡设计

$$
SS_{\text{处理}} = \frac{1}{n} \sum_{i=1}^{a} y_{i.}^{2} - \frac{1}{an} y_{..}^{2}
$$

而模型为

$$
y_{ij} = \mu + \tau_{i} + \varepsilon_{ij} \left\{ \begin{array}{l} i = 1, 2, \dots, a \\ j = 1, 2, \dots, n \end{array} \right.
$$

此外，我们还将用到以下有用的结果：

$$
E(\varepsilon_{ij}) = E(\varepsilon_{i.}) = E(\varepsilon_{..}) = 0, E(\varepsilon_{ij}^{2}) = \sigma^{2}, E(\varepsilon_{i.}^{2}) = n\sigma^{2}, E(\varepsilon_{..}^{2}) = an\sigma^{2}
$$

于是

$$
E\left(SS_{\text{处理}}\right) = E\left(\frac{1}{n} \sum_{i=1}^{a} y_{i.}^{2}\right) - E\left(\frac{1}{an} y_{..}^{2}\right)
$$

考虑上式右端的第一项：

$$
E\left(\frac{1}{n} \sum_{i=1}^{a} y_{i.}^{2}\right) = \frac{1}{n} \sum_{i=1}^{a} E\left(n\mu + n\tau_{i} + \varepsilon_{i.}\right)^{2}
$$

对括号中的表达式取平方再求期望，得到

$$
\begin{array}{c} E\left(\frac{1}{n} \sum_{i=1}^{a} y_{i.}^{2}\right) = \frac{1}{n} \left[ a(n\mu)^{2} + n^{2} \sum_{i=1}^{a} \tau_{i}^{2} + an\sigma^{2} \right] \\ = an\mu^{2} + n \sum_{i=1}^{a} \tau_{i}^{2} + a\sigma^{2} \end{array}
$$

因为三个交叉乘积项均为零。现在考虑 $E(SS_{\text{处理}})$ 右端的第二项：

$$
\begin{array}{r l} E\left(\frac{1}{an} y_{\cdot\cdot}^{2}\right) & = \frac{1}{an} E(an\mu + n \sum_{i=1}^{a} \tau_{i} + \varepsilon_{\cdot\cdot})^{2} \\ & = \frac{1}{an} E(an\mu + \varepsilon_{\cdot\cdot})^{2} \end{array}
$$

因为 $\sum_{i=1}^{a} \tau_{i} = 0$。对括号中的项取平方再求期望，可得

$$
\begin{array}{r l} E\left(\frac{1}{an} y_{\cdot\cdot}^{2}\right) & = \frac{1}{an}[(an\mu)^{2} + an\sigma^{2}] \\ & = an\mu^{2} + \sigma^{2} \end{array}
$$

因为交叉乘积项的期望为零。因此，

$$
\begin{array}{r l} E(SS_{\text{处理}}) & = E(\frac{1}{n} \sum_{i=1}^{a} y_{i.}^{2}) - E(\frac{1}{an} y_{..}^{2}) \\ & = an\mu^{2} + n \sum_{i=1}^{a} \tau_{i}^{2} + a\sigma^{2} - (an\mu^{2} + \sigma^{2}) \\ & = \sigma^{2}(a-1) + n \sum_{i=1}^{a} \tau_{i}^{2} \end{array}
$$

于是处理均方的期望值为

$$
\begin{array}{r l} E(MS_{\text{处理}}) & = E\left(\frac{SS_{\text{处理}}}{a-1}\right) \\ & = \frac{\sigma^{2}(a-1) + n \sum_{i=1}^{a} \tau_{i}^{2}}{a-1} \\ & = \sigma^{2} + \frac{n \sum_{i=1}^{a} \tau_{i}^{2}}{a-1} \end{array}
$$

这就是教材中给出的结果。

## S3.3 $\sigma^{2}$ 的置信区间

在建立方差分析（ANOVA）程序的过程中，我们注意到误差方差 $\sigma^{2}$ 由误差均方来估计，即

$$
\hat{\sigma}^{2} = \frac{SS_{E}}{N-a}
$$

下面我们给出 $\sigma^{2}$ 的**置信区间**（confidence interval）。由于我们假定观测服从正态分布，故

$$
\frac{SS_{E}}{\sigma^{2}}
$$

的分布为 $\chi_{N-a}^{2}$。因此，

$$
P\left(\chi_{1-\alpha/2, N-a}^{2} \leq \frac{SS_{E}}{\sigma^{2}} \leq \chi_{\alpha/2, N-a}^{2}\right) = 1-\alpha
$$

其中 $\chi_{1-\alpha/2, N-a}^{2}$ 和 $\chi_{\alpha/2, N-a}^{2}$ 分别是自由度为 $N-a$ 的 $\chi^{2}$ 分布的（下侧）$\alpha/2$ 分位点和（上侧）$\alpha/2$ 分位点。现在，如果我们对概率式中的表达式作变形，可得

$$
P\left(\frac{SS_{E}}{\chi_{\alpha/2, N-a}^{2}} \leq \sigma^{2} \leq \frac{SS_{E}}{\chi_{1-\alpha/2, N-a}^{2}}\right) = 1-\alpha
$$

因此，误差方差 $\sigma^{2}$ 的 $100(1-\alpha)$ 百分置信区间为

$$
\frac{SS_{E}}{\chi_{\alpha/2, N-a}^{2}} \leq \sigma^{2} \leq \frac{SS_{E}}{\chi_{1-\alpha/2, N-a}^{2}}
$$

这一置信区间表达式也在第 12 章关于随机效应试验的内容中给出。

有时试验者关心误差方差的上界，即 $\sigma^{2}$ 合理地最大能有多大？当已有关于 $\sigma^{2}$ 的先前试验信息，而试验者正在做计算以确定新试验的样本量时，这一点会很有用。$\sigma^{2}$ 的 $100(1-\alpha)$ 百分置信上限为

$$
\sigma^{2} \leq \frac{SS_{E}}{\chi_{1-\alpha, N-a}^{2}}
$$

如果转而希望得到标准差 $\sigma$ 的 $100(1-\alpha)$ 百分置信区间，则

$$
\sigma \leq \sqrt{\frac{SS_{E}}{\chi_{1-\alpha/2, N-a}^{2}}}
$$

## S3.4 处理均值的联合置信区间

在 3.3.3 节中，我们讨论了如何求某个处理均值的置信区间以及一对均值之差的置信区间。我们还说明了如何用 **Bonferroni 方法**（Bonferroni method）求一组处理均值或一组均值之差的**联合置信区间**（simultaneous confidence interval）。本质上，如果要构造一组 $r$ 个置信陈述，Bonferroni 方法只是把 $\alpha/2$ 替换为 $\alpha/(2r)$。这样得到一组 $r$ 个置信区间，其整体置信水平至少为 $100(1-\alpha)$ 个百分点。

为说明它为什么有效，考虑 $r = 2$ 的情形，即我们有两个 $100(1-\alpha)$ 百分置信区间。令 $E_{1}$ 表示第一个置信区间不正确（没有覆盖真均值）这一事件，$E_{2}$ 表示第二个置信区间不正确这一事件。于是

$$
P(E_{1}) = P(E_{2}) = \alpha
$$

两个区间中有一个或两个不正确的概率为

$$
P(E_{1} \cup E_{2}) = P(E_{1}) + P(E_{2}) - P(E_{1} \cap E_{2})
$$

由对立事件的概率，我们可以求出两个区间都正确的概率为

$$
\begin{array}{c} P(\overline{E}_{1} \cap \overline{E}_{2}) = 1 - P(E_{1} \cup E_{2}) \\ = 1 - P(E_{1}) - P(E_{2}) + P(E_{1} \cap E_{2}) \end{array}
$$

现在我们知道 $P(E_{1} \cap E_{2}) \geq 0$，所以由上式中最后一步可得 Bonferroni 不等式

$$
P(\overline{E}_{1} \cap \overline{E}_{2}) \geq 1 - P(E_{1}) - P(E_{2})
$$

在我们的例子中，该不等式的左端就是两个置信区间陈述都正确的概率，而 $P(E_{1}) = P(E_{2}) = \alpha$，所以

$$
\begin{array}{c} P(\overline{E}_{1} \cap \overline{E}_{2}) \geq 1 - \alpha - \alpha \\ \geq 1 - 2\alpha \end{array}
$$

因此，如果我们希望两个置信区间都正确的概率至少为 $1-\alpha$，只要构造 $100(1-\alpha/2)$ 百分的单个置信区间便可保证这一点。

如果有 $r$ 个我们关心的置信区间，可以用数学归纳法证明

$$
\begin{array}{c} P(\overline{E}_{1} \cap \overline{E}_{2} \cap \dots \cap \overline{E}_{r}) \geq 1 - \sum_{i=1}^{r} P(E_{i}) \\ \geq 1 - r\alpha \end{array}
$$

正如教材中所指出的，当希望构造的联合置信区间个数 $r$ 不太大时，Bonferroni 方法效果相当好。随着 $r$ 变大，各个置信区间的长度会增加。单个置信区间的长度可能变得很大，以致这些区间提供的信息不多。此外，各个置信陈述并不需要具有相同的置信水平。可以对一个陈述取 98%，对另一个取 92%，这样得到的两个置信区间的联合置信水平至少为 90%。

## S3.5 定量因子的回归模型

**回归模型**（regression model）将在第 10 章中详细讨论，但它们在本书中相当频繁地出现，因为用方程来表示响应与定量设计变量之间的关系十分方便。当只有一个定量设计因子时，把响应与该因子联系起来的线性回归模型为

$$
y = \beta_{0} + \beta_{1} x + \varepsilon
$$

其中 $x$ 表示设计因子的取值。在单因子试验中有 $N$ 个观测，每个观测都可以按该模型表示如下：

$$
y_{i} = \beta_{0} + \beta_{1} x_{i} + \varepsilon_{i}, i = 1, 2, \dots, N
$$

用**最小二乘**（least squares）方法估计该模型中未知的参数（即各 $\beta$）。这要求选择参数的取值使误差的平方和达到最小。最小二乘函数为

$$
L = \sum_{i=1}^{N} \varepsilon_{i}^{2} = \sum_{i=1}^{N} (y_{i} - \beta_{0} - \beta_{1} x_{i})^{2}
$$

为求最小二乘估计量，我们对 $L$ 关于各 $\beta$ 求偏导数并令其等于零：

$$
\frac{\partial L}{\partial \beta_{0}} = -2 \sum_{i=1}^{N} (y_{i} - \beta_{0} - \beta_{1} x_{i}) = 0
$$

$$
\frac{\partial L}{\partial \beta_{1}} = -2 \sum_{i=1}^{N} (y_{i} - \beta_{0} - \beta_{1} x_{i}) x_{i} = 0
$$

化简后，我们得到**最小二乘正规方程**（least squares normal equations）

$$
N\hat{\beta}_{0} + \hat{\beta}_{1} \sum_{i=1}^{N} x_{i} = \sum_{i=1}^{N} y_{i}
$$

$$
\hat{\beta}_{0} \sum_{i=1}^{N} x_{i} + \hat{\beta}_{1} \sum_{i=1}^{N} x_{i}^{2} = \sum_{i=1}^{N} x_{i} y_{i}
$$

其中 $\hat{\beta}_{0}$ 和 $\hat{\beta}_{1}$ 是模型参数的最小二乘估计量。因此，要用最小二乘方法把这一特定模型拟合到试验数据上，我们只需解正规方程即可。由于只有两个方程、两个未知数，这相当容易。

在教材中，我们对响应变量刻蚀速率（$y$）作为射频功率（$x$）的函数拟合了两个回归模型：上面给出的线性回归模型，以及一个二次模型

$$
y = \beta_{0} + \beta_{1} x + \beta_{2} x^{2} + \varepsilon
$$

二次模型的最小二乘正规方程为

$$
N\hat{\beta}_{0} + \hat{\beta}_{1} \sum_{i=1}^{N} x_{i} + \hat{\beta}_{2} \sum_{i=1}^{N} x_{i}^{2} = \sum_{i=1}^{N} y_{i}
$$

$$
\hat{\beta}_{0} \sum_{i=1}^{N} x_{i} + \hat{\beta}_{1} \sum_{i=1}^{N} x_{i}^{2} + \hat{\beta}_{2} \sum_{i=1}^{N} x_{i}^{3} = \sum_{i=1}^{N} x_{i} y_{i}
$$

$$
\hat{\beta}_{0} \sum_{i=1}^{N} x_{i}^{2} + \hat{\beta}_{1} \sum_{i=1}^{N} x_{i}^{3} + \hat{\beta}_{2} \sum_{i=1}^{N} x_{i}^{4} = \sum_{i=1}^{N} x_{i}^{2} y_{i}
$$

显然，随着模型阶数提高、需要估计的未知参数增多，正规方程会变得更复杂。在第 10 章中我们将用矩阵方法给出一般解法。大多数统计软件包都具有很好的回归模型拟合能力。

## S3.6 关于可估函数的更多内容

在 3.10.1 节中，我们用最小二乘方法估计单因子模型中的参数。假定试验设计是平衡的，我们得到的最小二乘正规方程即式 3-48，重列如下：

$$
\begin{array}{l} an\hat{\mu} + n\hat{\tau}_{1} + n\hat{\tau}_{2} + \dots + n\hat{\tau}_{a} = \sum_{i=1}^{a} \sum_{j=1}^{n} y_{ij} \\ n\hat{\mu} + n\hat{\tau}_{1} = \sum_{j=1}^{n} y_{1j} \\ n\hat{\mu} + n\hat{\tau}_{2} = \sum_{j=1}^{n} y_{2j} \\ \vdots \\ n\hat{\mu} + n\hat{\tau}_{a} = \sum_{j=1}^{n} y_{aj} \end{array}
$$

其中 $an = N$ 是观测总数。正如教材中所指出的，如果把这后 $a$ 个正规方程相加，就得到第一个正规方程。也就是说，这些正规方程不是线性无关的，因而没有唯一解。我们称效应模型是一个**过度参数化**（overparameterized）模型。

解决这一困难的一种做法是向正规方程组再增加一个线性无关的方程。最常见的做法是使用方程 $\sum_{i=1}^{a} \hat{\tau}_{i} = 0$。这与把因子效应定义为对总均值 $\mu$ 的偏离是一致的。如果施加这一约束，正规方程的解为

$$
\begin{array}{l} \hat{\mu} = \overline{y} \\ \hat{\tau}_{i} = \overline{y}_{i} - \overline{y}, i = 1, 2, \dots, a \end{array}
$$

也就是说，总均值由全部 $an$ 个样本观测的平均值来估计，而每个单独的因子效应则由该因子水平的样本均值与全部观测平均值之差来估计。

另一种可能的约束选择是把总均值设为一个常数，比如 $\hat{\mu} = 0$。由此得到的解为

$$
\begin{array}{l} \hat{\mu} = 0 \\ \hat{\tau}_{i} = \overline{y}_{i}, i = 1, 2, \dots, a \end{array}
$$

还有第三种选择是 $\hat{\tau}_{a} = 0$。例如，SAS 软件采用的就是这种做法。这一约束选择给出的解为

$$
\begin{array}{l} \hat{\mu} = \overline{y}_{a} \\ \hat{\tau}_{i} = \overline{y}_{i} - \overline{y}_{a}, i = 1, 2, \dots, a-1 \\ \hat{\tau}_{a} = 0 \end{array}
$$

可以用来求解正规方程的约束有无穷多种。幸运的是，正如书中所指出的，这其实无关紧要。对上面的三个解（事实上对正规方程组的任何解）我们都有

$$
\hat{\mu}_{i} = \hat{\mu} + \hat{\tau}_{i} = \overline{y}_{i}, i = 1, 2, \dots, a
$$

也就是说，第 $i$ 个因子水平均值的最小二乘估计量总是该因子水平下观测的样本均值。因此，即使我们无法得到效应模型中各参数的唯一估计，我们仍然可以对这些参数中我们所关心的某个函数得到唯一的估计量。

这就是**可估函数**（estimable function）的思想。模型参数的任意函数，只要不论选择何种约束来求解正规方程都能被唯一地估计，就是一个可估函数。

哪些函数是可估的？可以证明，任意观测的期望值都是可估的。由于

$$
E(y_{ij}) = \mu + \tau_{i}
$$

所以如上所示，第 $i$ 个处理的均值是可估的。任何作为正规方程左端线性组合的函数也是可估的。例如，用第二个正规方程减去第三个正规方程，得到 $\tau_{2} - \tau_{1}$。因此，任意两个处理效应之差是可估的。一般地，处理效应的任意**对照**（contrast） $\sum_{i=1}^{a} c_{i} \tau_{i}$（其中 $\sum_{i=1}^{a} c_{i} = 0$）都是可估的。注意各个模型参数 $\mu, \tau_{1}, \cdots, \tau_{a}$ 本身并不是可估的，因为不存在能分别产生这些参数的正规方程的线性组合。不过，这通常不成问题，因为如前所述，可估函数对应于试验者所关心的模型参数函数。

关于可估函数有一个出色且非常易读的论述，参见 Myers, R. H. and Milton, J. S. (1991), A First Course in the Theory of the Linear Model, PWS-Kent, Boston. MA.

## S3.7 回归与方差分析之间的关系

3.10 节探讨了方差分析（ANOVA）模型与回归模型之间的一些联系。我们说明了如何用最小二乘方法来估计模型参数，以及如何通过一种基于回归的程序——称为**一般回归显著性检验**（general regression significance test）——来导出方差分析，从而得到 ANOVA 检验统计量。每个 ANOVA 模型都可以显式地写成一个等价的线性回归模型。下面我们对有 $a = 3$ 个处理的单因子试验说明这一做法。

单因子平衡 ANOVA 模型为

$$
y_{ij} = \mu + \tau_{i} + \varepsilon_{ij} \left\{ \begin{array}{c} i = 1, 2, 3 \\ j = 1, 2, \dots, n \end{array} \right.
$$

等价的回归模型为

$$
y_{ij} = \beta_{0} + \beta_{1} x_{1j} + \beta_{2} x_{2j} + \varepsilon_{ij} \left\{ \begin{array}{l} i = 1, 2, 3 \\ j = 1, 2, \dots, n \end{array} \right.
$$

其中变量 $x_{1j}$ 和 $x_{2j}$ 定义如下：

$$
\begin{array}{l} x_{1j} = \left\{ \begin{array}{c} 1 \text{ 若观测 } j \text{ 来自处理 } 1 \\ 0 \text{ 其他} \end{array} \right. \\ x_{2j} = \left\{ \begin{array}{c} 1 \text{ 若观测 } j \text{ 来自处理 } 2 \\ 0 \text{ 其他} \end{array} \right. \end{array}
$$

回归模型中参数与 ANOVA 模型中参数之间的关系很容易确定。例如，如果观测来自处理 1，则 $x_{1j} = 1$ 且 $x_{2j} = 0$，回归模型为

$$
\begin{array}{c} y_{1j} = \beta_{0} + \beta_{1}(1) + \beta_{2}(0) + \varepsilon_{1j} \\ = \beta_{0} + \beta_{1} + \varepsilon_{1j} \end{array}
$$

由于在 ANOVA 模型中这些观测由 $y_{1j} = \mu + \tau_{1} + \varepsilon_{1j}$ 定义，这意味着

$$
\beta_{0} + \beta_{1} = \mu_{1} = \mu + \tau_{1}
$$

类似地，如果观测来自处理 2，则

$$
\begin{array}{c} y_{2j} = \beta_{0} + \beta_{1}(0) + \beta_{2}(1) + \varepsilon_{2j} \\ = \beta_{0} + \beta_{2} + \varepsilon_{2j} \end{array}
$$

参数之间的关系为

$$
\beta_{0} + \beta_{2} = \mu_{2} = \mu + \tau_{2}
$$

最后，考虑来自处理 3 的观测，其回归模型为

$$
\begin{array}{c} y_{3j} = \beta_{0} + \beta_{1}(0) + \beta_{2}(0) + \varepsilon_{3j} \\ = \beta_{0} + \varepsilon_{3j} \end{array}
$$

于是有

$$
\beta_{0} = \mu_{3} = \mu + \tau_{3}
$$

因此，在单因素方差分析模型的回归模型表述中，回归系数描述的是前两个处理均值与第三个处理均值的比较，即

$$
\begin{array}{l} \beta_{0} = \mu_{3} \\ \beta_{1} = \mu_{1} - \mu_{3} \\ \beta_{2} = \mu_{2} - \mu_{3} \end{array}
$$

一般地，如果有 $a$ 个处理，回归模型将有 $a-1$ 个回归自变量，即

$$
y_{ij} = \beta_{0} + \beta_{1} x_{1j} + \beta_{2} x_{2j} + \dots + \beta_{a-1} x_{a-1} + \varepsilon_{ij} \left\{ \begin{array}{l} i = 1, 2, \dots, a \\ j = 1, 2, \dots, n \end{array} \right.
$$

其中

$$
x_{ij} = \left\{ \begin{array}{c} 1 \text{ 若观测 } j \text{ 来自处理 } i \\ 0 \text{ 其他} \end{array} \right.
$$

由于这些回归自变量只取 0 和 1 两个值，它们常被称为**指示变量**（indicator variable）。ANOVA 模型参数与回归模型参数之间的关系为

$$
\begin{array}{l} \beta_{0} = \mu_{a} \\ \beta_{i} = \mu_{i} - \mu_{a}, i = 1, 2, \dots, a-1 \end{array}
$$

因此截距总是第 $a$ 个处理的均值，而回归系数 $\beta_{i}$ 估计的是第 $i$ 个处理均值与第 $a$ 个处理均值之差。

现在考虑假设检验。假设我们希望检验所有处理均值都相等（即通常的原假设）。如果该原假设为真，则回归模型中的参数变为

$$
\begin{array}{l} \beta_{0} = \mu_{a} \\ \beta_{i} = 0, i = 1, 2, \dots, a-1 \end{array}
$$

利用一般回归显著性检验程序，我们可以针对这一假设建立检验。它将与一元方差分析中的 $F$ 统计量检验完全相同。

大多数回归软件包会自动检验所有模型回归系数（截距除外）均为零这一假设。我们将用 Minitab 并结合例 3.1 中等离子体刻蚀试验的数据来说明这一点。回忆一下，在该例中工程师对确定 RF 功率对刻蚀速率的影响感兴趣，他实施了一个含有 4 个 RF 功率水平、每个水平 5 次重复的完全随机化试验。为方便起见，这里重复表 3.1 中的数据：

**表 3.1 等离子体刻蚀试验的刻蚀速率数据（单位：Å/min）**

| 功率 (W) | 观测值 1 | 2 | 3 | 4 | 5 |
| --- | --- | --- | --- | --- | --- |
| 160 | 575 | 542 | 530 | 539 | 570 |
| 180 | 565 | 593 | 590 | 579 | 610 |
| 200 | 600 | 651 | 610 | 637 | 629 |
| 220 | 725 | 700 | 715 | 685 | 710 |

按上述方式把数据转换成 $x_{ij}$ 的 0/1 指示变量。由于有 4 个处理，所以只有 3 个 $x$。用作 Minitab 输入的编码数据如下所示：

```txt
x1 x2 x3 Etch rate
1 0 0 575
1 0 0 542
1 0 0 530
1 0 0 539
1 0 0 570
0 1 0 565
0 1 0 593
0 1 0 590
0 1 0 579
0 1 0 610
0 0 1 600
0 0 1 651
0 0 1 610
0 0 1 637
0 0 1 629
0 0 0 725
0 0 0 700
0 0 0 715
0 0 0 685
```

使用上述数据表运行了 Minitab 的回归模块，其中 $x_{1}$ 至 $x_{3}$ 用作预测变量，变量“Etch rate”用作响应。输出如下所示。

Regression Analysis: Etch rate versus x1, x2, x3

```txt
The regression equation is
Etch rate = 707 - 156 x1 - 120 x2 - 81.6 x3

Predictor Coef SE Coef T P
Constant 707.000 8.169 86.54 0.000
x1 -155.80 11.55 -13.49 0.000
x2 -119.60 11.55 -10.35 0.000
x3 -81.60 11.55 -7.06 0.000

S = 18.2675 R-Sq = 92.6% R-Sq(adj) = 91.2%

Analysis of Variance

Source DF SS MS F P
Regression 3 66871 22290 66.80 0.000
Residual Error 16 5339 334
```

注意，这一回归输出中的方差分析表与表 3.4 中的方差分析显示完全相同（除舍入误差外）。因此，检验该回归模型中回归系数 $\beta_{1} = \beta_{2} = \beta_{3} = \beta_{4} = 0$ 这一假设，等价于检验原 ANOVA 模型表述中处理均值相等的原假设。

还要注意，上表中截距（即“常数”项）的估计值就是第 $4$ 个处理的均值。此外，每个回归系数都恰好是某个处理均值与第 $4$ 个处理均值之差。
