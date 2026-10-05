# 第 15 章补充材料

## S15.1 变换的形式

正文第 3.4.3 节介绍了用变换稳定响应方差，并在方差不齐与非正态同时发生时改善正态性。第 15.1.1 节的 Box–Cox 方法给出了解析选择变换的方法。不过，许多试验者仍通过尝试第 3 章表 3.9 或软件菜单中的简单幂变换，如 $\sqrt y,\ln y,1/y$，作经验选择。

幂变换可从理论上解释。设响应 y 的均值为 $\mu$，方差为 $\sigma^2=f(\mu)$，希望寻找 $x=h(y)$，使其方差为不依赖均值的常数。围绕 $\mu$ 作一阶 Taylor 展开：

$$
x=h(y)=h(\mu)+h'(\mu)(y-\mu)+R\approx h(\mu)+h'(\mu)(y-\mu).
$$

忽略余项 R 后，

$$
E(x)\approx h(\mu),\qquad
V(x)\approx E\{h'(\mu)(y-\mu)\}^2=\sigma^2[h'(\mu)]^2=f(\mu)[h'(\mu)]^2.
$$

令方差为 $c^2$，便有

$$
h'(\mu)=\frac c{\sqrt{f(\mu)}},\qquad
h(\mu)=c\int\frac{dt}{\sqrt{f(t)}}=cG(\mu)+k,
$$

其中 k 是积分常数。若均值等于方差，如 Poisson 分布，则 $f(t)=t$，于是

$$
h(\mu)=c\int t^{-1/2}\,dt=2c\sqrt\mu+k.
$$

因此，平方根变换可近似稳定 Poisson 数据或均值与方差接近的计数数据的方差，与正文建议一致。

再设标准差约与均值成正比，即 $\sigma^2=\mu^2$，则 $f(t)=t^2$，有

$$
h(\mu)=c\int\frac{dt}{t}=c\ln\mu+k,\qquad\mu>0.
$$

因此，对这类正响应，对数变换是适当的方差稳定化变换。[^1]

## S15.2 在 Box–Cox 方法中选择 λ

正文第 15.1.1 节指出 Box–Cox 方法以极大似然为依据。对归一化的变换响应，略去与 λ 无关的常数后，剖面对数似然为

$$
\ell(\lambda)=-\frac n2\ln SS_E(\lambda).
$$

最大化它等价于最小化误差平方和。近似 $100(1-\alpha)\%$ 置信区间由满足下式的 λ 组成：

$$
\ell(\hat\lambda)-\ell(\lambda)\le\frac12\chi_{\alpha,1}^2.
$$

n 为样本量，$\chi_{\alpha,1}^2$ 为自由度 1 的卡方分布上侧 α 分位数。在对数似然图中，于 $\ell(\hat\lambda)-\chi_{\alpha,1}^2/2$ 画水平线，与曲线两个交点确定区间端点。若绘制误差平方和图，对应高度为

$$
SS^*=SS_E(\hat\lambda)\exp(\chi_{\alpha,1}^2/n).
$$

正文式（15.2）改用 $1+t_{\alpha/2,\nu}^2/\nu$，其中 ν 为误差自由度。有的作者使用 $1+z_{\alpha/2}^2/\nu$，或将分母与 t 自由度换成 n。这些近似基于 $e^x\approx1+x$，以及 $\chi_1^2$ 与标准正态平方的关系、自由度较大时 t 与正态分布接近。n 与 ν 的选择有时会影响区间，通常差异较小。[^2]

## S15.3 广义线性模型

正文第 15.1.2 节讨论了正态性和恒定方差不成立时，替代响应变换的 GLM 方法；例 15.2～15.4 展示了它在设计试验中的应用。

GLM 将非正态响应和均值与线性预测子的非线性关系统一起来。响应分布属于指数族，包括正态、Poisson、二项、指数和 gamma 分布；正态理论线性模型是其特例。

先讨论 Logistic 回归，其响应只有两个可能结果，通常称为成功与失败，以 1、0 表示，名称本身是任意的。再考虑计数响应，例如例 15.3 的产品缺陷数，或一年中登陆美国的大西洋飓风数等较少发生事件的计数，最后说明 GLM 如何统一这些情形。

### S15.3.1 二元响应变量模型

响应只有 0、1 两个值，可以来自定性结果。例如，半导体器件功能电测中，成功表示正常工作，失败可能由短路、断路或其他问题造成。

先尝试模型

$$
y_i=\boldsymbol x_i'\boldsymbol\beta+\epsilon_i,
\qquad\boldsymbol x_i'=[1,x_{i1},\ldots,x_{ik}],\quad
\boldsymbol\beta'=[\beta_0,\ldots,\beta_k].
$$

$\boldsymbol x_i'\boldsymbol\beta$ 为线性预测子。假定 $y_i$ 为 Bernoulli 变量，$P(y_i=1)=\pi_i$，$P(y_i=0)=1-\pi_i$。若 $E(\epsilon_i)=0$，则

$$
E(y_i)=1\pi_i+0(1-\pi_i)=\pi_i=\boldsymbol x_i'\boldsymbol\beta.
$$

这意味着响应函数就是成功概率，但该模型有实质问题。首先，误差只能取

$$
\epsilon_i=\begin{cases}1-\boldsymbol x_i'\boldsymbol\beta,&y_i=1,\\-\boldsymbol x_i'\boldsymbol\beta,&y_i=0,\end{cases}
$$

因此不可能正态。其次，方差不恒定：

$$
V(y_i)=(1-\pi_i)^2\pi_i+\pi_i^2(1-\pi_i)=\pi_i(1-\pi_i)=E(y_i)[1-E(y_i)].
$$

误差方差与响应方差相同，均为均值的函数。最后，必须满足 $0\le E(y_i)=\pi_i\le1$，而线性响应函数不能一般地保证这一限制。

二元响应通常采用单调的 S 形或反 S 形函数，即 Logistic 响应函数：

$$
E(y)=\pi=\frac{\exp(\boldsymbol x'\boldsymbol\beta)}{1+\exp(\boldsymbol x'\boldsymbol\beta)}
=\frac1{1+\exp(-\boldsymbol x'\boldsymbol\beta)}.
$$

令

$$
\eta=\ln\left(\frac\pi{1-\pi}\right),
$$

即可把成功概率的函数表示为线性预测子。

因此 $\eta=\boldsymbol x'\boldsymbol\beta$。这称为概率的 logit 变换；$\pi/(1-\pi)$ 称为**优势**，其对数称为对数优势。其他类似形状的函数也可采用。例如，probit 链接使用标准正态分布函数的逆函数 $\Phi^{-1}(\pi)$；互补双对数链接为 $\ln[-\ln(1-\pi)]$，对应响应函数不关于 $\pi=0.5$ 对称。probit 同样可以包含多个预测变量。[^3]

### S15.3.2 Logistic 回归模型的参数估计

模型的一般形式为 $y_i=E(y_i)+\epsilon_i$，其中各观测为独立 Bernoulli 变量，

$$
E(y_i)=\pi_i=\frac{e^{\boldsymbol x_i'\boldsymbol\beta}}{1+e^{\boldsymbol x_i'\boldsymbol\beta}}.
$$

用极大似然估计线性预测子的参数。各观测的概率函数、似然及对数似然为

$$
\begin{aligned}
f_i(y_i)&=\pi_i^{y_i}(1-\pi_i)^{1-y_i},\qquad y_i=0,1,\\
L(\boldsymbol\beta)&=\prod_{i=1}^nf_i(y_i),\\
\ell(\boldsymbol\beta)&=\sum_i y_i\ln\frac{\pi_i}{1-\pi_i}+\sum_i\ln(1-\pi_i)\\
&=\sum_i y_i\boldsymbol x_i'\boldsymbol\beta-\sum_i\ln[1+e^{\boldsymbol x_i'\boldsymbol\beta}].
\end{aligned}
$$

设计试验中，各 x 水平下常有重复试验。令 $y_i$ 为第 i 组成功数，$n_i$ 为该组试验数，则忽略与参数无关的二项系数后，

$$
\ell(\boldsymbol\beta)=\sum_i[y_i\ln\pi_i+(n_i-y_i)\ln(1-\pi_i)]
=\sum_i[y_i\boldsymbol x_i'\boldsymbol\beta-n_i\ln(1+e^{\boldsymbol x_i'\boldsymbol\beta})].
$$

可用数值搜索求极大似然估计，也可使用**迭代重加权最小二乘（IRLS）**。极大似然估计满足得分方程 $\partial\ell/\partial\boldsymbol\beta=0$。由于

$$
\frac{\partial\ell}{\partial\pi_i}=\frac{y_i}{\pi_i}-\frac{n_i-y_i}{1-\pi_i},\qquad
\frac{\partial\pi_i}{\partial\boldsymbol\beta}=\pi_i(1-\pi_i)\boldsymbol x_i,
$$

由链式法则得到

$$
\frac{\partial\ell}{\partial\boldsymbol\beta}=\sum_i(y_i-n_i\pi_i)\boldsymbol x_i=\boldsymbol X'(\boldsymbol y-\boldsymbol\mu)=0,
\quad\boldsymbol\mu'=[n_1\pi_1,\ldots,n_n\pi_n].
$$

它与普通最小二乘的正规方程形式相似：线性回归中 $\boldsymbol\mu=\boldsymbol X\boldsymbol\beta$，$\boldsymbol X'\boldsymbol X\hat{\boldsymbol\beta}=\boldsymbol X'\boldsymbol y$，也可写成 $\boldsymbol X'(\boldsymbol y-\boldsymbol\mu)=0$。

Newton–Raphson 方法用局部一阶 Taylor 展开求解。令 $p_i=y_i/n_i$，$\eta_i=\boldsymbol x_i'\boldsymbol\beta$，在当前值附近近似为

$$
p_i-\pi_i\approx\frac{\partial\pi_i}{\partial\eta_i}(\eta_i^*-\eta_i),
\qquad\frac{\partial\eta_i}{\partial\boldsymbol\beta}=\boldsymbol x_i,
\qquad\frac{\partial\pi_i}{\partial\eta_i}=\pi_i(1-\pi_i).
$$

于是

$$
y_i-n_i\pi_i=n_i(p_i-\pi_i)\approx n_i\pi_i(1-\pi_i)(\eta_i^*-\eta_i).
$$

由比例响应的方差及一阶近似，工作响应在 η 尺度的方差为 $1/[n_i\pi_i(1-\pi_i)]$。将其组成对角矩阵 $\boldsymbol V$，则线性化得分方程为

$$
\boldsymbol X'\boldsymbol V^{-1}(\boldsymbol\eta^*-\boldsymbol X\boldsymbol\beta)=0,
$$

形式上的解为 $(\boldsymbol X'\boldsymbol V^{-1}\boldsymbol X)^{-1}\boldsymbol X'\boldsymbol V^{-1}\boldsymbol\eta^*$。但 $\boldsymbol\eta^*$ 未知，故用线性化工作响应

$$
z_i=\eta_i+(p_i-\pi_i)\frac{\partial\eta_i}{\partial\pi_i}
=\eta_i+\frac{p_i-\pi_i}{\pi_i(1-\pi_i)}
$$

替代。它的随机部分方差为

$$
\frac{\pi_i(1-\pi_i)}{n_i}\left[\frac1{\pi_i(1-\pi_i)}\right]^2
=\frac1{n_i\pi_i(1-\pi_i)}.
$$

因此更新公式为

$$
\hat{\boldsymbol\beta}_{\mathrm{new}}=(\boldsymbol X'\boldsymbol V^{-1}\boldsymbol X)^{-1}\boldsymbol X'\boldsymbol V^{-1}\boldsymbol z.
$$

IRLS 的步骤为：

1. 取得初始估计 $\hat{\boldsymbol\beta}_0$，原材料建议用普通最小二乘初始化。
2. 根据当前估计求 $\boldsymbol\pi$ 和 $\boldsymbol V$。
3. 求 $\boldsymbol\eta=\boldsymbol X\hat{\boldsymbol\beta}$。
4. 根据当前 η、π 构造工作响应 z。
5. 加权最小二乘更新 β，重复直至满足收敛条件。

在模型正确且通常正则条件成立时，最终估计渐近满足

$$
E(\hat{\boldsymbol\beta})\approx\boldsymbol\beta,\qquad
V(\hat{\boldsymbol\beta})\approx(\boldsymbol X'\boldsymbol V^{-1}\boldsymbol X)^{-1}.
$$

拟合概率为

$$
\hat\pi_i=\frac{e^{\boldsymbol x_i'\hat{\boldsymbol\beta}}}{1+e^{\boldsymbol x_i'\hat{\boldsymbol\beta}}}
=\frac1{1+e^{-\boldsymbol x_i'\hat{\boldsymbol\beta}}}.
$$

### S15.3.3 Logistic 回归模型的参数解释

先考虑一个预测变量：$\hat\eta(x_i)=\hat\beta_0+\hat\beta_1x_i$。在 $x_i+1$ 处，$\hat\eta(x_i+1)=\hat\beta_0+\hat\beta_1(x_i+1)$，两者之差为 $\hat\beta_1$，即

$$
\ln\frac{\operatorname{odds}_{x_i+1}}{\operatorname{odds}_{x_i}}=\hat\beta_1,
\qquad\widehat{OR}=\frac{\operatorname{odds}_{x_i+1}}{\operatorname{odds}_{x_i}}=e^{\hat\beta_1}.
$$

它表示预测变量增加一个单位时成功**优势**的乘法变化，不是成功概率增加的倍数。变化 d 个单位时，优势比为 $e^{d\hat\beta_1}$。多元 Logistic 回归中的解释相同：其他变量固定且不涉及随之变化的交互作用项时，$e^{\hat\beta_j}$ 是 $x_j$ 增加一个单位的优势比。[^4]

### S15.3.4 模型参数的假设检验

GLM 检验常基于似然比，是依赖渐近理论的大样本方法，产生称为**离差**的统计量。

**模型离差。** 比较拟合模型与饱和模型的对数似然。饱和模型有足够参数完全拟合各观测组。对未分组二元数据，令 $\hat\pi_i=y_i$ 可使饱和似然达到 1，对数似然为零。拟合模型参数较少，其对数似然不能超过饱和模型。定义

$$
D(\boldsymbol\beta)=2[\ell(\text{饱和模型})-\ell(\hat{\boldsymbol\beta})].
$$

在适当的大样本条件下，如分组二项数据的各组有足够信息，离差近似服从自由度 n−p 的卡方分布。大离差表明模型不适合，小离差表明拟合接近饱和模型；超过 $\chi_{\alpha,n-p}^2$ 时判为不适合。对每组只有一个 Bernoulli 观测的数据，不能仅凭总样本量大就直接采用这一拟合优度近似。[^5] 正态线性回归的离差是残差平方和除以 $\sigma^2$。

**参数子集检验。** 将线性预测子分为

$$
\boldsymbol\eta=\boldsymbol X\boldsymbol\beta=\boldsymbol X_1\boldsymbol\beta_1+\boldsymbol X_2\boldsymbol\beta_2,
$$

完整模型有 p 个参数，$\boldsymbol\beta_2$ 有 r 个。检验 $H_0:\boldsymbol\beta_2=0$ 时，简化模型为 $\boldsymbol\eta=\boldsymbol X_1\boldsymbol\beta_1$。简化模型离差不小于完整模型；若增加不多，不能拒绝原假设，若增量大，则至少一个被删参数可能非零。部分离差为

$$
D(\boldsymbol\beta_2\mid\boldsymbol\beta_1)=D(\boldsymbol\beta_1)-D(\boldsymbol\beta)
=2[\ell(\hat{\boldsymbol\beta})-\ell(\hat{\boldsymbol\beta}_1)].
$$

自由度差为 $n-(p-r)-(n-p)=r$。在原假设及正则条件下，大样本近似为 $\chi_r^2$；达到或超过 $\chi_{\alpha,r}^2$ 时拒绝原假设。它就是似然比统计量

$$
-2\ln\frac{L(\hat{\boldsymbol\beta}_1)}{L(\hat{\boldsymbol\beta})},
$$

因为相减时饱和模型对数似然抵消。

**单个系数检验。** 检验 $H_0:\beta_j=0$、$H_1:\beta_j\ne0$，可用上述离差差，也可采用 Wald 推断。极大似然估计在大样本下近似正态且偏差较小，其协方差可由对数似然二阶偏导求得。令 Hessian 矩阵

$$
G_{ij}=\frac{\partial^2\ell(\boldsymbol\beta)}{\partial\beta_i\partial\beta_j},\qquad i,j=0,\ldots,k.
$$

在估计处计算，有

$$
\hat{\boldsymbol\Sigma}=V(\hat{\boldsymbol\beta})\approx-\boldsymbol G(\hat{\boldsymbol\beta})^{-1}.
$$

对角元素平方根为标准误，统计量为

$$
Z_0=\frac{\hat\beta_j}{\operatorname{se}(\hat\beta_j)}.
$$

参考分布为标准正态。有些软件将其平方，与自由度 1 的卡方分布比较。也可构造 $\hat\beta_j\pm z_{\alpha/2}\operatorname{se}(\hat\beta_j)$ 的 Wald 置信区间。

### S15.3.5 Poisson 回归

另一种非正态响应是较少发生事件的计数，如产品缺陷、软件错误或环境粒子数。希望建立计数与预测变量的关系，例如产品缺陷数与生产条件的关系。响应取 $0,1,\ldots$，常采用 Poisson 概率模型

$$
f(y)=\frac{e^{-\mu}\mu^y}{y!},\qquad y=0,1,\ldots,\quad\mu>0,
$$

其均值与方差均为 μ。Poisson 回归写为 $y_i=E(y_i)+\epsilon_i$。

设 $E(y_i)=\mu_i$，链接函数满足

$$
g(\mu_i)=\beta_0+\beta_1x_{i1}+\cdots+\beta_kx_{ik}=\boldsymbol x_i'\boldsymbol\beta,
\qquad\mu_i=g^{-1}(\boldsymbol x_i'\boldsymbol\beta).
$$

常用的恒等链接为 $g(\mu_i)=\mu_i$，此时均值就是线性预测子。对数链接为 $g(\mu_i)=\ln\mu_i$，此时 $\mu_i=\exp(\boldsymbol x_i'\boldsymbol\beta)$，它保证预测均值为正，特别适用于 Poisson 回归。

独立观测的似然和对数似然为

$$
L(\boldsymbol\beta)=\prod_i\frac{e^{-\mu_i}\mu_i^{y_i}}{y_i!}
=\frac{\prod_i\mu_i^{y_i}\exp(-\sum_i\mu_i)}{\prod_iy_i!},
$$

$$
\ell(\boldsymbol\beta)=\sum_i y_i\ln\mu_i-\sum_i\mu_i-\sum_i\ln(y_i!).
$$

指定链接后，可像 Logistic 回归那样用 IRLS 求极大似然估计。拟合值为 $\hat y_i=g^{-1}(\boldsymbol x_i'\hat{\boldsymbol\beta})$；恒等链接给出 $\boldsymbol x_i'\hat{\boldsymbol\beta}$，对数链接给出其指数。推断方法也类似：离差用于总体拟合评价，离差差用于参数子集的似然比检验，Wald 方法用于单个参数的检验与置信区间。

### S15.3.6 广义线性模型

前述模型都属于 GLM 家族，将通常的正态线性回归与 Logistic、Poisson 等模型统一起来。关键假定是响应分布属于指数族，包括正态、二项、Poisson、逆高斯、指数和 gamma 分布。其一般形式为

$$
f(y_i;\theta_i,\phi)=\exp\left\{\frac{y_i\theta_i-b(\theta_i)}{a(\phi)}+h(y_i,\phi)\right\}.
$$

φ 是尺度参数，$\theta_i$ 是自然参数。指数族满足

$$
\mu_i=E(y_i)=b'(\theta_i),\qquad
V(y_i)=b''(\theta_i)a(\phi)=\frac{d\mu_i}{d\theta_i}a(\phi).
$$

令方差函数 $v(\mu_i)=V(y_i)/a(\phi)=d\mu_i/d\theta_i$，则 $d\theta_i/d\mu_i=1/v(\mu_i)$。

**正态分布。** 将密度重写为

$$
\begin{aligned}
f(y)&=\frac1{\sqrt{2\pi\sigma^2}}\exp\left[-\frac{(y-\mu)^2}{2\sigma^2}\right]\\
&=\exp\left[\frac{y\mu-\mu^2/2}{\sigma^2}-\frac{y^2}{2\sigma^2}-\frac12\ln(2\pi\sigma^2)\right].
\end{aligned}
$$

因此 $\theta=\mu$、$b(\theta)=\theta^2/2$、$a(\phi)=\sigma^2$、$h(y,\phi)=-y^2/(2\sigma^2)-\ln(2\pi\sigma^2)/2$，由导数得到 $E(y)=\mu$、$V(y)=\sigma^2$。

**二项分布。** 对 $y\sim\operatorname{Bin}(m,\pi)$，

$$
\begin{aligned}
f(y)&=\binom my\pi^y(1-\pi)^{m-y}\\
&=\exp\left[y\ln\frac\pi{1-\pi}+m\ln(1-\pi)+\ln\binom my\right].
\end{aligned}
$$

于是 $\theta=\ln[\pi/(1-\pi)]$、$\pi=e^\theta/(1+e^\theta)$、$b(\theta)=m\ln(1+e^\theta)=-m\ln(1-\pi)$、$a(\phi)=1$、$h(y,\phi)=\ln\binom my$。利用 $d\pi/d\theta=\pi(1-\pi)$，得 $E(y)=m\pi$ 和 $V(y)=m\pi(1-\pi)$。

**Poisson 分布。**

$$
f(y)=\frac{\lambda^ye^{-\lambda}}{y!}=\exp[y\ln\lambda-\lambda-\ln(y!)].
$$

所以 $\theta=\ln\lambda$、$b(\theta)=e^\theta=\lambda$、$a(\phi)=1$、$h(y,\phi)=-\ln(y!)$。因 $d\lambda/d\theta=\lambda$，均值与方差均为 λ。

### S15.3.7 链接函数与线性预测子

GLM 为平均响应的适当函数建立线性模型：

$$
\eta_i=g[E(y_i)]=g(\mu_i)=\boldsymbol x_i'\boldsymbol\beta,
\qquad E(y_i)=g^{-1}(\eta_i).
$$

若选择 $\eta_i=\theta_i$，称为**典则链接**。常见链接见表 S15.1。

表 S15.1 广义线性模型的典则链接

| 分布 | 链接 |
| --- | --- |
| 正态 | $\eta_i=\mu_i$，恒等链接 |
| 二项 | $\eta_i=\ln[\pi_i/(1-\pi_i)]$，logit 链接 |
| Poisson | $\eta_i=\ln\mu_i$，对数链接 |
| 指数 | $\eta_i=-1/\mu_i$，负倒数链接 |
| gamma | $\eta_i=-1/\mu_i$，负倒数链接（采用相应尺度参数化） |

原材料的倒数链接写为正号。常用软件也可能采用 $1/\mu$ 的符号约定；但若严格要求此处自然参数等于 η，应使用与密度参数化一致的负号。[^6]

其他选择包括 probit 链接 $\eta_i=\Phi^{-1}[E(y_i)]$（用于概率响应）、互补双对数链接 $\eta_i=\ln\{-\ln[1-E(y_i)]\}$，以及幂链接

$$
\eta_i=\begin{cases}[E(y_i)]^\lambda,&\lambda\ne0,\\\ln[E(y_i)],&\lambda=0.\end{cases}
$$

GLM 有两个基本组成部分：响应分布与链接函数。选择链接类似于选择变换，但链接作用于均值，并利用响应的自然分布，而不是直接变换观测值。不合适的链接也会造成明显拟合问题。

### S15.3.8 广义线性模型的参数估计

理论依据是极大似然，实际通常通过 IRLS 实现。先考虑典则链接，对数似然为

$$
\ell(\boldsymbol\beta)=\sum_i\left[\frac{y_i\theta_i-b(\theta_i)}{a(\phi)}+h(y_i,\phi)\right].
$$

因 $\theta_i=\eta_i=\boldsymbol x_i'\boldsymbol\beta$，

$$
\frac{\partial\ell}{\partial\boldsymbol\beta}
=\frac1{a(\phi)}\sum_i[y_i-b'(\theta_i)]\boldsymbol x_i
=\frac1{a(\phi)}\boldsymbol X'(\boldsymbol y-\boldsymbol\mu).
$$

这给出 p=k+1 个得分方程 $\boldsymbol X'(\boldsymbol y-\boldsymbol\mu)=0$，其中 μ 有 n 个元素。Logistic 情形中的元素为 $n_i\pi_i$。

一阶 Taylor 近似为

$$
y_i-\mu_i\approx\frac{d\mu_i}{d\eta_i}(\eta_i^*-\eta_i)
=v(\mu_i)(\eta_i^*-\eta_i).
$$

线性化工作响应的方差近似为

$$
\left(\frac{d\eta_i}{d\mu_i}\right)^2V(y_i)
=\frac{a(\phi)}{v(\mu_i)}.
$$

令 V 的对角元素为 $1/v(\mu_i)$（暂不含共同尺度 $a(\phi)$），则 $\boldsymbol y-\boldsymbol\mu\approx\boldsymbol V^{-1}(\boldsymbol\eta^*-\boldsymbol\eta)$，得分方程变成

$$
\boldsymbol X'\boldsymbol V^{-1}(\boldsymbol\eta^*-\boldsymbol X\boldsymbol\beta)=0.
$$

未知 η* 用工作响应替代：

$$
z_i=\hat\eta_i+(y_i-\hat\mu_i)\frac{d\eta_i}{d\mu_i},\qquad
\hat{\boldsymbol\beta}_{\mathrm{new}}=(\boldsymbol X'\boldsymbol V^{-1}\boldsymbol X)^{-1}\boldsymbol X'\boldsymbol V^{-1}\boldsymbol z.
$$

例如，对二项比例，$d\eta_i/d\pi_i=1/[\pi_i(1-\pi_i)]$，$V(p_i)=\pi_i(1-\pi_i)/n_i$，因此工作响应方差为 $1/[n_i\pi_i(1-\pi_i)]$，与前面结果相同。

迭代步骤仍是初始化 β、估计 V 和 μ、计算 η、构造 z、更新 β，直至收敛。大样本下，

$$
E(\hat{\boldsymbol\beta})\approx\boldsymbol\beta,\qquad
V(\hat{\boldsymbol\beta})\approx a(\phi)(\boldsymbol X'\boldsymbol V^{-1}\boldsymbol X)^{-1}.
$$

若不用典则链接，按链式法则求导：

$$
\frac{\partial\ell}{\partial\boldsymbol\beta}
=\sum_i\frac{y_i-\mu_i}{a(\phi)v(\mu_i)}\frac{d\mu_i}{d\eta_i}\boldsymbol x_i,
$$

其中用了 $\partial\ell/\partial\theta_i=(y_i-\mu_i)/a(\phi)$、$d\theta_i/d\mu_i=1/v(\mu_i)$、$\partial\eta_i/\partial\boldsymbol\beta=\boldsymbol x_i$。仍使用相同的工作响应，但 V 的对角元素一般为

$$
V_{ii}=v(\mu_i)\left(\frac{d\eta_i}{d\mu_i}\right)^2.
$$

线性化得分方程和加权最小二乘更新形式不变。非典则链接的一般 IRLS 对应 Fisher 得分迭代，不能一般地与使用观测 Hessian 的 Newton–Raphson 完全等同。[^7]

关于 GLM，有以下几点：

1. 通常的变换分析在变换尺度使用普通最小二乘。
2. GLM 明确承认响应方差不恒定，以加权最小二乘作为参数估计的计算基础。
3. 若变换后仍有方差不齐，合适的 GLM 可能优于标准变换分析。
4. Logistic 回归中的离差、似然比及 Wald 推断也适用于一般 GLM，但应满足相应近似条件。

### S15.3.9 广义线性模型的预测与估计

在关注点 $\boldsymbol x_0$，平均响应估计为

$$
\hat y_0=\hat\mu_0=g^{-1}(\boldsymbol x_0'\hat{\boldsymbol\beta}).
$$

若模型含交互作用等项，$\boldsymbol x_0$ 应展开为相应模型向量。线性预测子的估计方差为 $\boldsymbol x_0'\hat{\boldsymbol\Sigma}\boldsymbol x_0$，其中 Σ̂ 为参数协方差估计。先构造 η 尺度区间

$$
\eta_L=\boldsymbol x_0'\hat{\boldsymbol\beta}-z_{\alpha/2}\sqrt{\boldsymbol x_0'\hat{\boldsymbol\Sigma}\boldsymbol x_0},\qquad
\eta_U=\boldsymbol x_0'\hat{\boldsymbol\beta}+z_{\alpha/2}\sqrt{\boldsymbol x_0'\hat{\boldsymbol\Sigma}\boldsymbol x_0}.
$$

再用逆链接变回响应尺度。逆链接递增时，$L=g^{-1}(\eta_L)$、$U=g^{-1}(\eta_U)$；递减时交换两端，始终以较小者为下限。原材料直接给出的端点顺序只适用于递增逆链接。

SAS PROC GENMOD 用这一方式报告平均响应区间。极大似然估计具有不变性，其函数仍是相应函数的极大似然估计。也可直接用响应尺度上的 Wald 近似，详见 Myers 和 Montgomery（1997）。

### S15.3.10 广义线性模型的残差分析

残差分析用于评价模型、检查假定及链接是否适当。原始残差为 $e_i=y_i-\hat y_i=y_i-\hat\mu_i$。通常推荐**离差残差**：

$$
r_{Di}=\operatorname{sign}(y_i-\hat\mu_i)\sqrt{d_i},
$$

其中 $d_i$ 是第 i 个观测对总离差的贡献。二项 Logistic 模型中，

$$
d_i=2\left[y_i\ln\frac{y_i}{n_i\hat\pi_i}+(n_i-y_i)\ln\frac{1-y_i/n_i}{1-\hat\pi_i}\right],
\qquad\hat\pi_i=\frac1{1+e^{-\boldsymbol x_i'\hat{\boldsymbol\beta}}}.
$$

Poisson 对数链接模型中，

$$
d_i=2\left[y_i\ln\frac{y_i}{\hat\mu_i}-(y_i-\hat\mu_i)\right],
\qquad\hat\mu_i=e^{\boldsymbol x_i'\hat{\boldsymbol\beta}}.
$$

规定 $0\ln0=0$。拟合值接近观测时，离差残差接近零。原补充材料的两条贡献公式缺少与离差定义一致的因子 2，已补全。[^8]

可绘制离差残差的概率图及残差对拟合值图。离差残差不一定正态，尤其是稀疏或二元数据，图形应结合响应分布判断。为使拟合值处于恒定信息尺度，原材料建议：正态用 $\hat y_i$；二项用 $2\arcsin\sqrt{\hat\pi_i}$；Poisson 用 $2\sqrt{\hat y_i}$；gamma 用 $2\ln\hat y_i$。

## S15.4 析因设计中的不平衡数据

正文讨论了几种近似分析法，但也有精确方法，通常利用方差分析与回归的关系。可回顾第 3、5 章及其补充材料。

用例 5.1 的电池寿命试验的修改数据说明。A 为三种材料，B 为三种温度，响应为寿命。表 S15.2 从原数据删除了材料 1 在三个温度下各自最小的观测，又从另外两个单元各随机删除一个观测。

### S15.4.1 回归模型方法

将方差分析模型写成回归模型，采用一般回归显著性检验或额外平方和法。当各单元至少有一个观测时，这一方法较容易应用。

表 S15.2 例 5.1 的修改数据

| 材料 | 15°C | 70°C | 125°C |
| --- | --- | --- | --- |
| 1 | 130，155，180 | 40，80，75 | 70，82，58 |
| 2 | 150，188，159，126 | 136，122，106，115 | 25，70，45 |
| 3 | 138，110，168，160 | 120，150，139 | 96，104，82，60 |

以指示变量编码：材料 1、2、3 分别为 $(x_1,x_2)=(0,0),(1,0),(0,1)$；温度 15、70、125 分别为 $(x_3,x_4)=(0,0),(1,0),(0,1)$。回归模型为

$$
\begin{aligned}
y_{ijk}={}&\beta_0+\beta_1x_1+\beta_2x_2+\beta_3x_3+\beta_4x_4\\
&+\beta_5x_1x_3+\beta_6x_1x_4+\beta_7x_2x_3+\beta_8x_2x_4+\epsilon_{ijk}.
\end{aligned}
$$

$i,j=1,2,3$，$k=1,\ldots,n_{ij}$。五个单元 $(1,1),(1,2),(1,3),(2,3),(3,2)$ 有三次观测，其余有四次，共 31 次。材料和温度主效应各有两个系数、两个自由度，交互作用有四个系数、四个自由度。令 $x_5=x_1x_3,x_6=x_1x_4,x_7=x_2x_3,x_8=x_2x_4$，表 S15.3 给出回归形式的数据。

表 S15.3 回归模型形式的数据

| y | x1 | x2 | x3 | x4 | x5 | x6 | x7 | x8 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 130 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 150 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 136 | 1 | 0 | 1 | 0 | 1 | 0 | 0 | 0 |
| 25 | 1 | 0 | 0 | 1 | 0 | 1 | 0 | 0 |
| 138 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 |
| 96 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 1 |
| 155 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 40 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 |
| 70 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 |
| 188 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 122 | 1 | 0 | 1 | 0 | 1 | 0 | 0 | 0 |
| 70 | 1 | 0 | 0 | 1 | 0 | 1 | 0 | 0 |
| 110 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 |
| 120 | 0 | 1 | 1 | 0 | 0 | 0 | 1 | 0 |
| 104 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 1 |
| 80 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 |
| 82 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 |
| 159 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 106 | 1 | 0 | 1 | 0 | 1 | 0 | 0 | 0 |
| 58 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 |
| 168 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 |
| 150 | 0 | 1 | 1 | 0 | 0 | 0 | 1 | 0 |
| 82 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 1 |
| 180 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 75 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 |
| 126 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 115 | 1 | 0 | 1 | 0 | 1 | 0 | 0 | 0 |
| 45 | 1 | 0 | 0 | 1 | 0 | 1 | 0 | 0 |
| 160 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 |
| 139 | 0 | 1 | 1 | 0 | 0 | 0 | 1 | 0 |
| 60 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 1 |

将含全部八个预测项的模型称为完整模型。其 Minitab 拟合方程为

$$
\hat y=155+0.75x_1-11x_2-90x_3-85x_4+54x_5-24.08x_6+82.33x_7+26.50x_8.
$$

| 参数 | 系数 | 标准误 | t | P 值 |
| --- | ---: | ---: | ---: | ---: |
| 截距 | 155.00 | 12.03 | 12.88 | <0.001 |
| X1 | 0.75 | 15.92 | 0.05 | 0.963 |
| X2 | -11.00 | 15.92 | -0.69 | 0.497 |
| X3 | -90.00 | 17.01 | -5.29 | <0.001 |
| X4 | -85.00 | 17.01 | -5.00 | <0.001 |
| X5 | 54.00 | 22.51 | 2.40 | 0.025 |
| X6 | -24.08 | 23.30 | -1.03 | 0.313 |
| X7 | 82.33 | 23.30 | 3.53 | 0.002 |
| X8 | 26.50 | 22.51 | 1.18 | 0.252 |

残差标准差 S=20.84，$R^2=83.1\%$，调整 $R^2=76.9\%$。

| 来源 | 自由度 | 平方和 | 均方 | F | P 值 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 回归 | 8 | 46814.0 | 5851.8 | 13.48 | <0.001 |
| 残差误差 | 22 | 9553.8 | 434.3 | | |
| 总计 | 30 | 56367.9 | | | |

先检验交互作用：$H_0:\beta_5=\beta_6=\beta_7=\beta_8=0$，备择为至少一个不为零。简化模型只含截距及 $x_1,x_2,x_3,x_4$。拟合为

$$
\hat y=138.02+12.53x_1+23.92x_2-41.91x_3-82.14x_4.
$$

| 参数 | 系数 | 标准误 | t | P 值 |
| --- | ---: | ---: | ---: | ---: |
| 截距 | 138.02 | 11.02 | 12.53 | <0.001 |
| X1 | 12.53 | 11.89 | 1.05 | 0.302 |
| X2 | 23.92 | 11.89 | 2.01 | 0.055 |
| X3 | -41.91 | 11.56 | -3.62 | 0.001 |
| X4 | -82.14 | 11.56 | -7.10 | <0.001 |

S=26.43，$R^2=67.8\%$，调整 $R^2=62.8\%$。

| 来源 | 自由度 | 平方和 | 均方 | F | P 值 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 回归 | 4 | 38212.5 | 9553.1 | 13.68 | <0.001 |
| 残差误差 | 26 | 18155.3 | 698.3 | | |
| 总计 | 30 | 56367.9 | | | |

交互作用的额外平方和为 $46814.0-38212.5=8601.5$，四个自由度，故

$$
F_0=\frac{8601.5/4}{434.3}=4.95.
$$

P 值约为 0.0045，有交互作用证据。

再检验指示变量参数中的 $H_0:\beta_1=\beta_2=0$。简化模型删除 $x_1,x_2$，保留温度和四个交互作用项，拟合为

$$
\hat y=151.273-86.27x_3-81.27x_4+54.75x_5-23.33x_6+71.33x_7+15.50x_8.
$$

| 参数 | 系数 | 标准误 | t | P 值 |
| --- | ---: | ---: | ---: | ---: |
| 截距 | 151.273 | 6.120 | 24.72 | <0.001 |
| X3 | -86.27 | 13.22 | -6.53 | <0.001 |
| X4 | -81.27 | 13.22 | -6.15 | <0.001 |
| X5 | 54.75 | 15.50 | 3.53 | 0.002 |
| X6 | -23.33 | 16.57 | -1.41 | 0.172 |
| X7 | 71.33 | 16.57 | 4.30 | <0.001 |
| X8 | 15.50 | 15.50 | 1.00 | 0.327 |

S=20.30，$R^2=82.5\%$，调整 $R^2=78.1\%$。

| 来源 | 自由度 | 平方和 | 均方 | F | P 值 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 回归 | 6 | 46480.6 | 7746.8 | 18.80 | <0.001 |
| 残差误差 | 24 | 9887.3 | 412.0 | | |
| 总计 | 30 | 56367.9 | | | |

额外平方和为 $46814.0-46480.6=333.4$，$F_0=(333.4/2)/434.3=0.38$，不显著。注意，这一假设比较的是**基准温度下**材料效应，并不等于下一节的第 III 类平均主效应检验。[^9]

对 $H_0:\beta_3=\beta_4=0$，删除温度的两个指示变量主项，保留材料与交互作用，得到

$$
\hat y=96.67+59.08x_1+47.33x_2-36x_5-109.08x_6-7.67x_7-58.50x_8.
$$

| 参数 | 系数 | 标准误 | t | P 值 |
| --- | ---: | ---: | ---: | ---: |
| 截距 | 96.67 | 10.74 | 9.00 | <0.001 |
| X1 | 59.08 | 19.36 | 3.05 | 0.005 |
| X2 | 47.33 | 19.36 | 2.45 | 0.022 |
| X5 | -36.00 | 22.78 | -1.58 | 0.127 |
| X6 | -109.08 | 24.60 | -4.43 | <0.001 |
| X7 | -7.67 | 24.60 | -0.31 | 0.758 |
| X8 | -58.50 | 22.78 | -2.57 | 0.017 |

S=32.21，$R^2=55.8\%$，调整 $R^2=44.8\%$。

| 来源 | 自由度 | 平方和 | 均方 | F | P 值 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 回归 | 6 | 31464 | 5244 | 5.05 | 0.002 |
| 残差误差 | 24 | 24904 | 1038 | | |
| 总计 | 30 | 56368 | | | |

额外平方和为 $46814.0-31464.0=15350.0$，$F_0=(15350.0/2)/434.3=17.67$，P<0.0001。因此，在基准材料下温度效应显著，加上显著交互作用，总体上仍会得到与原平衡数据类似的实际结论。

### S15.4.2 第 III 类分析

另一方法是直接使用第 III 类平方和，即调整平方和。许多软件能完成，例如 Minitab 一般线性模型程序。下面采用每个单元都有观测的情形，材料和温度均为固定因子。

| 来源 | 自由度 | 顺序平方和 | 调整平方和 | 调整均方 | F | P 值 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 材料 | 2 | 2910.4 | 3202.4 | 1601.2 | 3.69 | 0.042 |
| 温度 | 2 | 35302.1 | 36588.7 | 18294.3 | 42.13 | <0.001 |
| 材料 × 温度 | 4 | 8601.5 | 8601.5 | 2150.4 | 4.95 | 0.005 |
| 误差 | 22 | 9553.8 | 9553.8 | 434.3 | | |
| 总计 | 30 | 56367.9 | | | | |

调整平方和是第 III 类平方和，用于 F 检验分子。所检验的等权边际均值假设与平衡情形相对应。误差与交互作用平方和与前述回归分析相同，但主效应平方和不同，因为假设不同。

在各单元均有观测的不平衡设计中，第 III 类分析是常用方法，尤其在研究问题要求各水平等权比较时。参见 Freund、Littell 和 Spector（1988）及 SAS/STAT 手册。方法选择仍应依据所关注假设，而不是只依据平方和名称。

### S15.4.3 第 I、II、III、IV 类平方和

许多软件报告第 I 和第 III 类平方和，SAS 还报告第 II 和第 IV 类。详见 Driscoll 和 Borror（1999）。

**第 I 类**是按模型项进入次序作顺序分解。交互作用应在相应主效应之后加入，嵌套因子按嵌套层次进入。

**第 II 类**衡量某效应在调整其他不包含它的效应后的贡献；例如 AB 包含 A、B，检验 A 时不先调整 AB。不平衡数据下，检验可能依赖单元观测数，未必等于平衡设计的等权假设。对不存在方差分析式过参数化的普通回归模型，第 II 类平方和常是适当的，SAS PROC REG 等报告第 I、II 类。

**第 III、IV 类**常称为偏平方和。平衡设计中四类相同，不平衡时可能不同。对两因子固定效应模型：

- 比例数据中，主效应的第 I=第 II 类，第 III=第 IV 类；交互作用四类相同。若样本数反映目标总体中各水平的比例，按样本数加权的第 I 类可能合适；若希望等权比较，可用第 III 类。
- 不平衡但无空单元，若顺序为 A、B、AB，A 的第 I、II 类通常不同，B 的第 I、II 类相同；主效应的第 III、IV 类相同，AB 的四类相同。等权主效应问题常采用第 III 类。
- 有空单元时，不同类型主效应假设的差别更大。原材料建议考虑第 IV 类，但具体假设由哪些单元缺失决定。某些参数不存在或不能估计，只有可估计函数才能形成可检验假设。必须明确软件实际检验的函数，不能仅凭“第 IV 类”名称判断其含义。SAS PROC GLM 可显示各类平方和的可估计函数，详见 Driscoll 和 Borror（1999）。

### S15.4.4 用均值模型分析不平衡数据

有时放弃通常的效应模型

$$
y_{ijk}=\mu+\tau_i+\beta_j+(\tau\beta)_{ij}+\epsilon_{ijk}
$$

而采用**均值模型** $y_{ijk}=\mu_{ij}+\epsilon_{ijk}$ 更方便，其中 $\mu_{ij}=\mu+\tau_i+\beta_j+(\tau\beta)_{ij}$，$k=1,\ldots,n_{ij}$。有空单元时尤其有用。设空单元数为 m，可将其看作含 ab−m 个处理的单因子模型，每个“处理”对应原析因设计的一个已观测组合。

例如，在表 S15.2 的基础上再删除 $(3,3)$ 单元，形成表 S15.4。材料 3 从未在最高温度下试验，因此没有该组合的信息。

表 S15.4 含一个空单元的电池寿命数据

| 材料 | 15°C | 70°C | 125°C |
| --- | --- | --- | --- |
| 1 | 130，155，180 | 40，80，75 | 70，82，58 |
| 2 | 150，188，159，126 | 136，122，106，115 | 25，70，45 |
| 3 | 138，110，168，160 | 120，150，139 | 无观测 |

以八个单元作为处理的 Minitab 单因素方差分析为

| 来源 | 自由度 | 平方和 | 均方 | F | P 值 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 单元 | 7 | 43843 | 6263 | 14.10 | <0.001 |
| 误差 | 19 | 8439 | 444 | | |
| 总计 | 26 | 52282 | | | |

合并标准差为 21.07。各单元统计量为

| 单元 | n | 均值 | 标准差 |
| --- | ---: | ---: | ---: |
| m11 | 3 | 155.00 | 25.00 |
| m12 | 3 | 65.00 | 21.79 |
| m13 | 3 | 70.00 | 12.00 |
| m21 | 4 | 155.75 | 25.62 |
| m22 | 4 | 119.75 | 12.66 |
| m23 | 3 | 46.67 | 22.55 |
| m31 | 4 | 144.00 | 25.97 |
| m32 | 3 | 136.33 | 15.18 |

原输出还给出基于合并标准差的各均值置信区间，以及 Fisher 成对比较。单次比较错误率 0.05，报告的族错误率 0.453，临界 t 值 2.093。下面完整保留其成对均值差置信区间（前者减后者）。

| 比较 | 下限 | 上限 |
| --- | ---: | ---: |
| m11−m12 | 53.98 | 126.02 |
| m11−m13 | 48.98 | 121.02 |
| m12−m13 | −41.02 | 31.02 |
| m11−m21 | −34.44 | 32.94 |
| m12−m21 | −124.44 | −57.06 |
| m13−m21 | −119.44 | −52.06 |
| m11−m22 | 1.56 | 68.94 |
| m12−m22 | −88.44 | −21.06 |
| m13−m22 | −83.44 | −16.06 |
| m21−m22 | 4.81 | 67.19 |
| m11−m23 | 72.32 | 144.35 |
| m12−m23 | −17.68 | 54.35 |
| m13−m23 | −12.68 | 59.35 |
| m21−m23 | 75.39 | 142.77 |
| m22−m23 | 39.39 | 106.77 |
| m11−m31 | −22.69 | 44.69 |
| m12−m31 | −112.69 | −45.31 |
| m13−m31 | −107.69 | −40.31 |
| m21−m31 | −19.44 | 42.94 |
| m22−m31 | −55.44 | 6.94 |
| m23−m31 | −131.02 | −63.64 |
| m11−m32 | −17.35 | 54.68 |
| m12−m32 | −107.35 | −35.32 |
| m13−m32 | −102.35 | −30.32 |
| m21−m32 | −14.27 | 53.11 |
| m22−m32 | −50.27 | 17.11 |
| m23−m32 | −125.68 | −53.65 |
| m31−m32 | −26.02 | 41.36 |

F=14.10 且 P 很小，说明处理均值有显著差异。原材料列举了七对显著差异：$\mu_{11}\ne\mu_{12},\mu_{11}\ne\mu_{13},\mu_{11}\ne\mu_{22},\mu_{11}\ne\mu_{23},\mu_{21}\ne\mu_{22},\mu_{21}\ne\mu_{23},\mu_{22}\ne\mu_{23}$；完整输出还显示其他显著比较。较长寿命主要来自低温下三种材料，以及中温下材料 2、3。

下一步可构造所关注的单元均值对比。完整九单元设计有四个交互作用自由度，缺一个单元后只剩三个。可选择三个线性独立对比：

$$
\begin{aligned}
C_1&=\mu_{11}-\mu_{13}-\mu_{21}+\mu_{23},\\
C_2&=\mu_{21}-\mu_{22}-\mu_{31}+\mu_{32},\\
C_3&=\mu_{11}-\mu_{12}-\mu_{31}+\mu_{32}.
\end{aligned}
$$

可联合检验三个对比为零，但需更多线性模型知识；这里分别作 t 检验。以样本均值替代 μ，

$$
\hat C_1=155-70-155.75+46.67=-24.08,
$$

$$
V(\hat C_1)=\sigma^2\left(\frac1{n_{11}}+\frac1{n_{13}}+\frac1{n_{21}}+\frac1{n_{23}}\right)
=\frac54\sigma^2.
$$

用 $MS_E=444$，

$$
t_1=\frac{-24.08}{\sqrt{444(5/4)}}=-1.02,\qquad
t_2=\frac{28.33}{\sqrt{444(13/12)}}=1.29,\qquad
t_3=\frac{82.33}{\sqrt{444(5/4)}}=3.49.
$$

只有 $C_3$ 显著，双侧 P 约为 0.0024，表明有材料与温度交互作用的证据。[^10]

结论与第 5 章平衡数据类似：低温下材料差异小，中温下材料 2、3 接近，而材料 1 寿命明显较低；交互作用说明各材料随温度变化的表现不同。原试验能比较高温下三种材料，这里缺少材料 3 的信息。只能说高温下材料 1、2 没有显著差异，不能推断未观测的材料 3。对于温度变化，应按各材料实际已观测均值和相应对比解释。

## 补充参考文献

Driscoll, M. F. and Borror, C. M.（1999），*Sums of Squares and Expected Mean Squares in SAS*，Technical Report，Department of Industrial Engineering，Arizona State University，Tempe，AZ。

Freund, R. J., Littell, R. C., and Spector, P. C.（1988），*The SAS System for Linear Models*，SAS Institute，Cary，NC。

[^1]: 译者注：原补充材料第二个例子的文字写为 $\mu^{1/2}=\sigma^2$，与后面的 $f(t)=t^2$ 和对数变换推导矛盾。对数稳定化要求标准差与均值成正比，即方差与均值平方成正比，现按推导修正。

[^2]: 译者注：原材料的对数似然差界线多印了除以 n；若对数似然定义为 $-n\ln SS_E/2$，界线应为 $\chi_{\alpha,1}^2/2$。这里采用与后续指数界线一致的形式，并将旧版式号对应到正文式（15.2）。

[^3]: 译者注：原材料称 probit 不易包含多个预测变量，这不是该模型的限制；probit 可像 Logistic 一样采用多元线性预测子。

[^4]: 译者注：原材料把优势比解释为成功概率增加量。优势比是优势 $\pi/(1-\pi)$ 的比值，与概率比不同。

[^5]: 译者注：离差的卡方拟合优度近似依赖观测组的信息量；未分组 Bernoulli 数据不满足仅靠总样本量增大即可取得这一近似的条件。参数子集似然比检验有其自身的大样本正则条件。

[^6]: 译者注：补充材料的互补双对数链接遗漏内部负号；指数和 gamma 的倒数链接符号也需与自然参数约定一致。这里均按所列指数族表达式修正。

[^7]: 译者注：非典则链接的一般 IRLS 使用期望信息矩阵，对应 Fisher 得分法；它不一定等于 Newton–Raphson。原材料还将若干工作响应方差近似称为估计线性预测子的方差，这里区分工作响应方差与参数拟合后的预测子方差。

[^8]: 译者注：原材料的二项和 Poisson 观测离差贡献缺少因子 2；为与 $D=2(\ell_{\mathrm{sat}}-\ell_{\mathrm{fit}})$ 一致，已补全。

[^9]: 译者注：含交互作用并采用 0/1 基准编码时，直接删除某主效应系数检验的是基准水平下的简单效应。第 III 类检验的是相应等权边际均值假设，不能把两者当作同一个主效应检验。

[^10]: 译者注：原材料将 $C_1$ 的末项下标误写为 32，实际应为 23；其计算使用的是 46.67，亦即单元 23。$t=3.49$、19 个自由度的双侧 P 约为 0.0024，原文 0.0012 为单侧值。另外，材料 1 的高温均值 70 大于中温均值 65，原文“两个材料高温寿命均低于中、低温”不成立，结论已据数据调整。
