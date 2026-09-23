CHAPTER 3

# Experiments with a Single Factor: The Analysis of Variance

## CHAPTER LEARNING OBJECTIVES

1. Understand how to set up and run a completely randomized experiment.

2. Understand how to perform a single-factor analysis of variance for a completely randomized design.

3. Know the assumptions underlying the ANOVA and how to check for departures from these assumptions.

4. Know how to apply methods for post-ANOVA comparisons for individual differences between means.

5. Know how to interpret computer output from some standard statistics packages.

6. Understand several approaches for determining appropriate sample sizes in designed experiments.

In Chapter 2, we discussed methods for comparing two conditions or treatments. For example, the Portland cement tension bond experiment involved two different mortar formulations. Another way to describe this experiment is as a single-factor experiment with two levels of the factor, where the factor is mortar formulation and the two levels are the two different formulation methods. Many experiments of this type involve more than two levels of the factor. This chapter focuses on methods for the design and analysis of single-factor experiments with an arbitrary number $a$ levels of the factor (or $a$ treatments). We will assume that the experiment has been completely randomized.

## 3.1 An Example

In many integrated circuit manufacturing steps, wafers are completely coated with a layer of material such as silicon dioxide or a metal. The unwanted material is then selectively removed by etching through a mask, thereby creating circuit patterns, electrical interconnects, and areas in which diffusions or metal depositions are to be made. A plasma etching process is widely used for this operation, particularly in small geometry applications. Figure 3.1 shows the important features of a typical single-wafer etching tool. Energy is supplied by a radio-frequency (RF) generator causing plasma to be generated in the gap between the electrodes. The chemical species in the plasma are determined by the particular gases used. Fluorocarbons, such as $CF_{4}$ (tetrafluoromethane) or $C_{2}F_{6}$ (hexafluoroethane), are often used in plasma etching, but other gases and mixtures of gases are relatively common, depending on the application.

■ FIGURE 3.1 A single-wafer plasma etching tool

![](images/figure3.1.jpg)

An engineer is interested in investigating the relationship between the RF power setting and the etch rate for this tool. The objective of an experiment like this is to model the relationship between etch rate and RF power and to specify the power setting that will give a desired target etch rate. She is interested in a particular gas $\left(\mathrm{C}_{2}\mathrm{F}_{6}\right)$ and gap (0.80 cm) and wants to test four levels of RF power: 160, 180, 200, and 220 W. She decided to test five wafers at each level of RF power.

This is an example of a single-factor experiment with $a = 4$ levels of the factor and $n = 5$ replicates. The 20 runs should be made in random order. A very efficient way to generate the run order is to enter the 20 runs in a spreadsheet (Excel), generate a column of random numbers using the RAND() function, and then sort by that column.

Suppose that the test sequence obtained from this process is given as below:

<table><tr><td>Test Sequence</td><td>Excel Random Number (Sorted)</td><td>Power</td></tr><tr><td>1</td><td>12417</td><td>200</td></tr><tr><td>2</td><td>18369</td><td>220</td></tr><tr><td>3</td><td>21238</td><td>220</td></tr><tr><td>4</td><td>24621</td><td>160</td></tr><tr><td>5</td><td>29337</td><td>160</td></tr><tr><td>6</td><td>32318</td><td>180</td></tr><tr><td>7</td><td>36481</td><td>200</td></tr><tr><td>8</td><td>40062</td><td>160</td></tr><tr><td>9</td><td>43289</td><td>180</td></tr><tr><td>10</td><td>49271</td><td>200</td></tr><tr><td>11</td><td>49813</td><td>220</td></tr><tr><td>12</td><td>52286</td><td>220</td></tr><tr><td>13</td><td>57102</td><td>160</td></tr><tr><td>14</td><td>63548</td><td>160</td></tr><tr><td>15</td><td>67710</td><td>220</td></tr><tr><td>16</td><td>71834</td><td>180</td></tr><tr><td>17</td><td>77216</td><td>180</td></tr><tr><td>18</td><td>84675</td><td>180</td></tr><tr><td>19</td><td>89323</td><td>200</td></tr><tr><td>20</td><td>94037</td><td>200</td></tr></table>

This randomized test sequence is necessary to prevent the effects of unknown nuisance variables, perhaps varying out of control during the experiment, from contaminating the results. To illustrate this, suppose that we were to run the 20 test wafers in the original nonrandomized order (that is, all five 160 W power runs are made first, all five 180 W power runs are made next, and so on). If the etching tool exhibits a warm-up effect such that the longer it is on, the lower the observed etch rate readings will be, the warm-up effect will potentially contaminate the data and destroy the validity of the experiment.

Suppose that the engineer runs the experiment that we have designed in the indicated random order. The observations that she obtains on etch rate are shown in Table 3.1.

It is always a good idea to examine experimental data graphically. Figure 3.2a presents box plots for etch rate at each level of RF power and Figure 3.2b presents a scatter diagram of etch rate versus RF power. Both graphs indicate that etch rate increases as the power setting increases. There is no strong evidence to suggest that the variability in etch rate around the average depends on the power setting. On the basis of this simple graphical analysis, we strongly suspect that (1) RF power setting affects the etch rate and (2) higher power settings result in increased etch rate.

Suppose that we wish to be more objective in our analysis of the data. Specifically, suppose that we wish to test for differences between the mean etch rates at all a = 4 levels of RF power. Thus, we are interested in testing the equality of all four means. It might seem that this problem could be solved by performing a t-test for all six possible pairs of means. However, this is not the best solution to this problem. First of all, performing all six pairwise t-tests is inefficient. It takes a lot of effort. Second, conducting all these pairwise comparisons inflates the type I error. Suppose that all four means are equal, so if we select $\alpha = 0.05$ , the probability of reaching the correct decision on any single comparison is 0.95. However, the probability of reaching the correct conclusion on all six comparisons is considerably less than 0.95, so the type I error is inflated.

The appropriate procedure for testing the equality of several means is the analysis of variance. However, the analysis of variance has a much wider application than the problem above. It is probably the most useful technique in the field of statistical inference.

## TABLE 3.1

Etch Rate Data (in Å/min) from the Plasma Etching Experiment

<table><tr><td rowspan="2">Power (W)</td><td colspan="5">Observations</td><td rowspan="2">Totals</td><td rowspan="2">Averages</td></tr><tr><td>1</td><td>2</td><td>3</td><td>4</td><td>5</td></tr><tr><td>160</td><td>575</td><td>542</td><td>530</td><td>539</td><td>570</td><td>2756</td><td>551.2</td></tr><tr><td>180</td><td>565</td><td>593</td><td>590</td><td>579</td><td>610</td><td>2937</td><td>587.4</td></tr><tr><td>200</td><td>600</td><td>651</td><td>610</td><td>637</td><td>629</td><td>3127</td><td>625.4</td></tr><tr><td>220</td><td>725</td><td>700</td><td>715</td><td>685</td><td>710</td><td>3535</td><td>707.0</td></tr></table>

![](images/figure3.2.jpg)  
■ FIGURE 3.2 Box plots and scatter diagram of the etch rate data

![](images/figure3.2.jpg)

## 3.2 The Analysis of Variance

Suppose we have a treatments or different levels of a single factor that we wish to compare. The observed response from each of the a treatments is a random variable. The data would appear as in Table 3.2. An entry in Table 3.2 (e.g., $y_{ij}$ ) represents the jth observation taken under factor level or treatment i. There will be, in general, n observations under the ith treatment. Notice that Table 3.2 is the general case of the data from the plasma etching experiment in Table 3.1.

Models for the Data. We will find it useful to describe the observations from an experiment with a model. One way to write this model is

$$
y _ {i j} = \mu_ {i} + \epsilon_ {i j} \left\{ \begin{array}{l} i = 1, 2, \ldots , a \\ j = 1, 2, \ldots , n \end{array} \right.\tag{3.1}
$$

where $y_{ij}$ is the ijth observation, $\mu_{i}$ is the mean of the ith factor level or treatment, and $\epsilon_{ij}$ is a random error component that incorporates all other sources of variability in the experiment including measurement, variability arising from uncontrolled factors, differences between the experimental units (such as test material) to which the treatments are applied, and the general background noise in the process (such as variability over time, effects of environmental variables). It is convenient to think of the errors as having mean zero, so that $E(y_{ij}) = \mu_{i}$ .

Equation 3.1 is called the means model. An alternative way to write a model for the data is to define

$$
\mu_ {i} = \mu + \tau_ {i}, \quad i = 1, 2, \dots , a
$$

so that Equation 3.1 becomes

$$
y _ {i j} = \mu + \tau_ {i} + \epsilon_ {i j} \left\{ \begin{array}{l} i = 1, 2, \ldots , a \\ j = 1, 2, \ldots , n \end{array} \right.\tag{3.2}
$$

In this form of the model, $\mu$ is a parameter common to all treatments called the overall mean, and $\tau_{i}$ is a parameter unique to the ith treatment called the ith treatment effect. Equation 3.2 is usually called the effects model.

Both the means model and the effects model are linear statistical models; that is, the response variable $y_{ij}$ is a linear function of the model parameters. Although both forms of the model are useful, the effects model is more widely encountered in the experimental design literature. It has some intuitive appeal in that $\mu$ is a constant and the treatment effects $\tau_{i}$ represent deviations from this constant when the specific treatments are applied.

Equation 3.2 (or 3.1) is also called the one-way or single-factor analysis of variance (ANOVA) model because only one factor is investigated. Furthermore, we will require that the experiment be performed in random order so that the environment in which the treatments are applied (often called the experimental units) is as uniform as possible. Thus, the experimental design is a completely randomized design. Our objectives will be to test appropriate hypotheses about the treatment means and to estimate them. For hypothesis testing, the model errors are assumed to be normally and independently distributed random variables with mean zero and variance $\sigma^{2}$ . The variance $\sigma^{2}$ is assumed to be constant for all levels of the factor. This implies that the observations

$$
y _ {i j} \sim N (\mu + \tau_ {i}, \sigma^ {2})
$$

and that the observations are mutually independent.

## TABLE 3.2

Typical Data for a Single-Factor Experiment

<table><tr><td>Treatment (Level)</td><td></td><td>Observations</td><td></td><td></td><td>Totals</td><td>Averages</td></tr><tr><td>1</td><td> $y_{11}$ </td><td> $y_{12}$ </td><td> $\cdots$ </td><td> $y_{1n}$ </td><td> $y_{1.}$ </td><td> $\overline{y}_{1.}$ </td></tr><tr><td>2</td><td> $y_{21}$ </td><td> $y_{22}$ </td><td> $\cdots$ </td><td> $y_{2n}$ </td><td> $y_{2.}$ </td><td> $\overline{y}_{2.}$ </td></tr><tr><td> $\vdots$ </td><td> $\vdots$ </td><td> $\vdots$ </td><td> $\cdots$ </td><td> $\vdots$ </td><td> $\vdots$ </td><td> $\vdots$ </td></tr><tr><td>a</td><td> $y_{a1}$ </td><td> $y_{a2}$ </td><td> $\cdots$ </td><td> $y_{an}$ </td><td> $\frac{y_{a.}}{y..}$ </td><td> $\frac{\overline{y}_{a.}}{y..}$ </td></tr></table>

Fixed or Random Factor? The statistical model, Equation 3.2, describes two different situations with respect to the treatment effects. First, the a treatments could have been specifically chosen by the experimenter. In this situation, we wish to test hypotheses about the treatment means, and our conclusions will apply only to the factor levels considered in the analysis. The conclusions cannot be extended to similar treatments that were not explicitly considered. We may also wish to estimate the model parameters $(\mu, \tau_{i}, \sigma^{2})$ . This is called the fixed effects model. Alternatively, the a treatments could be a random sample from a larger population of treatments. In this situation, we should like to be able to extend the conclusions (which are based on the sample of treatments) to all treatments in the population, whether or not they were explicitly considered in the analysis. Here, the $\tau_{i}$ are random variables, and knowledge about the particular ones investigated is relatively useless. Instead, we test hypotheses about the variability of the $\tau_{i}$ and try to estimate this variability. This is called the random effects model or components of variance model. We discuss the single-factor random effects model in Section 3.9. However, we will defer a more complete discussion of experiments with random factors to Chapter 13.

## 3.3 Analysis of the Fixed Effects Model

In this section, we develop the single-factor analysis of variance for the fixed effects model. Recall that $y_{i}$ represents the total of the observations under the ith treatment. Let $\overline{y}_{i}$ represent the average of the observations under the ith treatment. Similarly, let $y_{..}$ represent the grand total of all the observations and $\overline{y}_{..}$ represent the grand average of all the observations. Expressed symbolically,

$$
y _ {i.} = \sum_ {j = 1} ^ {n} y _ {i j} \quad \overline {{{{y}}}} _ {i.} = y _ {i.} / n \quad i = 1, 2, \dots , a\tag{3.3}
$$

$$
y _ {\cdot \cdot} = \sum_ {i = 1} ^ {a} \sum_ {j = 1} ^ {n} y _ {i j} \quad \overline {{y}} _ {\cdot \cdot} = y _ {\cdot \cdot} / N
$$

where N = an is the total number of observations. We see that the “dot” subscript notation implies summation over the subscript that it replaces.

We are interested in testing the equality of the $a$ treatment means; that is, $E(y_{ij}) = \mu + \tau_i = \mu_i$ , $i = 1, 2, \ldots, a$ . The appropriate hypotheses are

$$
\begin{array}{l} H _ {0}: \mu_ {1} = \mu_ {2} = \dots = \mu_ {a} \\ H _ {1}: \mu_ {i} \neq \mu_ {j} \quad \text { for   at   least   one   pair } (i, j) \end{array}\tag{3.4}
$$

In the effects model, we break the ith treatment mean $\mu_{i}$ into two components such that $\mu_{i} = \mu + \tau_{i}$ . We usually think of $\mu$ as an overall mean so that

$$
\frac {\sum_ {i = 1} ^ {a} \mu_ {i}}{a} = \mu
$$

This definition implies that

$$
\sum_ {i = 1} ^ {a} \tau_ {i} = 0
$$

That is, the treatment or factor effects can be thought of as deviations from the overall mean. $^{1}$ Consequently, an equivalent way to write the above hypotheses is in terms of the treatment effects $\tau_{i}$ , say

$$
\begin{array}{l} H _ {0}: \tau_ {1} = \tau_ {2} = \dots \tau_ {a} = 0 \\ H _ {1}: \tau_ {i} \neq 0 \quad \text { for   at   least   one } i \end{array}
$$

Thus, we speak of testing the equality of treatment means or testing that the treatment effects (the $\tau_{i}$ ) are zero. The appropriate procedure for testing the equality of a treatment means is the analysis of variance.

## 3.3.1 Decomposition of the Total Sum of Squares

The name analysis of variance is derived from a partitioning of total variability into its component parts. The total corrected sum of squares

$$
S S _ {T} = \sum_ {i = 1} ^ {a} \sum_ {j = 1} ^ {n} (y _ {i j} - \overline {{y}} _ {\cdot \cdot}) ^ {2}
$$

is used as a measure of overall variability in the data. Intuitively, this is reasonable because if we were to divide $SS_{T}$ by the appropriate number of degrees of freedom (in this case, an - 1 = N - 1), we would have the sample variance of the y's. The sample variance is, of course, a standard measure of variability.

Note that the total corrected sum of squares $SS_{T}$ may be written as

$$
\sum_ {i = 1} ^ {a} \sum_ {j = 1} ^ {n} (y _ {i j} - \overline {{y}} _ {\cdot \cdot}) ^ {2} = \sum_ {i = 1} ^ {a} \sum_ {j = 1} ^ {n} [ (\overline {{y}} _ {i.} - \overline {{y}} _ {\cdot \cdot}) + (y _ {i j} - \overline {{y}} _ {i.}) ] ^ {2}\tag{3.5}
$$

or

$$
\begin{array}{r} \sum_ {i = 1} ^ {a} \sum_ {j = 1} ^ {n} (y _ {i j} - \overline {{y}} _ {\cdot \cdot}) ^ {2} = n \sum_ {i = 1} ^ {a} (\overline {{y}} _ {i.} - \overline {{y}} _ {\cdot \cdot}) ^ {2} + \sum_ {i = 1} ^ {a} \sum_ {j = 1} ^ {n} (y _ {i j} - \overline {{y}} _ {i.}) ^ {2} \\ + 2 \sum_ {i = 1} ^ {a} \sum_ {j = 1} ^ {n} (\overline {{y}} _ {i.} - \overline {{y}} _ {\cdot \cdot}) (y _ {i j} - \overline {{y}} _ {i.}) \end{array}
$$

However, the cross-product term in this last equation is zero, because

$$
\sum_ {j = 1} ^ {n} (y _ {i j} - \overline {{{{y}}}} _ {i.}) = y _ {i.} - n \overline {{{{y}}}} _ {i.} = y _ {i.} - n (y _ {i.} / n) = 0
$$

Therefore, we have

$$
\sum_ {i = 1} ^ {a} \sum_ {j = 1} ^ {n} (y _ {i j} - \bar {y} _ {\cdot}) ^ {2} = n \sum_ {i = 1} ^ {a} (\bar {y} _ {i.} - \bar {y} _ {\cdot}) ^ {2} + \sum_ {i = 1} ^ {a} \sum_ {j = 1} ^ {n} (y _ {i j} - \bar {y} _ {i.}) ^ {2}\tag{3.6}
$$

Equation 3.6 is the fundamental ANOVA identity. It states that the total variability in the data, as measured by the total corrected sum of squares, can be partitioned into a sum of squares of the differences between the treatment averages and the grand average plus a sum of squares of the differences of observations within treatments from the treatment average. Now, the difference between the observed treatment averages and the grand average is a measure of the differences between treatment means, whereas the differences of observations within a treatment from the treatment average can be due to only random error. Thus, we may write Equation 3.6 symbolically as

$$
S S _ {T} = S S _ {\mathrm{Treatments}} + S S _ {E}
$$

where $SS_{Treatments}$ is called the sum of squares due to treatments (i.e., between treatments) and $SS_{E}$ is called the sum of squares due to error (i.e., within treatments). There are an = N total observations; thus, $SS_{T}$ has N - 1 degrees of freedom. There are a levels of the factor (and a treatment means), so $SS_{Treatments}$ has a - 1 degrees of freedom. Finally, there are n replicates within any treatment providing n - 1 degrees of freedom with which to estimate the experimental error. Because there are a treatments, we have $a(n - 1) = an - a = N - a$ degrees of freedom for error.

It is instructive to examine explicitly the two terms on the right-hand side of the fundamental ANOVA identity. Consider the error sum of squares

$$
S S _ {E} = \sum_ {i = 1} ^ {a} \sum_ {j = 1} ^ {n} (y _ {i j} - \overline {{y}} _ {i.}) ^ {2} = \sum_ {i = 1} ^ {a} \left[ \sum_ {j = 1} ^ {n} (y _ {i j} - \overline {{y}} _ {i.}) ^ {2} \right]
$$

$$
S _ {i} ^ {2} = \frac {\sum_ {j = 1} ^ {n} (y _ {i j} - \overline {{y}} _ {i .}) ^ {2}}{n - 1} \qquad i = 1, 2, \ldots , a
$$

In this form, it is easy to see that the term within square brackets, if divided by $n - 1$ , is the sample variance in the $i$ th treatment, or

Now a sample variances may be combined to give a single estimate of the common population variance as follows:

$$
\begin{array}{r l} \frac {(n - 1) S _ {1} ^ {2} + (n - 1) S _ {2} ^ {2} + \cdots + (n - 1) S _ {a} ^ {2}}{(n - 1) + (n - 1) + \cdots + (n - 1)} & = \frac {\sum_ {i = 1} ^ {a} \left[ \sum_ {j = 1} ^ {n} (y _ {i j} - \overline {{y}} _ {i .}) ^ {2} \right]}{\sum_ {i = 1} ^ {a} (n - 1)} \\ & = \frac {S S _ {E}}{(N - a)} \end{array}
$$

Thus, $SS_{E} / (N - a)$ is a pooled estimate of the common variance within each of the $a$ treatments.

Similarly, if there were no differences between the a treatment means, we could use the variation of the treatment averages from the grand average to estimate $\sigma^{2}$ . Specifically,

$$
\frac {S S _ {\text {Treatments}}}{a - 1} = \frac {n \sum_ {i = 1} ^ {a} (\overline {{y}} _ {i .} - \overline {{y}} _ {. .}) ^ {2}}{a - 1}
$$

is an estimate of $\sigma^2$ if the treatment means are equal. The reason for this may be intuitively seen as follows: The quantity $\sum_{i=1}^{a} (\overline{y}_i - \overline{y}_{\cdot})^2 / (a-1)$ estimates $\sigma^2 / n$ , the variance of the treatment averages, so $n \sum_{i=1}^{a} (\overline{y}_{i.} - \overline{y}_{\cdot})^2 / (a-1)$ must estimate $\sigma^2$ if there are no differences in treatment means.

We see that the ANOVA identity (Equation 3.6) provides us with two estimates of $\sigma^{2}$ —one based on the inherent variability within treatments and the other based on the variability between treatments. If there are no differences in the treatment means, these two estimates should be very similar, and if they are not, we suspect that the observed difference must be caused by differences in the treatment means. Although we have used an intuitive argument to develop this result, a somewhat more formal approach can be taken.

The quantities

$$
M S _ {\text { Treatments }} = \frac {S S _ {\text { Treatments }}}{a - 1}
$$

and

$$
M S _ {E} = \frac {S S _ {E}}{N - a}
$$

are called mean squares. We now examine the expected values of these mean squares. Consider

$$
\begin{array}{r l} & E (M S _ {E}) = E \left(\frac {S S _ {E}}{N - a}\right) = \frac {1}{N - a} E \left[ \sum_ {i = 1} ^ {a} \sum_ {j = 1} ^ {n} (y _ {i j} - \overline {{y}} _ {i.}) ^ {2} \right] \\ & \qquad = \frac {1}{N - a} E \left[ \sum_ {i = 1} ^ {a} \sum_ {j = 1} ^ {n} (y _ {i j} ^ {2} - 2 y _ {i j} \overline {{y}} _ {i.} + \overline {{y}} _ {i.} ^ {2}) \right] \\ & \qquad = \frac {1}{N - a} E \left[ \sum_ {i = 1} ^ {a} \sum_ {j = 1} ^ {n} y _ {i j} ^ {2} - 2 n \sum_ {i = 1} ^ {a} \overline {{y}} _ {i.} ^ {2} + n \sum_ {i = 1} ^ {a} \overline {{y}} _ {i.} ^ {2} \right] \\ & \qquad = \frac {1}{N - a} E \left[ \sum_ {i = 1} ^ {a} \sum_ {j = 1} ^ {n} y _ {i j} ^ {2} - \frac {1}{n} \sum_ {i = 1} ^ {a} \overline {{y}} _ {i.} ^ {2} \right] \end{array}
$$

Substituting the model (Equation 3.1) into this equation, we obtain

$$
E (M S _ {E}) = \frac {1}{N - a} E \left[ \sum_ {i = 1} ^ {a} \sum_ {j = 1} ^ {n} (\mu + \tau_ {i} + \epsilon_ {i j}) ^ {2} - \frac {1}{n} \sum_ {i = 1} ^ {a} \left(\sum_ {i = 1} ^ {n} \mu + \tau_ {i} + \epsilon_ {i j}\right) ^ {2} \right]
$$

Now when squaring and taking expectation of the quantity within the brackets, we see that terms involving $\epsilon_{ij}^{2}$ and $\epsilon_{i}^{2}$ are replaced by $\sigma^{2}$ and $n\sigma^{2}$ , respectively, because $E(\epsilon_{ij}) = 0$ . Furthermore, all cross products involving $\epsilon_{ij}$ have zero expectation. Therefore, after squaring and taking expectation, the last equation becomes

$$
E (M S _ {E}) = \frac {1}{N - a} \left[ N \mu^ {2} + n \sum_ {i = 1} ^ {a} \tau_ {i} ^ {2} + N \sigma^ {2} - N \mu^ {2} - n \sum_ {i = 1} ^ {a} \tau_ {i} ^ {2} - a \sigma^ {2} \right]
$$

or

$$
E (M S _ {E}) = \sigma^ {2}
$$

By a similar approach, we may also show that $^{2}$

$$
E (M S _ {\text { Treatments }}) = \sigma^ {2} + \frac {n \sum_ {i = 1} ^ {a} \tau_ {i} ^ {2}}{a - 1}
$$

Thus, as we argued heuristically, $MS_{E} = SS_{E}/(N - a)$ estimates $\sigma^{2}$ , and, if there are no differences in treatment means (which implies that $\tau_{i} = 0$ ), $MS_{Treatments} = SS_{Treatments}/(a - 1)$ also estimates $\sigma^{2}$ . However, note that if treatment means do differ, the expected value of the treatment mean square is greater than $\sigma^{2}$ .

It seems clear that a test of the hypothesis of no difference in treatment means can be performed by comparing $MS_{Treatments}$ and $MS_{E}$ . We now consider how this comparison may be made.

## 3.3.2 Statistical Analysis

We now investigate how a formal test of the hypothesis of no differences in treatment means $(H_{0}:\mu_{1}=\mu_{2}=\cdots=\mu_{a},$ or equivalently, $H_{0}:\tau_{1}=\tau_{2}=\cdots=\tau_{a}=0)$ can be performed. Because we have assumed that the errors $\epsilon_{ij}$ are normally and independently distributed with mean zero and variance $\sigma^{2}$ , the observations $y_{ij}$ are normally and independently distributed with mean $\mu+\tau_{i}$ and variance $\sigma^{2}$ . Thus, $SS_{T}$ is a sum of squares in normally distributed random variables; consequently, it can be shown that $SS_{T}/\sigma^{2}$ is distributed as chi-square with N-1 degrees of freedom. Furthermore, we can show that $SS_{E}/\sigma^{2}$ is chi-square with N-a degrees of freedom and that $SS_{Treatments}/\sigma^{2}$ is chi-square with a-1 degrees of freedom if the null hypothesis $H_{0}:\tau_{i}=0$ is true. However, all three sums of squares are not necessarily independent because $SS_{Treatments}$ and $SS_{E}$ add to $SS_{T}$ . The following theorem, which is a special form of one attributed to William G. Cochran, is useful in establishing the independence of $SS_{E}$ and $SS_{Treatments}$ .

## THEOREM 3-1 Cochran's Theorem

Let $Z_{i}$ be NID(0, 1) for $i = 1, 2, \ldots, v$ and

$$
\sum_ {i = 1} ^ {v} Z _ {i} ^ {2} = Q _ {1} + Q _ {2} + \dots + Q _ {s}
$$

where $s \leq v$ , and $Q_{i}$ has $v_{i}$ degrees of freedom ( $i = 1, 2, \ldots, s$ ). Then $Q_{1}, Q_{2}, \ldots, Q_{s}$ are independent chi-square random variables with $v_{1}, v_{2}, \ldots, v_{s}$ degrees of freedom, respectively, if and only if

$$
v = v _ {1} + v _ {2} + \dots + v _ {s}
$$

Because the degrees of freedom for $SS_{\text{Treatments}}$ and $SS_E$ add to $N - 1$ , the total number of degrees of freedom, Cochran's theorem implies that $SS_{\text{Treatments}} / \sigma^2$ and $SS_E / \sigma^2$ are independently distributed chi-square random variables. Therefore, if the null hypothesis of no difference in treatment means is true, the ratio

$$
F _ {0} = \frac {S S _ {\text { Treatments }} / (a - 1)}{S S _ {E} / (N - a)} = \frac {M S _ {\text { Treatments }}}{M S _ {E}}\tag{3.7}
$$

is distributed as F with a - 1 and N - a degrees of freedom. Equation 3.7 is the test statistic for the hypothesis of no differences in treatment means.

From the expected mean squares we see that, in general, $MS_{E}$ is an unbiased estimator of $\sigma^{2}$ . Also, under the null hypothesis, $MS_{Treatments}$ is an unbiased estimator of $\sigma^{2}$ . However, if the null hypothesis is false, the expected value of $MS_{Treatments}$ is greater than $\sigma^{2}$ . Therefore, under the alternative hypothesis, the expected value of the numerator of the test statistic (Equation 3.7) is greater than the expected value of the denominator, and we should reject $H_{0}$ on values of the test statistic that are too large. This implies an upper-tail, one-tail critical region. Therefore, we should reject $H_{0}$ and conclude that there are differences in the treatment means if

$$
F _ {0} > F _ {\alpha , a - 1, N - a}
$$

where $F_{0}$ is computed from Equation 3.7. Alternatively, we could use the P-value approach for decision making. The table of F percentages in the Appendix (Table IV) can be used to find bounds on the P-value.

The sums of squares may be computed in several ways. One direct approach is to make use of the definition

$$
y _ {i j} - \bar {y} _ {\cdot \cdot} = (\bar {y} _ {\cdot} - \bar {y} _ {\cdot \cdot}) + (y _ {i j} - \bar {y} _ {i.})
$$

Use a spreadsheet to compute these three terms for each observation. Then, sum up the squares to obtain $SS_{T}$ , $SS_{Treatments}$ , and $SS_{E}$ . Another approach is to rewrite and simplify the definitions of $SS_{Treatments}$ and $SS_{T}$ in Equation 3.6, which results in

$$
S S _ {T} = \sum_ {i = 1} ^ {a} \sum_ {j = 1} ^ {n} y _ {i j} ^ {2} - \frac {y _ {. .} ^ {2}}{N}\tag{3.8}
$$

$$
S S _ {\text { Treatments }} = \frac {1}{n} \sum_ {i = 1} ^ {a} y _ {i.} ^ {2} - \frac {y _ {. .} ^ {2}}{N}\tag{3.9}
$$

and

$$
S S _ {E} = S S _ {T} - S S _ {\mathrm{Treatments}}\tag{3.10}
$$

This approach is nice because some calculators are designed to accumulate the sum of entered numbers in one register and the sum of the squares of those numbers in another, so each number only has to be entered once. In practice, we use computer software to do this.

The test procedure is summarized in Table 3.3. This is called an analysis of variance (or ANOVA) table.

## TABLE 3.3

The Analysis of Variance Table for the Single-Factor, Fixed Effects Model

<table><tr><td>Source of Variation</td><td>Sum of Squares</td><td>Degrees of Freedom</td><td>Mean Square</td><td> $F_0$ </td></tr><tr><td>Between treatments</td><td> $SS_{Treatments}=n\sum_{i=1}^{a}(\overline{y}_{i.}-\overline{y}_{..})^2$ </td><td>a-1</td><td> $MS_{Treatments}$ </td><td> $F_0=\frac{MS_{Treatments}}{MS_E}$ </td></tr><tr><td>Error (within treatments)</td><td> $SS_E=SS_T-SS_{Treatments}$ </td><td>N-a</td><td> $MS_E$ </td><td></td></tr><tr><td>Total</td><td> $SS_T=\sum_{i=1}^{a}\sum_{j=1}^{n}(y_{ij}-\overline{y}_{..})^2$ </td><td>N-1</td><td></td><td></td></tr></table>

## EXAMPLE 3.1 The Plasma Etching Experiment

To illustrate the analysis of variance, return to the first example discussed in Section 3.1. Recall that the engineer is interested in determining if the RF power setting affects the etch rate, and she has run a completely randomized experiment with four levels of RF power and five replicates. For convenience, we repeat here the data from Table 3.1:

Usually, these calculations would be performed on a computer, using a software package with the capability to analyze data from designed experiments.

<table><tr><td rowspan="2">RF Power (W)</td><td colspan="5">Observed Etch Rate (Å/min)</td><td rowspan="2">Totals $y_{i.}$ </td><td rowspan="2">Averages $\overline{y}_{i.}$ </td></tr><tr><td>1</td><td>2</td><td>3</td><td>4</td><td>5</td></tr><tr><td>160</td><td>575</td><td>542</td><td>530</td><td>539</td><td>570</td><td>2756</td><td>551.2</td></tr><tr><td>180</td><td>565</td><td>593</td><td>590</td><td>579</td><td>610</td><td>2937</td><td>587.4</td></tr><tr><td>200</td><td>600</td><td>651</td><td>610</td><td>637</td><td>629</td><td>3127</td><td>625.4</td></tr><tr><td>220</td><td>725</td><td>700</td><td>715</td><td>685</td><td>710</td><td>3535</td><td>707.0</td></tr><tr><td></td><td></td><td></td><td></td><td></td><td></td><td> $y_{i.}=12,355$ </td><td> $\overline{y}_{..}=617.75$ </td></tr></table>

We will use the analysis of variance to test $H_{0}:\mu_{1}=\mu_{2}=\mu_{3}=\mu_{4}$ against the alternative $H_{1}$ : some means are different. The sums of squares required are computed using Equations 3.8, 3.9, and 3.10 as follows:

$$
\begin{array}{r l} S S _ {T} & = \sum_ {i = 1} ^ {4} \sum_ {j = 1} ^ {5} y _ {i j} ^ {2} - \frac {y _ {. .} ^ {2}}{N} \\ & = (5 7 5) ^ {2} + (5 4 2) ^ {2} + \dots + (7 1 0) ^ {2} - \frac {(1 2 , 3 5 5) ^ {2}}{2 0} \\ & = 7 2, 2 0 9. 7 5 \end{array}
$$

$$
\begin{array}{r l} S S _ {\text {Treatments}} & = \frac {1}{n} \sum_ {i = 1} ^ {4} y _ {i.} ^ {2} - \frac {y _ {. .} ^ {2}}{N} \\ & = \frac {1}{5} [ (2 7 5 6) ^ {2} + \dots + (3 5 3 5) ^ {2} ] - \frac {(1 2 , 3 5 5) ^ {2}}{2 0} \\ & = 6 6, 8 7 0. 5 5 \\ S S _ {E} & = S S _ {T} - S S _ {\text {Treatments}} \\ & = 7 2, 2 0 9. 7 5 - 6 6, 8 7 0. 5 5 = 5 3 3 9. 2 0 \end{array}
$$

The ANOVA is summarized in Table 3.4. Note that the RF power or between-treatment mean square (22,290.18) is many times larger than the within-treatment or error mean square (333.70). This indicates that it is unlikely that the treatment means are equal. More formally, we can compute the F ratio $F_{0}=22,290.18/333.70=66.80$ and compare this to an appropriate upper-tail percentage point of the $F_{3,16}$ distribution. To use a fixed significance level approach, suppose that the experimenter has selected $\alpha=0.05$ . From Appendix Table IV, we find that $F_{0.05,3,16}=3.24$ . Because $F_{0}=66.80>3.24$ , we reject $H_{0}$ and conclude that the treatment means differ; that is, the RF power setting significantly affects the mean etch rate. We could also compute a P-value for this test statistic. Figure 3.3 shows the reference distribution ( $F_{3,16}$ ) for the test statistic $F_{0}$ . Clearly, the P-value is very small in this case. From Appendix Table A-4, we find that $F_{0.01,3,16}=5.29$ and because $F_{0}>5.29$ , we can conclude that an upper bound for the P-value is 0.01; that is, P<0.01 (the exact P-value is $P=2.88\times10^{-9}$ ).

TABLE 3.4  
ANOVA for the Plasma Etching Experiment

<table><tr><td>Source of Variation</td><td>Sum of Square</td><td>Degrees of Freedom</td><td>Mean Squares</td><td> $F_0$ </td><td>P-Value</td></tr><tr><td>RF Power</td><td>66,870.55</td><td>3</td><td>22,290.18</td><td> $F_0 = 66.80$ </td><td>&lt; 0.01</td></tr><tr><td>Error</td><td>5339.20</td><td>16</td><td>333.70</td><td></td><td></td></tr><tr><td>Total</td><td>72,209.75</td><td>19</td><td></td><td></td><td></td></tr></table>

![](images/figure3.3.jpg)  
■ FIGURE 3.3 The reference distribution $(F_{3,16})$ for the test statistic $F_{0}$ in Example 3.1

Coding the Data. Generally, we need not be too concerned with computing because there are many widely available computer programs for performing the calculations. These computer programs are also helpful in performing many other analyses associated with experimental design (such as residual analysis and model adequacy checking). In many cases, these programs will also assist the experimenter in setting up the design.

However, when hand calculations are necessary, it is sometimes helpful to code the observations. This is illustrated in Example 3.2.

## EXAMPLE 3.2 Coding the Observations

The ANOVA calculations may often be made more easily or accurately by coding the observations. For example, consider the plasma etching data in Example 3.1. Suppose that we subtract 600 from each observation. The coded data are shown in Table 3.5. It is easy to verify that

$$
\begin{array}{r l} & S S _ {T} = (- 2 5) ^ {2} + (- 5 8) ^ {2} + \dots \\ & \qquad + (1 1 0) ^ {2} - \frac {(3 5 5) ^ {2}}{2 0} = 7 2, 2 0 9. 7 5 \end{array}
$$

$$
\begin{array}{r l} S S _ {\text { Treatment }} & = \frac {(- 2 4 4) ^ {2} + (- 6 3) ^ {2} + (1 2 7) ^ {2} + (5 3 5) ^ {2}}{5} \\ & - \frac {(3 5 5) ^ {2}}{2 0} = 6 6, 8 7 0. 5 5 \end{array}
$$

and

$$
S S _ {E} = 5 3 3 9. 2 0
$$

Comparing these sums of squares to those obtained in Example 3.1, we see that subtracting a constant from the original data does not change the sums of squares.

Now suppose that we multiply each observation in Example 3.1 by 2. It is easy to verify that the sums of squares for the transformed data are $SS_{T}=288,839.00$ , $SS_{Treatments}=267,482.20$ , and $SS_{E}=21,356.80$ . These sums of squares appear to differ considerably from those obtained in Example 3.1. However, if they are divided by 4 (i.e., $2^{2}$ ), the results are identical. For example, for the treatment sum of squares 267,482.20/4=66,870.55. Also, for the coded data, the F-ratio is $F=(267,482.20/3)/(21,356.80/16)=66.80$ , which is identical to the F-ratio for the original data. Thus, the ANOVAs are equivalent.

TABLE 3.5  
Coded Etch Rate Data for Example 3.2

<table><tr><td rowspan="2">RF Power (W)</td><td colspan="5">Observations</td><td rowspan="2">Totals  $y_{i.}$ </td></tr><tr><td>1</td><td>2</td><td>3</td><td>4</td><td>5</td></tr><tr><td>160</td><td>-25</td><td>-58</td><td>-70</td><td>-61</td><td>-30</td><td>-244</td></tr><tr><td>180</td><td>-35</td><td>-7</td><td>-10</td><td>-21</td><td>10</td><td>-63</td></tr><tr><td>200</td><td>0</td><td>51</td><td>10</td><td>37</td><td>29</td><td>127</td></tr><tr><td>220</td><td>125</td><td>100</td><td>115</td><td>85</td><td>110</td><td>535</td></tr></table>

Randomization Tests and Analysis of Variance. In our development of the ANOVA F-test, we have used the assumption that the random errors $\epsilon_{ij}$ are normally and independently distributed random variables. The F-test can also be justified as an approximation to a randomization test. To illustrate this, suppose that we have five observations on each of two treatments and that we wish to test the equality of treatment means. The data would look like this:

<table><tr><td>Treatment 1</td><td>Treatment 2</td></tr><tr><td> $y_{11}$ </td><td> $y_{21}$ </td></tr><tr><td> $y_{12}$ </td><td> $y_{22}$ </td></tr><tr><td> $y_{13}$ </td><td> $y_{23}$ </td></tr><tr><td> $y_{14}$ </td><td> $y_{24}$ </td></tr><tr><td> $y_{15}$ </td><td> $y_{25}$ </td></tr></table>

We could use the ANOVA F-test to test $H_{0}:\mu_{1}=\mu_{2}$ . Alternatively, we could use a somewhat different approach. Suppose that we consider all the possible ways of allocating the 10 numbers in the above sample to the two treatments. There are $10!/5!5!=252$ possible arrangements of the 10 observations. If there is no difference in treatment means, all 252 arrangements are equally likely. For each of the 252 arrangements, we calculate the value of the F-statistic using Equation 3.7. The distribution of these F values is called a randomization distribution, and a large value of F indicates that the data are not consistent with the hypothesis $H_{0}:\mu_{1}=\mu_{2}$ . For example, if the value of F actually observed was exceeded by only five of the values of the randomization distribution, this would correspond to rejection of $H_{0}:\mu_{1}=\mu_{2}$ at a significance level of $\alpha=5/252=0.0198$ (or 1.98 percent). Notice that no normality assumption is required in this approach.

The difficulty with this approach is that, even for relatively small problems, it is computationally prohibitive to enumerate the exact randomization distribution. However, numerous studies have shown that the exact randomization distribution is well approximated by the usual normal-theory F distribution. Thus, even without the normality assumption, the ANOVA F-test can be viewed as an approximation to the randomization test. For further reading on randomization tests in the analysis of variance, see Box, Hunter, and Hunter (2005).

## 3.3.3 Estimation of the Model Parameters

We now present estimators for the parameters in the single-factor model

$$
y _ {i j} = \mu + \tau_ {i} + \epsilon_ {i j}.
$$

and confidence intervals on the treatment means. We will prove later that reasonable estimates of the overall mean and the treatment effects are given by

$$
\begin{array}{l} \hat {\mu} = \overline {{y}} _ {..} \\ \hat {\tau} _ {i} = \overline {{y}} _ {i.} - \overline {{y}} _ {..}, i = 1, 2, \ldots , a \end{array}\tag{3.11}
$$

These estimators have considerable intuitive appeal; note that the overall mean is estimated by the grand average of the observations and that any treatment effect is just the difference between the treatment average and the grand average.

A confidence interval estimate of the ith treatment mean may be easily determined. The mean of the ith treatment is

$$
\mu_ {i} = \mu + \tau_ {i}
$$

A point estimator of $\mu_{i}$ would be $\hat{\mu}_{i} = \hat{\mu} + \hat{\tau}_{i} = \overline{y}_{i}$ . Now, if we assume that the errors are normally distributed, each treatment average $\overline{y}_{i}$ is distributed $\mathrm{NID}(\mu_{i}, \sigma^{2}/n)$ . Thus, if $\sigma^{2}$ were known, we could use the normal distribution to define the confidence interval. Using the $MS_{E}$ as an estimator of $\sigma^{2}$ , we would base the confidence interval on the t distribution. Therefore, a $100(1 - \alpha)$ percent confidence interval on the ith treatment mean $\mu_{i}$ is

$$
\overline {{y}} _ {i.} - t _ {\alpha / 2, N - a} \sqrt {\frac {M S _ {E}}{n}} \leq \mu_ {i} \leq \overline {{y}} _ {i.} + t _ {\alpha / 2, N - a} \sqrt {\frac {M S _ {E}}{n}}\tag{3.12}
$$

Differences in treatments are frequently of great practical interest. A $100(1-\alpha)$ percent confidence interval on the difference in any two treatment means, say $\mu_{i}-\mu_{j}$ , would be

$$
\overline {{y}} _ {i.} - \overline {{y}} _ {j.} - t _ {\alpha / 2, N - a} \sqrt {\frac {2 M S _ {E}}{n}} \leq \mu_ {i} - \mu_ {j} \leq \overline {{y}} _ {i.} - \overline {{y}} _ {j.} + t _ {\alpha / 2, N - a} \sqrt {\frac {2 M S _ {E}}{n}}\tag{3.13}
$$

## EXAMPLE 3.3

Using the data in Example 3.1, we may find the estimates of the overall mean and the treatment effects as $\hat{\mu} = 12,355 / 20 = 617.75$ and

Equation 3.12 as

$$
\hat {\tau} _ {1} = \overline {{{y}}} _ {1.} - \overline {{{y}}} _ {..} = 5 5 1. 2 0 - 6 1 7. 7 5 = - 6 6. 5 5
$$

$$
7 0 7. 0 0 - 2. 1 2 0 \sqrt {\frac {3 3 3 . 7 0}{5}} \leq \mu_ {4} \leq 7 0 7. 0 0 + 2. 1 2 0 \sqrt {\frac {3 3 3 . 7 0}{5}}
$$

$$
\hat {\tau} _ {2} = \bar {y} _ {2.} - \bar {y} _ {.} = 5 8 7. 4 0 - 6 1 7. 7 5 = - 3 0. 3 5
$$

or

$$
\hat {\tau} _ {3} = \overline {{{y}}} _ {3.} - \overline {{{y}}} _ {..} = 6 2 5. 4 0 - 6 1 7. 7 5 = 7. 6 5
$$

$$
\hat {\tau} _ {4} = \overline {{{y}}} _ {4.} - \overline {{{y}}} _ {..} = 7 0 7. 0 0 - 6 1 7. 7 5 = 8 9. 2 5
$$

$$
7 0 7. 0 0 - 1 7. 3 2 \leq \mu_ {4} \leq 7 0 7. 0 0 + 1 7. 3 2
$$

A 95 percent confidence interval on the mean of treatment 4 (220 W of RF power) is computed from

Thus, the desired 95 percent confidence interval is $689.68 \leq \mu_{4} \leq 724.32$ .

Simultaneous Confidence Intervals. The confidence interval expressions given in Equations 3.12 and 3.13 are one-at-a-time confidence intervals. That is, the confidence level $1 - \alpha$ applies to only one particular estimate. However, in many problems, the experimenter may wish to calculate several confidence intervals, one for each of a number of means or differences between means. If there are $r$ such $100(1 - \alpha)$ percent confidence intervals of interest, the probability that the $r$ intervals will simultaneously be correct is at least $1 - r\alpha$ . The probability $r\alpha$ is often called the experimentwise error rate or overall confidence coefficient. The number of intervals $r$ does not have to be large before the set of confidence intervals becomes relatively uninformative. For example, if there are $r = 5$ intervals and $\alpha = 0.05$ (a typical choice), the simultaneous confidence level for the set of five confidence intervals is at least 0.75, and if $r = 10$ and $\alpha = 0.05$ , the simultaneous confidence level is at least 0.50.

One approach to ensuring that the simultaneous confidence level is not too small is to replace $\alpha/2$ in the one-at-a-time confidence interval Equations 3.12 and 3.13 with $\alpha/(2r)$ . This is called the Bonferroni method, and it allows the experimenter to construct a set of r simultaneous confidence intervals on treatment means or differences in treatment means for which the overall confidence level is at least $100(1-\alpha)$ percent. When r is not too large, this is a very nice method that leads to reasonably short confidence intervals. For more information, refer to the supplemental text material for Chapter 3.

## 3.3.4 Unbalanced Data

In some single-factor experiments, the number of observations taken within each treatment may be different. We then say that the design is unbalanced. The analysis of variance described may still be used, but slight modifications must be made in the sum of squares formulas. Let $n_{i}$ observations be taken under treatment i ( $i = 1, 2, \ldots, a$ ) and $N = \sum_{i=1}^{a} n_{i}$ . The manual computational formulas for $SS_{T}$ and $SS_{Treatments}$ become

$$
S S _ {T} = \sum_ {i = 1} ^ {a} \sum_ {j = 1} ^ {n _ {i}} y _ {i j} ^ {2} - \frac {y _ {. .} ^ {2}}{N}\tag{3.14}
$$

and

$$
S S _ {\text { Treatments }} = \sum_ {i = 1} ^ {a} \frac {y _ {i .} ^ {2}}{n _ {i}} - \frac {y _ {. .} ^ {2}}{N}\tag{3.15}
$$

No other changes are required in the analysis of variance.

There are two advantages in choosing a balanced design. First, the test statistic is relatively insensitive to small departures from the assumption of equal variances for the a treatments if the sample sizes are equal. This is not the case for unequal sample sizes. Second, the power of the test is maximized if the samples are of equal size.

## 3.4 Model Adequacy Checking

The decomposition of the variability in the observations through an analysis of variance identity (Equation 3.6) is a purely algebraic relationship. However, the use of the partitioning to test formally for no differences in treatment means requires that certain assumptions be satisfied. Specifically, these assumptions are that the observations are adequately described by the model

$$
y _ {i j} = \mu + \tau_ {i} + \epsilon_ {i j}
$$

and that the errors are normally and independently distributed with mean zero and constant but unknown variance $\sigma^{2}$ . If these assumptions are valid, the analysis of variance procedure is an exact test of the hypothesis of no difference in treatment means.

In practice, however, these assumptions will usually not hold exactly. Consequently, it is usually unwise to rely on the analysis of variance until the validity of these assumptions has been checked. Violations of the basic assumptions and model adequacy can be easily investigated by the examination of residuals. We define the residual for observation j in treatment i as

$$
e _ {i j} = y _ {i j} - \hat {y} _ {i j}\tag{3.16}
$$

where $\hat{y}_{ij}$ is an estimate of the corresponding observation $y_{ij}$ obtained as follows:

$$
\begin{array}{r l} & {\hat {y} _ {i j} = \hat {\mu} + \hat {\tau} _ {i}} \\ & {\quad = \overline {{y}} _ {..} + (\overline {{y}} _ {i.} - \overline {{y}} _ {..})} \\ & {\quad = \overline {{y}} _ {i.}} \end{array}\tag{3.17}
$$

Equation 3.17 gives the intuitively appealing result that the estimate of any observation in the ith treatment is just the corresponding treatment average.

Examination of the residuals should be an automatic part of any analysis of variance. If the model is adequate, the residuals should be structureless; that is, they should contain no obvious patterns. Through analysis of residuals, many types of model inadequacies and violations of the underlying assumptions can be discovered. In this section, we show how model diagnostic checking can be done easily by graphical analysis of residuals and how to deal with several commonly occurring abnormalities.

## 3.4.1 The Normality Assumption

A check of the normality assumption could be made by plotting a histogram of the residuals. If the $\mathrm{NID}(0,\sigma^{2})$ assumption on the errors is satisfied, this plot should look like a sample from a normal distribution centered at zero. Unfortunately, with small samples, considerable fluctuation in the shape of a histogram often occurs, so the appearance of a moderate departure from normality does not necessarily imply a serious violation of the assumptions. Gross deviations from normality are potentially serious and require further analysis.

An extremely useful procedure is to construct a normal probability plot of the residuals. Recall from Chapter 2 that we used a normal probability plot of the raw data to check the assumption of normality when using the t-test. In the analysis of variance, it is usually more effective (and straightforward) to do this with the residuals. If the underlying error distribution is normal, this plot will resemble a straight line. In visualizing the straight line, place more emphasis on the central values of the plot than on the extremes.

Table 3.6 shows the original data and the residuals for the etch rate data in Example 3.1. The normal probability plot is shown in Figure 3.4. The general impression from examining this display is that the error distribution is approximately normal. The tendency of the normal probability plot to bend down slightly on the left side and upward slightly on the right side implies that the tails of the error distribution are somewhat thinner than would be anticipated in a normal distribution; that is, the largest residuals are not quite as large (in absolute value) as expected. This plot is not grossly nonnormal, however.

In general, moderate departures from normality are of little concern in the fixed effects analysis of variance (recall our discussion of randomization tests in Section 3.3.2). An error distribution that has considerably thicker or thinner tails than the normal is of more concern than a skewed distribution. Because the F-test is only slightly affected, we say that the analysis of variance (and related procedures such as multiple comparisons) is robust to the normality assumption. Departures from normality usually cause both the true significance level and the power to differ slightly from the advertised values, with the power generally being lower. The random effects model that we will discuss in Section 3.9 and Chapter 13 is more severely affected by nonnormality.

## TABLE 3.6

Etch Rate Data and Residuals from Example 3.1 $^{a}$

<table><tr><td rowspan="2">Power (w)</td><td colspan="8">Observations (j)</td><td rowspan="2"> $\hat{y}_{ij}=\bar{y}_{i}$ .</td></tr><tr><td>1</td><td>2</td><td>3</td><td>4</td><td>5</td><td></td><td></td><td></td></tr><tr><td rowspan="2">160</td><td>23.8</td><td>-9.2</td><td>-21.2</td><td>-12.2</td><td>18.8</td><td></td><td></td><td></td><td></td></tr><tr><td>575 (13)</td><td>542 (14)</td><td>530 (8)</td><td>539 (5)</td><td>570 (4)</td><td>551.2</td><td></td><td></td><td></td></tr><tr><td rowspan="2">180</td><td>-22.4</td><td>5.6</td><td>2.6</td><td>-8.4</td><td>22.6</td><td></td><td></td><td></td><td></td></tr><tr><td>565 (18)</td><td>593 (9)</td><td>590 (6)</td><td>579 (16)</td><td>610 (17)</td><td>587.4</td><td></td><td></td><td></td></tr><tr><td rowspan="2">200</td><td>-25.4</td><td>25.6</td><td>-15.4</td><td>11.6</td><td>3.6</td><td></td><td></td><td></td><td></td></tr><tr><td>600 (7)</td><td>651 (19)</td><td>610 (10)</td><td>637 (20)</td><td>629 (1)</td><td>625.4</td><td></td><td></td><td></td></tr><tr><td rowspan="2">220</td><td>18.0</td><td>-7.0</td><td>8.0</td><td>-22.0</td><td>3.0</td><td></td><td></td><td></td><td></td></tr><tr><td>725 (2)</td><td>700 (3)</td><td>715 (15)</td><td>685 (11)</td><td>710 (12)</td><td>707.0</td><td></td><td></td><td></td></tr></table>

$^{a}$ The residuals are shown in the box in each cell. The numbers in parentheses indicate the order in which each experimental run was made.

■ FIGURE 3.4 Normal probability plot of residuals for Example 3.1

![](images/figure3.4.jpg)

A very common defect that often shows up on normal probability plots is one residual that is very much larger than any of the others. Such a residual is often called an outlier. The presence of one or more outliers can seriously distort the analysis of variance, so when a potential outlier is located, careful investigation is called for. Frequently, the cause of the outlier is a mistake in calculations or a data coding or copying error. If this is not the cause, the experimental circumstances surrounding this run must be carefully studied. If the outlying response is a particularly desirable value (high strength, low cost, etc.), the outlier may be more informative than the rest of the data. We should be careful not to reject or discard an outlying observation unless we have reasonably nonstatistical grounds for doing so. At worst, you may end up with two analyses: one with the outlier and one without.

Several formal statistical procedures may be used for detecting outliers [e.g., see Stefansky (1972), John and Prescott (1975), and Barnett and Lewis (1994)]. Some statistical software packages report the results of a statistical test for normality (such as the Anderson–Darling test) on the normal probability plot of residuals. This should be viewed with caution as those tests usually assume that the data to which they are applied are independent and residuals are not independent.

A rough check for outliers may be made by examining the standardized residuals

$$
d _ {i j} = \frac {e _ {i j}}{\sqrt {M S _ {E}}}\tag{3.18}
$$

If the errors $\epsilon_{ij}$ are $N(0,\sigma^{2})$ , the standardized residuals should be approximately normal with mean zero and unit variance. Thus, about 68 percent of the standardized residuals should fall within the limits $\pm1$ , about 95 percent of them should fall within $\pm2$ , and virtually all of them should fall within $\pm3$ . A residual bigger than 3 or 4 standard deviations from zero is a potential outlier.

For the tensile strength data of Example 3.1, the normal probability plot gives no indication of outliers. Furthermore, the largest standardized residual is

$$
d _ {1} = \frac {e _ {1}}{\sqrt {M S _ {E}}} = \frac {2 5 . 6}{\sqrt {3 3 3 . 7 0}} = \frac {2 5 . 6}{1 8 . 2 7} = 1. 4 0
$$

which should cause no concern.

## 3.4.2 Plot of Residuals in Time Sequence

Plotting the residuals in time order of data collection is helpful in detecting strong correlation between the residuals. A tendency to have runs of positive and negative residuals indicates positive correlation. This would imply that the independence assumption on the errors has been violated. This is a potentially serious problem and one that is difficult to correct, so it is important to prevent the problem if possible when the data are collected. Proper randomization of the experiment is an important step in obtaining independence.

Sometimes the skill of the experimenter (or the subjects) may change as the experiment progresses, or the process being studied may “drift” or become more erratic. This will often result in a change in the error variance over time. This condition often leads to a plot of residuals versus time that exhibits more spread at one end than at the other. Nonconstant variance is a potentially serious problem. We will have more to say on the subject in Sections 3.4.3 and 3.4.4.

Table 3.6 displays the residuals and the time sequence of data collection for the tensile strength data. A plot of these residuals versus run order or time is shown in Figure 3.5. There is no reason to suspect any violation of the independence or constant variance assumptions.

## 3.4.3 Plot of Residuals Versus Fitted Values

If the model is correct and the assumptions are satisfied, the residuals should be structureless; in particular, they should be unrelated to any other variable including the predicted response. A simple check is to plot the residuals versus the fitted values $\hat{y}_{ij}$ . (For the single-factor experiment model, remember that $\hat{y}_{ij} = \overline{y}_{i}$ , the ith treatment average.) This plot should not reveal any obvious pattern. Figure 3.6 plots the residuals versus the fitted values for the tensile strength data of Example 3.1. No unusual structure is apparent.

A defect that occasionally shows up on this plot is nonconstant variance. Sometimes the variance of the observations increases as the magnitude of the observation increases. This would be the case if the error or background noise in the experiment was a constant percentage of the size of the observation. (This commonly happens with many measuring instruments—error is a percentage of the scale reading.) If this were the case, the residuals would get larger as $y_{ij}$ gets larger, and the plot of residuals versus $\hat{y}_{ij}$ would look like an outward-opening funnel or megaphone. Nonconstant variance also arises in cases where the data follow a nonnormal, skewed distribution because in skewed distributions the variance tends to be a function of the mean.

![](images/figure3.5.jpg)  
■ FIGURE 3.5 Plot of residuals versus run order or time

![](images/figure3.6.jpg)  
■ FIGURE 3.6 Plot of residuals versus fitted values

If the assumption of homogeneity of variances is violated, the F-test is only slightly affected in the balanced (equal sample sizes in all treatments) fixed effects model. However, in unbalanced designs or in cases where one variance is very much larger than the others, the problem is more serious. Specifically, if the factor levels having the larger variances also have the smaller sample sizes, the actual type I error rate is larger than anticipated (or confidence intervals have lower actual confidence levels than were specified). Conversely, if the factor levels with larger variances also have the larger sample sizes, the significance levels are smaller than anticipated (confidence levels are higher). This is a good reason for choosing equal sample sizes whenever possible. For random effects models, unequal error variances can significantly disturb inferences on variance components even if balanced designs are used.

Inequality of variance also shows up occasionally on the plot of residuals versus run order. An outward-opening funnel pattern indicates that variability is increasing over time. This could result from operator/subject fatigue, accumulated stress on equipment, changes in material properties such as catalyst degradation, or tool wear, or any of a number of causes.

The usual approach to dealing with nonconstant variance when it occurs for the aforementioned reasons is to apply a variance-stabilizing transformation and then to run the analysis of variance on the transformed data. In this approach, one should note that the conclusions of the analysis of variance apply to the transformed populations.

Considerable research has been devoted to the selection of an appropriate transformation. If experimenters know the theoretical distribution of the observations, they may utilize this information in choosing a transformation. For example, if the observations follow the Poisson distribution, the square root transformation $y_{ij}^{*} = \sqrt{y_{ij}}$ or $y_{ij}^{*} = \sqrt{1 + y_{ij}}$ would be used. If the data follow the lognormal distribution, the logarithmic transformation $y_{ij}^{*} = \log y_{ij}$ is appropriate. For binomial data expressed as fractions, the arcsin transformation $y_{ij}^{*} = \arcsin \sqrt{y_{ij}}$ is useful. When there is no obvious transformation, the experimenter usually empirically seeks a transformation that equalizes the variance regardless of the value of the mean. We offer some guidance on this at the conclusion of this section. In factorial experiments, which we introduce in Chapter 5, another approach is to select a transformation that minimizes the interaction mean square, resulting in an experiment that is easier to interpret. In Chapter 15, we discuss methods for analytically selecting the form of the transformation in more detail. Transformations made for inequality of variance also affect the form of the error distribution. In most cases, the transformation brings the error distribution closer to normal. For more discussion of transformations, refer to Bartlett (1947), Dolby (1963), Box and Cox (1964), and Draper and Hunter (1969).

Statistical Tests for Equality of Variance. Although residual plots are frequently used to diagnose inequality of variance, several statistical tests have also been proposed. These tests may be viewed as formal tests of the hypotheses

$$
\begin{array}{l} H _ {0}: \sigma_ {1} ^ {2} = \sigma_ {2} ^ {2} = \dots = \sigma_ {a} ^ {2} \\ H _ {1}: \text { above   not   true   for   at   least   one } \sigma_ {i} ^ {2} \end{array}
$$

A widely used procedure is Bartlett's test. The procedure involves computing a statistic whose sampling distribution is closely approximated by the chi-square distribution with $a - 1$ degrees of freedom when the $a$ random samples are from independent normal populations. The test statistic is

$$
\chi_ {0} ^ {2} = 2. 3 0 2 6 \frac {q}{c}\tag{3.19}
$$

where

$$
q = (N - a) \mathrm{log} _ {1 0} S _ {p} ^ {2} - \sum_ {i = 1} ^ {a} (n _ {i} - 1) \mathrm{log} _ {1 0} S _ {i} ^ {2}
$$

$$
c = 1 + \frac {1}{3 (a - 1)} \left(\sum_ {i = 1} ^ {a} (n _ {i} - 1) ^ {- 1} - (N - a) ^ {- 1}\right)
$$

$$
S _ {p} ^ {2} = \frac {\sum_ {i = 1} ^ {a} (n _ {i} - 1) S _ {i} ^ {2}}{N - a}
$$

and $S_{i}^{2}$ is the sample variance of the ith population.

The quantity $q$ is large when the sample variances $S_{i}^{2}$ differ greatly and is equal to zero when all $S_{i}^{2}$ are equal. Therefore, we should reject $H_{0}$ on values of $\chi_0^2$ that are too large; that is, we reject $H_{0}$ only when

$$
\chi_ {0} ^ {2} > \chi_ {\alpha , a - 1} ^ {2}
$$

where $\chi_{\alpha,a-1}^{2}$ is the upper $\alpha$ percentage point of the chi-square distribution with a-1 degrees of freedom. The P-value approach to decision making could also be used.

Bartlett's test is very sensitive to the normality assumption. Consequently, when the validity of this assumption is doubtful, Bartlett's test should not be used.

## EXAMPLE 3.4

In the plasma etch experiment, the normality assumption is not in question, so we can apply Bartlett's test to the etch rate data. We first compute the sample variances in each treatment and find that $S_{1}^{2}=400.7$ , $S_{2}^{2}=280.3$ , $S_{3}^{2}=421.3$ , and $S_{4}^{2}=232.5$ . Then

$$
S _ {p} ^ {2} = \frac {4 (4 0 0 . 7) + 4 (2 8 0 . 3) + 4 (4 2 1 . 3) + 4 (2 3 2 . 5)}{1 6} = 3 3 3. 7
$$

$$
\begin{array}{r} q = 1 6 \log_ {1 0} (3 3 3. 7) - 4 [ \log_ {1 0} 4 0 0. 7 + \log_ {1 0} 2 8 0. 3 \\ + \log_ {1 0} 4 2 1. 3 + \log_ {1 0} 2 3 2. 5 ] = 0. 2 1 \end{array}
$$

$$
c = 1 + \frac {1}{3 (3)} \left(\frac {4}{4} - \frac {1}{1 6}\right) = 1. 1 0
$$

and the test statistic is

$$
\chi_ {0} ^ {2} = 2. 3 0 2 6 \frac {(0 . 2 1)}{(1 . 1 0)} = 0. 4 3
$$

From Appendix Table III, we find that $\chi^{2}_{0.05,3}=7.81$ (the P-value is P=0.934), so we cannot reject the null hypothesis. There is no evidence to counter the claim that all five variances are the same. This is the same conclusion reached by analyzing the plot of residuals versus fitted values.

Because Bartlett's test is sensitive to the normality assumption, there may be situations where an alternative procedure would be useful. Anderson and McLean (1974) present a useful discussion of statistical tests for equality of variance. The modified Levene test [see Levene (1960) and Conover, Johnson, and Johnson (1981)] is a very nice procedure that is robust to departures from normality. To test the hypothesis of equal variances in all treatments, the modified Levene test uses the absolute deviation of the observations $y_{ij}$ in each treatment from the treatment median, say, $\tilde{y}_i$ . Denote these deviations by

$$
d _ {i j} = | y _ {i j} - \tilde {y} _ {i}. | \left\{ \begin{array}{l l} i = 1, 2, \ldots , a \\ j = 1, 2, \ldots n _ {i} \end{array} \right.
$$

The modified Levene test then evaluates whether or not the means of these deviations are equal for all treatments. It turns out that if the mean deviations are equal, the variances of the observations in all treatments will be the same. The test statistic for Levene's test is simply the usual ANOVA $F$ -statistic for testing equality of means applied to the absolute deviations.

## EXAMPLE 3.5

A civil engineer is interested in determining whether four different methods of estimating flood flow frequency produce equivalent estimates of peak discharge when applied to the same watershed. Each procedure is used six times on the watershed, and the resulting discharge data (in cubic feet per second) are shown in the upper panel of Table 3.7. The analysis of variance for the data, summarized in Table 3.8, implies that there is a difference in mean peak discharge estimates given by the four procedures. The plot of residuals versus fitted values, shown in Figure 3.7, is disturbing because the outward-opening funnel shape indicates that the constant variance assumption is not satisfied.

We will apply the modified Levene test to the peak discharge data. The upper panel of Table 3.7 contains the treatment medians $\tilde{y}_i$ and the lower panel contains the deviations $d_{ij}$ around the medians. Levene's test consists of conducting a standard analysis of variance on the $d_{ij}$ .

The F-test statistic that results from this is $F_{0} = 4.55$ , for which the P-value is P = 0.0137. Therefore, Levene's test rejects the null hypothesis of equal variances, essentially confirming the diagnosis we made from visual examination of Figure 3.7. The peak discharge data are a good candidate for data transformation.

## TABLE 3.7

Peak Discharge Data

<table><tr><td colspan="2">Estimation Method</td><td colspan="5">Observations</td><td> $\bar{y}_{i.}$ </td><td> $\tilde{y}_{i}$ </td><td> $S_{i}$ </td></tr><tr><td>1</td><td>0.34</td><td>0.12</td><td>1.23</td><td>0.70</td><td>1.75</td><td>0.12</td><td>0.71</td><td>0.520</td><td>0.66</td></tr><tr><td>2</td><td>0.91</td><td>2.94</td><td>2.14</td><td>2.36</td><td>2.86</td><td>4.55</td><td>2.63</td><td>2.610</td><td>1.09</td></tr><tr><td>3</td><td>6.31</td><td>8.37</td><td>9.75</td><td>6.09</td><td>9.82</td><td>7.24</td><td>7.93</td><td>7.805</td><td>1.66</td></tr><tr><td>4</td><td>17.15</td><td>11.82</td><td>10.95</td><td>17.20</td><td>14.35</td><td>16.82</td><td>14.72</td><td>15.59</td><td>2.77</td></tr><tr><td>Estimation Method</td><td colspan="9">Deviations  $d_{ij}$  for the Modified Levene Test</td></tr><tr><td>1</td><td>0.18</td><td>0.40</td><td>0.71</td><td>0.18</td><td>1.23</td><td>0.40</td><td></td><td></td><td></td></tr><tr><td>2</td><td>1.70</td><td>0.33</td><td>0.47</td><td>0.25</td><td>0.25</td><td>1.94</td><td></td><td></td><td></td></tr><tr><td>3</td><td>1.495</td><td>0.565</td><td>1.945</td><td>1.715</td><td>2.015</td><td>0.565</td><td></td><td></td><td></td></tr><tr><td>4</td><td>1.56</td><td>3.77</td><td>4.64</td><td>1.61</td><td>1.24</td><td>1.23</td><td></td><td></td><td></td></tr></table>

TABLE 3.8  
Analysis of Variance for Peak Discharge Data

<table><tr><td>Source of Variation</td><td>Sum of Squares</td><td>Degrees of Freedom</td><td>Mean Square</td><td> $F_0$ </td><td>P-Value</td></tr><tr><td>Methods</td><td>708.3471</td><td>3</td><td>236.1157</td><td>76.07</td><td>&lt; 0.001</td></tr><tr><td>Error</td><td>62.0811</td><td>20</td><td>3.1041</td><td></td><td></td></tr><tr><td>Total</td><td>770.4282</td><td>23</td><td></td><td></td><td></td></tr></table>

![](images/figure3.7.jpg)  
■ FIGURE 3.7 Plot of residuals versus $\hat{y}_{ij}$ for Example 3.5

Empirical Selection of a Transformation. We observed above that if experimenters knew the relationship between the variance of the observations and the mean, they could use this information to guide them in selecting the form of the transformation. We now elaborate on this point and show one method for empirically selecting the form of the required transformation from the data.

Let $E(y) = \mu$ be the mean of $y$ , and suppose that the standard deviation of $y$ is proportional to a power of the mean of $y$ such that

$$
\sigma_ {y} \propto \mu^ {\alpha}
$$

We want to find a transformation on $y$ that yields a constant variance. Suppose that the transformation is a power of the original data, say

$$
y ^ {*} = y ^ {\lambda}\tag{3.20}
$$

Then it can be shown that

$$
\sigma_ {y *} \propto \mu^ {\lambda + \alpha - 1}\tag{3.21}
$$

Clearly, if we set $\lambda = 1 - \alpha$ , the variance of the transformed data $y^{*}$ is constant.

Several of the common transformations discussed previously are summarized in Table 3.9. Note that $\lambda = 0$ implies the log transformation. These transformations are arranged in order of increasing strength. By the strength of a transformation, we mean the amount of curvature it induces. A mild transformation applied to data spanning a narrow range has little effect on the analysis, whereas a strong transformation applied over a large range may have dramatic results. Transformations often have little effect unless the ratio $y_{max}/y_{min}$ is larger than 2 or 3.

In many experimental design situations where there is replication, we can empirically estimate $\alpha$ from the data. Because in the ith treatment combination $\sigma_{y_{i}} \propto \mu_{i}^{\alpha} = \theta \mu_{i}^{\alpha}$ , where $\theta$ is a constant of proportionality, we may take logs to obtain

$$
\log \sigma_ {y _ {i}} = \log \theta + \alpha \log \mu_ {i}\tag{3.22}
$$

Therefore, a plot of $\log \sigma_{y_i}$ versus $\log \mu_i$ would be a straight line with slope $\alpha$ . Because we don't know $\sigma_{y_i}$ and $\mu_i$ , we may substitute reasonable estimates of them in Equation 3.22 and use the slope of the resulting straight line fit as an estimate of $\alpha$ . Typically, we would use the standard deviation $S_i$ and the average $\overline{y}_i$ . of the $i$ th treatment (or, more generally, the $i$ th treatment combination or set of experimental conditions) to estimate $\sigma_{y_i}$ and $\mu_i$ .

To investigate the possibility of using a variance-stabilizing transformation on the peak discharge data from Example 3.5, we plot $\log S_{i}$ versus $\log \overline{y}_{i}$ . in Figure 3.8. The slope of a straight line passing through these four points is close to 1/2 and from Table 3.9 this implies that the square root transformation may be appropriate. The analysis of variance for the transformed data $y^{*} = \sqrt{y}$ is presented in Table 3.10, and a plot of residuals versus the predicted response is shown in Figure 3.9. This residual plot is much improved in comparison to Figure 3.7, so we conclude that the square root transformation has been helpful. Note that in Table 3.10 we have reduced the degrees of freedom for error and total by one to account for the use of the data to estimate the transformation parameter $\alpha$ .

## TABLE 3.9

Variance-Stabilizing Transformations

<table><tr><td>Relationship Between σy and μ</td><td>α</td><td>λ = 1 - α</td><td>Transformation</td><td>Comment</td></tr><tr><td>σy ∝ constant</td><td>0</td><td>1</td><td>No transformation</td><td></td></tr><tr><td>σy ∝ μ1/2</td><td>1/2</td><td>1/2</td><td>Square root</td><td>Poisson (count) data</td></tr><tr><td>σy ∝ μ</td><td>1</td><td>0</td><td>Log</td><td></td></tr><tr><td>σy ∝ μ3/2</td><td>3/2</td><td>-1/2</td><td>Reciprocal square root</td><td></td></tr><tr><td>σy ∝ μ2</td><td>2</td><td>-1</td><td>Reciprocal</td><td></td></tr></table>

![](images/figure3.8.jpg)  
■ FIGURE 3.8 Plot of $\log S_{i}$ versus $\log \bar{y}_{i}$ for the peak discharge data from Example 3.5

![](images/figure3.9.jpg)  
■ FIGURE 3.9 Plot of residuals from transformed data versus $\hat{y}_{ij}^{*}$ for the peak discharge data in Example 3.5

## TABLE 3.10

Analysis of Variance for Transformed Peak Discharge Data, $y^{*} = \sqrt{y}$

<table><tr><td>Source of Variation</td><td>Sum of Squares</td><td>Degrees of Freedom</td><td>Mean Square</td><td> $F_0$ </td><td>P-Value</td></tr><tr><td>Methods</td><td>32.6842</td><td>3</td><td>10.8947</td><td>76.99</td><td>&lt; 0.001</td></tr><tr><td>Error</td><td>2.6884</td><td>19</td><td>0.1415</td><td></td><td></td></tr><tr><td>Total</td><td>35.3726</td><td>22</td><td></td><td></td><td></td></tr></table>

In practice, many experimenters select the form of the transformation by simply trying several alternatives and observing the effect of each transformation on the plot of residuals versus the predicted response. The transformation that produced the most satisfactory residual plot is then selected. Alternatively, there is a formal method called the Box-Cox Method for selecting a variance-stability transformation. In Chapter 15, we discuss and illustrate this procedure. It is widely used and implemented in many software packages.

## 3.4.4 Plots of Residuals Versus Other Variables

If data have been collected on any other variables that might possibly affect the response, the residuals should be plotted against these variables. For example, in the tensile strength experiment of Example 3.1, strength may be significantly affected by the thickness of the fiber, so the residuals should be plotted versus fiber thickness. If different testing machines were used to collect the data, the residuals should be plotted against machines. Patterns in such residual plots imply that the variable affects the response. This suggests that the variable should be either controlled more carefully in future experiments or included in the analysis.

## 3.5 Practical Interpretation of Results

After conducting the experiment, performing the statistical analysis, and investigating the underlying assumptions, the experimenter is ready to draw practical conclusions about the problem he or she is studying. Often this is relatively easy, and certainly in the simple experiments we have considered so far, this might be done somewhat informally, perhaps by inspection of graphical displays such as the box plots and scatter diagram in Figures 3.1 and 3.2. However, in some cases, more formal techniques need to be applied. We present some of these techniques in this section.

## 3.5.1 A Regression Model

The factors involved in an experiment can be either quantitative or qualitative. A quantitative factor is one whose levels can be associated with points on a numerical scale, such as temperature, pressure, or time. Qualitative factors, on the other hand, are factors for which the levels cannot be arranged in order of magnitude. Operators, batches of raw material, and shifts are typical qualitative factors because there is no reason to rank them in any particular numerical order.

Insofar as the initial design and analysis of the experiment are concerned, both types of factors are treated identically. The experimenter is interested in determining the differences, if any, between the levels of the factors. In fact, the analysis of variance treats the design factor as if it were qualitative or categorical. If the factor is really qualitative, such as operators, it is meaningless to consider the response for a subsequent run at an intermediate level of the factor. However, with a quantitative factor such as time, the experimenter is usually interested in the entire range of values used, particularly the response from a subsequent run at an intermediate factor level. That is, if the levels 1.0, 2.0, and 3.0 hours are used in the experiment, we may wish to predict the response at 2.5 hours. Thus, the experimenter is frequently interested in developing an interpolation equation for the response variable in the experiment. This equation is an empirical model of the process that has been studied.

The general approach to fitting empirical models is called regression analysis, which is discussed extensively in Chapter 10. See also the supplemental text material for this chapter. This section briefly illustrates the technique using the etch rate data of Example 3.1.

Figure 3.10 presents scatter diagrams of etch rate y versus the power x for the experiment in Example 3.1. From examining the scatter diagram, it is clear that there is a strong relationship between the etch rate and power. As a first approximation, we could try fitting a linear model to the data, say

$$
y = \beta_ {0} + \beta_ {1} x + \epsilon
$$

where $\beta_{0}$ and $\beta_{1}$ are unknown parameters to be estimated and $\epsilon$ is a random error term. The method often used to estimate the parameters in a model such as this is the method of least squares. This consists of choosing estimates of the $\beta$ 's such that the sum of the squares of the errors (the $\epsilon$ 's) is minimized. The least squares fit in our example is

$$
\hat {y} = 1 3 7. 6 2 + 2. 5 2 7 x
$$

(If you are unfamiliar with regression methods, see Chapter 10 and the supplemental text material for this chapter.)

This linear model is shown in Figure 3.10a. It does not appear to be very satisfactory at the higher power settings. Perhaps an improvement can be obtained by adding a quadratic term in x. The resulting quadratic model fit is

$$
\hat {y} = 1 1 4 7. 7 7 - 8. 2 5 5 5 x + 0. 0 2 8 3 7 5 x ^ {2}
$$

This quadratic fit is shown in Figure 3.10b. The quadratic model appears to be superior to the linear model because it provides a better fit at the higher power settings.

In general, we would like to fit the lowest order polynomial that adequately describes the system or process. In this example, the quadratic polynomial seems to fit better than the linear model, so the extra complexity of the quadratic model is justified. Selecting the order of the approximating polynomial is not always easy, however, and it is relatively easy to overfit, that is, to add high-order polynomial terms that do not really improve the fit but increase the complexity of the model and often damage its usefulness as a predictor or interpolation equation.

In this example, the empirical model could be used to predict etch rate at power settings within the region of experimentation. In other cases, the empirical model could be used for process optimization, that is, finding the levels of the design variables that result in the best values of the response. We will discuss and illustrate these problems extensively later in the book.

![](images/figure3.10.jpg)  

■ FIGURE 3.10 Scatter diagrams and regression models for the etch rate data of Example 3.1

## 3.5.2 Comparisons Among Treatment Means

Suppose that in conducting an analysis of variance for the fixed effects model the null hypothesis is rejected. Thus, there are differences between the treatment means but exactly which means differ is not specified. Sometimes in this situation, further comparisons and analysis among groups of treatment means may be useful. The ith treatment mean is defined as $\mu_{i} = \mu + \tau_{i}$ , and $\mu_{i}$ is estimated by $\overline{y}_{i}$ . Comparisons between treatment means are made in terms of either the treatment totals $\{y_{i}\}$ or the treatment averages $\{\overline{y}_{i}\}$ . The procedures for making these comparisons are usually called multiple comparison methods. In the next several sections, we discuss methods for making comparisons among individual treatment means or groups of these means.

## 3.5.3 Graphical Comparisons of Means

It is very easy to develop a graphical procedure for the comparison of means following an analysis of variance. Suppose that the factor of interest has $a$ levels and that $\overline{y}_{1},\overline{y}_{2},\ldots ,\overline{y}_{a}$ , are the treatment averages. If we know $\sigma$ , any treatment average would have a standard deviation $\sigma /\sqrt{n}$ . Consequently, if all factor level means are identical, the observed sample means $\overline{y}_i$ . would behave as if they were a set of observations drawn at random from a normal distribution with mean $\overline{y}_i$ and standard deviation $\sigma /\sqrt{n}$ . Visualize a normal distribution capable of being slid along an axis below which the $\overline{y}_{1},\overline{y}_{2},\ldots ,\overline{y}_{a}$ , are plotted. If the treatment means are all equal, there should be some position for this distribution that makes it obvious that the $\overline{y}_i$ . values were drawn from the same distribution. If this is not the case, the $\overline{y}_i$ . values that appear not to have been drawn from this distribution are associated with factor levels that produce different mean responses.

The only flaw in this logic is that $\sigma$ is unknown. Box, Hunter, and Hunter (2005) point out that we can replace $\sigma$ with $\sqrt{MS_E}$ from the analysis of variance and use a $t$ distribution with a scale factor $\sqrt{MS_E / n}$ instead of the normal. Such an arrangement for the etch rate data of Example 3.1 is shown in Figure 3.11. Focus on the $t$ distribution shown as a solid line curve in the middle of the display.

To sketch the t distribution in Figure 3.11, simply multiply the abscissa t value by the scale factor

$$
\sqrt {M S _ {E} / n} = \sqrt {3 3 0 . 7 0 / 5} = 8. 1 3
$$

![](images/figure3.11.jpg)  
■ FIGURE 3.11 Etch rate averages from Example 3.1 in relation to a t distribution with scale factor $\sqrt{MS_{E}/n} = \sqrt{330.70/5} = 8.13$

and plot this against the ordinate of t at that point. Because the t distribution looks much like the normal, except that it is a little flatter near the center and has longer tails, this sketch is usually easily constructed by eye. If you wish to be more precise, there is a table of abscissa t values and the corresponding ordinates in Box, Hunter, and Hunter (2005). The distribution can have an arbitrary origin, although it is usually best to choose one in the region of the $\overline{y}_{i}$ values to be compared. In Figure 3.11, the origin is 615 Å/min.

Now visualize sliding the t distribution in Figure 3.11 along the horizontal axis as indicated by the dashed lines and examine the four means plotted in the figure. Notice that there is no location for the distribution such that all four averages could be thought of as typical, randomly selected observations from the distribution. This implies that all four means are not equal; thus, the figure is a graphical display of the ANOVA results. Furthermore, the figure indicates that all four levels of power (160, 180, 200, 220 W) produce mean etch rates that differ from each other. In other words, $\mu_{1} \neq \mu_{2} \neq \mu_{3} \neq \mu_{4}$ .

This simple procedure is a rough but effective technique for many multiple comparison problems. However, there are more formal methods. We now give a brief discussion of some of these procedures.

## 3.5.4 Contrasts

Many multiple comparison methods use the idea of a contrast. Consider the plasma etching experiment of Example 3.1. Because the null hypothesis was rejected, we know that some power settings produce different etch rates than others, but which ones actually cause this difference? We might suspect at the outset of the experiment that 200 W and 220 W produce the same etch rate, implying that we would like to test the hypothesis

$$
\begin{array}{l} {H _ {0}: \mu_ {3} = \mu_ {4}} \\ {H _ {1}: \mu_ {3} \neq \mu_ {4}} \end{array}
$$

or equivalently

$$
\begin{array}{l} H _ {0}: \mu_ {3} - \mu_ {4} = 0 \\ H _ {1}: \mu_ {3} - \mu_ {4} \neq 0 \end{array}\tag{3.23}
$$

If we had suspected at the start of the experiment that the average of the lowest levels of power did not differ from the average of the highest levels of power, then the hypothesis would have been

$$
\begin{array}{l} {H _ {0}: \mu_ {1} + \mu_ {2} = \mu_ {3} + \mu_ {4}} \\ {H _ {1}: \mu_ {1} + \mu_ {2} \neq \mu_ {3} + \mu_ {4}} \end{array}
$$

or

$$
\begin{array}{l} H _ {0}: \mu_ {1} + \mu_ {2} - \mu_ {3} - \mu_ {4} = 0 \\ H _ {1}: \mu_ {1} + \mu_ {2} - \mu_ {3} - \mu_ {4} \neq 0 \end{array}\tag{3.24}
$$

In general, a contrast is a linear combination of parameters of the form

$$
\Gamma = \sum_ {i = 1} ^ {a} c _ {i} \mu_ {i}
$$

where the contrast constants $c_{1}, c_{2}, \ldots, c_{a}$ sum to zero; that is, $\sum_{i=1}^{a} c_{i} = 0$ . Both of the above hypotheses can be expressed in terms of contrasts:

$$
H _ {0}: \sum_ {i = 1} ^ {a} c _ {i} \mu_ {i} = 0
$$

$$
H _ {1}: \sum_ {i = 1} ^ {a} c _ {i} \mu_ {i} \neq 0\tag{3.25}
$$

The contrast constants for the hypotheses in Equation 3.23 are $c_{1}=c_{2}=0$ , $c_{3}=+1$ , and $c_{4}=-1$ , whereas for the hypotheses in Equation 3.24, they are $c_{1}=c_{2}=+1$ and $c_{3}=c_{4}=-1$ .

Testing hypotheses involving contrasts can be done in two basic ways. The first method uses a t-test. Write the contrast of interest in terms of the treatment averages, giving

$$
C = \sum_ {i = 1} ^ {a} c _ {i} \overline {{y}} _ {i.}
$$

The variance of $C$ is

$$
V (C) = \frac {\sigma^ {2}}{n} \sum_ {i = 1} ^ {a} c _ {i} ^ {2}\tag{3.26}
$$

when the sample sizes in each treatment are equal. If the null hypothesis in Equation 3.25 is true, the ratio

$$
\frac {\sum_ {i = 1} ^ {a} c _ {i} \overline {{y}} _ {i .}}{\sqrt {\frac {\sigma^ {2}}{n} \sum_ {i = 1} ^ {a} c _ {i} ^ {2}}}
$$

has the $N(0,1)$ distribution. Now we would replace the unknown variance $\sigma^{2}$ by its estimate, the mean square error $MS_{E}$ and use the statistic

$$
t _ {0} = \frac {\sum_ {i = 1} ^ {a} c _ {i} \overline {{y}} _ {i .}}{\sqrt {\frac {M S _ {E}}{n} \sum_ {i = 1} ^ {a} c _ {i} ^ {2}}}\tag{3.27}
$$

to test the hypotheses in Equation 3.25. The null hypothesis would be rejected if $|t_0|$ in Equation 3.27 exceeds $t_{\alpha / 2,N - a}$ .

The second approach uses an F-test. Now the square of a t random variable with v degrees of freedom is an F random variable with 1 numerator and v denominator degrees of freedom. Therefore, we can obtain

$$
F _ {0} = t _ {0} ^ {2} = \frac {\left(\sum_ {i = 1} ^ {a} c _ {i} \overline {{y}} _ {i .}\right) ^ {2}}{\frac {M S _ {E}}{n} \sum_ {i = 1} ^ {a} c _ {i} ^ {2}}\tag{3.28}
$$

$$
F _ {0} = \frac {M S _ {C}}{M S _ {E}} = \frac {S S _ {C} / 1}{M S _ {E}}
$$

as an $F$ -statistic for testing Equation 3.25. The null hypothesis would be rejected if $F_0 > F_{\alpha,1,N - a}$ . We can write the test statistic of Equation 3.28 as

where the single-degree-of-freedom contrast sum of squares is

$$
S S _ {C} = \frac {\left(\sum_ {i = 1} ^ {a} c _ {i} \overline {{y}} _ {i .}\right) ^ {2}}{\frac {1}{n} \sum_ {i = 1} ^ {a} c _ {i} ^ {2}}\tag{3.29}
$$

Confidence Interval for a Contrast. Instead of testing hypotheses about a contrast, it may be more useful to construct a confidence interval. Suppose that the contrast of interest is

$$
\Gamma = \sum_ {i = 1} ^ {a} c _ {i} \mu_ {i}
$$

Replacing the treatment means with the treatment averages yields

$$
C = \sum_ {i = 1} ^ {a} c _ {i} \overline {{y}} _ {i.}
$$

Because

$$
E \left(\sum_ {i = 1} ^ {a} c _ {i} \bar {y} _ {i.}\right) = \sum_ {i = 1} ^ {a} c _ {i} \mu_ {i} \quad \text { and } \quad V (C) = \sigma^ {2} / n \sum_ {i = 1} ^ {a} c _ {i} ^ {2}
$$

the $100(1 - \alpha)$ percent confidence interval on the contrast $\Sigma_{i=1}^{a} c_i \mu_i$ is

$$
\sum_ {i = 1} ^ {a} c _ {i} \overline {{y}} _ {i.} - t _ {\alpha / 2, N - a} \sqrt {\frac {M S _ {E}}{n} \sum_ {i = 1} ^ {a} c _ {i} ^ {2}} \leq \sum_ {i = 1} ^ {a} c _ {i} \mu_ {i} \leq \sum_ {i = 1} ^ {a} c _ {i} \overline {{y}} _ {i.} + t _ {\alpha / 2, N - a} \sqrt {\frac {M S _ {E}}{n} \sum_ {i = 1} ^ {a} c _ {i} ^ {2}}\tag{3.30}
$$

Note that we have used $MS_{E}$ to estimate $\sigma^{2}$ . Clearly, if the confidence interval in Equation 3.30 includes zero, we would be unable to reject the null hypothesis in Equation 3.25.

Standardized Contrast. When more than one contrast is of interest, it is often useful to evaluate them on the same scale. One way to do this is to standardize the contrast so that it has variance $\sigma^{2}$ . If the contrast $\Sigma_{i=1}^{a}c_{i}\mu_{i}$ is written in terms of treatment averages as $\Sigma_{i=1}^{a}c_{i}\overline{y}_{i}$ , dividing it by $\sqrt{(1/n)\Sigma_{i=1}^{a}c_{i}^{2}}$ will produce a standardized contrast with variance $\sigma^{2}$ . Effectively, then, the standardized contrast is

$$
\sum_ {i = 1} ^ {a} c _ {i} ^ {*} \overline {{y}} _ {i}.
$$

where

$$
c _ {i} ^ {*} = \frac {c _ {i}}{\sqrt {\frac {1}{n} \sum_ {i = 1} ^ {a} c _ {i} ^ {2}}}
$$

Unequal Sample Sizes. When the sample sizes in each treatment are different, minor modifications are made in the above results. First, note that the definition of a contrast now requires that

$$
\sum_ {i = 1} ^ {a} n _ {i} c _ {i} = 0
$$

Other required changes are straightforward. For example, the t statistic in Equation 3.27 becomes

$$
t _ {0} = \frac {\sum_ {i = 1} ^ {a} c _ {i} \overline {{y}} _ {i .}}{\sqrt {M S _ {E} \sum_ {i = 1} ^ {a} \frac {c _ {i} ^ {2}}{n _ {i}}}}
$$

and the contrast sum of squares from Equation 3.29 becomes

$$
S S _ {C} = \frac {\left(\sum_ {i = 1} ^ {a} c _ {i} \overline {{y}} _ {i .}\right) ^ {2}}{\sum_ {i = 1} ^ {a} \frac {c _ {i} ^ {2}}{n _ {i}}}
$$

## 3.5.5 Orthogonal Contrasts

A useful special case of the procedure in Section 3.5.4 is that of orthogonal contrasts. Two contrasts with coefficients $\{c_i\}$ and $\{d_i\}$ are orthogonal if

$$
\sum_ {i = 1} ^ {a} c _ {i} d _ {i} = 0
$$

or, for an unbalanced design, if

$$
\sum_ {i = 1} ^ {a} c _ {i} d _ {i} / n _ {i} = 0
$$

For a treatments, the set of a - 1 orthogonal contrasts partition the sum of squares due to treatments into a - 1 independent single-degree-of-freedom components. Thus, tests performed on orthogonal contrasts are independent.

There are many ways to choose the orthogonal contrast coefficients for a set of treatments. Usually, something in the nature of the experiment should suggest which comparisons will be of interest. For example, if there are a = 3 treatments, with treatment 1 a control and treatments 2 and 3 actual levels of the factor of interest to the experimenter, appropriate orthogonal contrasts might be as follows:

<table><tr><td>Treatment</td><td colspan="2">Coefficients for Orthogonal Contrasts</td></tr><tr><td>1 (control)</td><td>-2</td><td>0</td></tr><tr><td>2 (level 1)</td><td>1</td><td>-1</td></tr><tr><td>3 (level 2)</td><td>1</td><td>1</td></tr></table>

Note that contrast 1 with $c_{i} = -2, 1, 1$ compares the average effect of the factor with the control, whereas contrast 2 with $d_{i} = 0, -1, 1$ compares the two levels of the factor of interest.

Generally, the method of contrasts (or orthogonal contrasts) is useful for what are called preplanned comparisons. That is, the contrasts are specified prior to running the experiment and examining the data. The reason for this is that if comparisons are selected after examining the data, most experimenters would construct tests that correspond to large observed differences in means. These large differences could be the result of the presence of real effects, or they could be the result of random error. If experimenters consistently pick the largest differences to compare, they will inflate the type I error of the test because it is likely that, in an unusually high percentage of the comparisons selected, the observed differences will be the result of error. Examining the data to select comparisons of potential interest is often called data snooping. The Scheffé method for all comparisons, discussed in the next section, permits data snooping.

## EXAMPLE 3.6

Consider the plasma etching experiment in Example 3.1. There are four treatment means and three degrees of freedom between these treatments. Suppose that prior to running the experiment the following set of comparisons among the treatment means (and their associated contrasts) were specified:

<table><tr><td>Hypothesis</td><td>Contrast</td></tr><tr><td> $H_0: \mu_1 = \mu_2$ </td><td> $C_1 = \overline{y}_{1.} - \overline{y}_{2.}$ </td></tr><tr><td> $H_0: \mu_1 + \mu_2 = \mu_3 + \mu_4$ </td><td> $C_2 = \overline{y}_{1.} + \overline{y}_{2.} - \overline{y}_{3.} - \overline{y}_{4.}$ </td></tr><tr><td> $H_0: \mu_3 = \mu_4$ </td><td> $C_3 = \overline{y}_{3.} - \overline{y}_{4.}$ </td></tr></table>

Notice that the contrast coefficients are orthogonal. Using the data in Table 3.4, we find the numerical values of the contrasts and the sums of squares to be as follows:

$$
C _ {1} = + 1 (5 5 1. 2) - 1 (5 8 7. 4) = - 3 6. 2
$$

$$
S S _ {C _ {1}} = \frac {(- 3 6 . 2) ^ {2}}{\frac {1}{5} (2)} = 3 2 7 6. 1 0
$$

$$
C _ {2} = \begin{array}{l} + 1 (5 5 1. 2) + 1 (5 8 7. 4) \\ - 1 (6 2 5. 4) - 1 (7 0 7. 0) \end{array} = - 1 9 3. 8
$$

$$
S S _ {C _ {2}} = \frac {(- 1 9 3 . 8) ^ {2}}{\frac {1}{5} (4)} = 4 6, 9 4 8. 0 5
$$

$$
C _ {3} = + 1 (6 2 5. 4) - 1 (7 0 7. 6) = - 8 1. 6
$$

$$
S S _ {C _ {3}} = \frac {(- 8 1 . 6) ^ {2}}{\frac {1}{5} (2)} = 1 6, 6 4 6. 4 0
$$

These contrast sums of squares completely partition the treatment sum of squares. The tests on such orthogonal contrasts are usually incorporated in the ANOVA, as shown in Table 3.11. We conclude from the P-values that there are significant differences in mean etch rates between levels 1 and 2 and between levels 3 and 4 of the power settings, and that the average of levels 1 and 2 does differ significantly from the average of levels 3 and 4 at the $\alpha = 0.05$ level.

## TABLE 3.11

Analysis of Variance for the Plasma Etching Experiment

<table><tr><td>Source of Variation</td><td>Sum of Squares</td><td>Degrees of Freedom</td><td>Mean Square</td><td> $F_0$ </td><td>P-Value</td></tr><tr><td>Power setting</td><td>66,870.55</td><td>3</td><td>22,290.18</td><td>66.80</td><td>&lt; 0.001</td></tr><tr><td colspan="6">Orthogonal contrasts</td></tr><tr><td> $C_1: \mu_1 = \mu_2$ </td><td>(3276.10)</td><td>1</td><td>3276.10</td><td>9.82</td><td>&lt; 0.01</td></tr><tr><td> $C_2: \mu_1 + \mu_3 = \mu_3 + \mu_4$ </td><td>(46,948.05)</td><td>1</td><td>46,948.05</td><td>140.69</td><td>&lt; 0.001</td></tr><tr><td> $C_3: \mu_3 = \mu_4$ </td><td>(16,646.40)</td><td>1</td><td>16,646.40</td><td>49.88</td><td>&lt; 0.001</td></tr><tr><td>Error</td><td>5,339.20</td><td>16</td><td>333.70</td><td></td><td></td></tr><tr><td>Total</td><td>72,209.75</td><td>19</td><td></td><td></td><td></td></tr></table>

## 3.5.6 Scheffé's Method for Comparing All Contrasts

In many situations, experimenters may not know in advance which contrasts they wish to compare, or they may be interested in more than a - 1 possible comparisons. In many exploratory experiments, the comparisons of interest are discovered only after preliminary examination of the data. Scheffé (1953) has proposed a method for comparing any and all possible contrasts between treatment means. In the Scheffé method, the type I error is at most $\alpha$ for any of the possible comparisons.

Suppose that a set of m contrasts in the treatment means

$$
\Gamma_ {u} = c _ {1 u} \mu_ {1} + c _ {2 u} \mu_ {2} + \dots + c _ {a u} \mu_ {a} \quad u = 1, 2, \ldots , m\tag{3.31}
$$

of interest have been determined. The corresponding contrast in the treatment averages $\overline{y}_{i}$ is

$$
C _ {u} = c _ {1 u} \overline {{y}} _ {1.} + c _ {2 u} \overline {{y}} _ {2.} + \dots + c _ {a u} \overline {{y}} _ {a.} \quad u = 1, 2, \ldots , m\tag{3.32}
$$

and the standard error of this contrast is

$$
S _ {C _ {u}} = \sqrt {M S _ {E} \sum_ {i = 1} ^ {a} (c _ {i u} ^ {2} / n _ {i})}\tag{3.33}
$$

where $n_{i}$ is the number of observations in the ith treatment. It can be shown that the critical value against which $C_{u}$ should be compared is

$$
S _ {\alpha , u} = S _ {C _ {u}} \sqrt {(a - 1) F _ {\alpha , a - 1 , N - a}}\tag{3.34}
$$

To test the hypothesis that the contrast $\Gamma_{u}$ differs significantly from zero, refer $C_{u}$ to the critical value. If $|C_{u}| > S_{\alpha,u}$ , the hypothesis that the contrast $\Gamma_{u}$ equals zero is rejected.

The Scheffé procedure can also be used to form confidence intervals for all possible contrasts among treatment means. The resulting intervals, say $C_{u} - S_{\alpha,u} \leq \Gamma_{u} \leq C_{u} + S_{\alpha,u}$ , are simultaneous confidence intervals in that the probability that all of them are simultaneously true is at least $1 - \alpha$ .

To illustrate the procedure, consider the data in Example 3.1 and suppose that the contrasts of interests are

$$
\Gamma_ {1} = \mu_ {1} + \mu_ {2} - \mu_ {3} - \mu_ {4}
$$

and

$$
\Gamma_ {2} = \mu_ {1} - \mu_ {4}
$$

The numerical values of these contrasts are

$$
\begin{array}{r l} C _ {1} & = \overline {{y}} _ {1.} + \overline {{y}} _ {2.} - \overline {{y}} _ {3.} - \overline {{y}} _ {4.} \\ & = 5 5 1. 2 + 5 8 7. 4 - 6 2 5. 4 - 7 0 7. 0 = - 1 9 3. 8 0 \end{array}
$$

and

$$
\begin{array}{r l} C _ {2} & = \overline {{y}} _ {1.} - \overline {{y}} _ {4.} \\ & = 5 5 1. 2 - 7 0 7. 0 = - 1 5 5. 8 \end{array}
$$

The standard errors are found from Equation 3.33 as

$$
S _ {C _ {1}} = \sqrt {M S _ {E} \sum_ {i = 1} ^ {5} (c _ {i 1} ^ {2} / n _ {i})} = \sqrt {3 3 3 . 7 0 (1 + 1 + 1 + 1) / 5} = 1 6. 3 4
$$

and

$$
S _ {C _ {2}} = \sqrt {M S _ {E} \sum_ {i = 1} ^ {5} (c _ {i 2} ^ {2} / n _ {i})} = \sqrt {3 3 3 . 7 0 (1 + 1) / 5} = 1 1. 5 5
$$

From Equation 3.34, the 1 percent critical values are

$$
S _ {0. 0 1, 1} = S _ {C _ {1}} \sqrt {(a - 1) F _ {0 . 0 1 , a - 1 , N - a}} = 1 6. 3 4 \sqrt {3 (5 . 2 9)} = 6 5. 0 9
$$

and

$$
S _ {0. 0 1, 2} = S _ {C _ {2}} \sqrt {(a - 1) F _ {0 . 0 1 , a - 1 , N - a}} = 1 1. 5 5 \sqrt {3 (5 . 2 9)} = 4 5. 9 7
$$

Because $|C_{1}| > S_{0.01,1}$ , we conclude that the contrast $\Gamma_{1} = \mu_{1} + \mu_{2} - \mu_{3} - \mu_{4}$ does not equal zero; that is, we conclude that the mean etch rates of power settings 1 and 2 as a group differ from the means of power settings 3 and 4 as a group. Furthermore, because $|C_{2}| > S_{0.01,2}$ , we conclude that the contrast $\Gamma_{2} = \mu_{1} - \mu_{4}$ does not equal zero; that is, the mean etch rates of treatments 1 and 4 differ significantly.

## 3.5.7 Comparing Pairs of Treatment Means

In many practical situations, we will wish to compare only pairs of means. Frequently, we can determine which means differ by testing the differences between all pairs of treatment means. Thus, we are interested in contrasts of the form $\Gamma = \mu_{j} - \mu_{j}$ for all $i \neq j$ . Although the Scheffé method described in the previous section could be easily applied to this problem, it is not the most sensitive procedure for such comparisons. We now turn to a consideration of methods specifically designed for pairwise comparisons between all a population means.

Suppose that we are interested in comparing all pairs of a treatment means and that the null hypotheses that we wish to test are $H_{0}: \mu_{i} = \mu_{j}$ for all $i \neq j$ . There are numerous procedures available for this problem. We now present two popular methods for making such comparisons.

Tukey's Test. Suppose that, following an ANOVA in which we have rejected the null hypothesis of equal treatment means, we wish to test all pairwise mean comparisons:

$$
\begin{array}{l} H _ {0}: \mu_ {i} = \mu_ {j} \\ H _ {1}: \mu_ {i} \neq \mu_ {j} \end{array}
$$

for all $i \neq j$ . Tukey (1953) proposed a procedure for testing hypotheses for which the overall significance level is exactly $\alpha$ when the sample sizes are equal and at most $\alpha$ when the sample sizes are unequal. His procedure can also be used to construct confidence intervals on the differences in all pairs of means. For these intervals, the simultaneous confidence level is $100(1 - \alpha)$ percent when the sample sizes are equal and at least $100(1 - \alpha)$ percent when sample sizes are unequal. In other words, the Tukey procedure controls the experimentwise or “family” error rate at the selected level $\alpha$ . This is an excellent data snooping procedure when interest focuses on pairs of means.

Tukey's procedure makes use of the distribution of the studentized range statistic

$$
q = \frac {\overline {{{y}}} _ {\max} - \overline {{{y}}} _ {\min}}{\sqrt {M S _ {E} / n}}
$$

where $\overline{y}_{\mathrm{max}}$ and $\overline{y}_{\mathrm{min}}$ are the largest and smallest sample means, respectively, out of a group of $p$ sample means. Appendix Table V contains values of $q_{\alpha}(p,f)$ , the upper $\alpha$ percentage points of $q$ , where $f$ is the number of degrees of freedom associated with the $MS_E$ . For equal sample sizes, Tukey's test declares two means significantly different if the absolute value of their sample differences exceeds

$$
T _ {\alpha} = q _ {\alpha} (a, f) \sqrt {\frac {M S _ {E}}{n}}\tag{3.35}
$$

Equivalently, we could construct a set of $100(1 - \alpha)$ percent confidence intervals for all pairs of means as follows:

$$
\begin{array}{r l} & {\overline {{y}} _ {i.} - \overline {{y}} _ {j.} - q _ {\alpha} (a, f) \sqrt {\frac {M S _ {E}}{n}} \leq \mu_ {i} - \mu_ {j}} \\ & {\qquad \leq \overline {{y}} _ {i.} - \overline {{y}} _ {j.} + q _ {\alpha} (a, f) \sqrt {\frac {M S _ {E}}{n}}, i \neq j.} \end{array}\tag{3.36}
$$

When sample sizes are not equal, Equations 3.35 and 3.36 become

$$
T _ {\alpha} = \frac {q _ {\alpha} (a , f)}{\sqrt {2}} \sqrt {M S _ {E} \left(\frac {1}{n _ {i}} + \frac {1}{n _ {j}}\right)}\tag{3.37}
$$

and

$$
\begin{array}{r l} & {\overline {{y}} _ {i.} - \overline {{y}} _ {j.} - \frac {q _ {\alpha} (a , f)}{\sqrt {2}} \sqrt {M S _ {E} \left(\frac {1}{n _ {i}} + \frac {1}{n _ {j}}\right)} \leq \mu_ {i} - \mu_ {j}} \\ & {\qquad \leq \overline {{y}} _ {i.} - \overline {{y}} _ {j.} + \frac {q _ {\alpha} (a , f)}{\sqrt {2}} \sqrt {M S _ {E} \left(\frac {1}{n _ {i}} + \frac {1}{n _ {j}}\right)}, i \neq j} \end{array}\tag{3.38}
$$

respectively. The unequal sample size version is sometimes called the Tukey–Kramer procedure.

## EXAMPLE 3.7

To illustrate Tukey's test, we use the data from the plasma etching experiment in Example 3.1. With $\alpha = 0.05$ and $f = 16$ degrees of freedom for error, Appendix Table V gives $q_{0.05}(4,16) = 4.05$ . Therefore, from Equation 3.35,

$$
T _ {0. 0 5} = q _ {0. 0 5} (4, 1 6) \sqrt {\frac {M S _ {E}}{n}} = 4. 0 5 \sqrt {\frac {3 3 3 . 7 0}{5}} = 3 3. 0 9
$$

Thus, any pairs of treatment averages that differ in absolute value by more than 33.09 would imply that the corresponding pair of population means are significantly different. The four treatment averages are

$$
\begin{array}{r l} \overline {{y}} _ {1.} = 5 5 1. 2 & \overline {{y}} _ {2.} = 5 8 7. 4 \\ \overline {{y}} _ {3.} = 6 2 5. 4 & \overline {{y}} _ {4.} = 7 0 7. 0 \end{array}
$$

and the differences in averages are

$$
\overline {{{y}}} _ {1.} - \overline {{{y}}} _ {2.} = 5 5 1. 2 - 5 8 7. 4 = - 3 6. 2 0 ^ {*}
$$

$$
\overline {{{y}}} _ {1.} - \overline {{{y}}} _ {3.} = 5 5 1. 2 - 6 2 5. 4 = - 7 4. 2 0 ^ {*}
$$

$$
\overline {{{y}}} _ {1.} - \overline {{{y}}} _ {4.} = 5 5 1. 2 - 7 0 7. 0 = - 1 5 5. 8 ^ {*}
$$

$$
\overline {{{y}}} _ {2.} - \overline {{{y}}} _ {3.} = 5 8 7. 4 - 6 2 5. 4 = - 3 8. 0 ^ {*}
$$

$$
\overline {{{y}}} _ {2.} - \overline {{{y}}} _ {4.} = 5 8 7. 4 - 7 0 7. 0 = - 1 1 9. 6 ^ {*}
$$

$$
\overline {{y}} _ {3.} - \overline {{y}} _ {4.} = 6 2 5. 4 - 7 0 7. 0 = - 8 1. 6 0 ^ {*}
$$

The starred values indicate the pairs of means that are significantly different. Note that the Tukey procedure indicates that all pairs of means differ. Therefore, each power setting results in a mean etch rate that differs from the mean etch rate at any other power setting.

When using any procedure for pairwise testing of means, we occasionally find that the overall F-test from the ANOVA is significant, but the pairwise comparison of means fails to reveal any significant differences. This situation occurs because the F-test is simultaneously considering all possible contrasts involving the treatment means, not just pairwise comparisons. That is, in the data at hand, the significant contrasts may not be of the form $\mu_{i}-\mu_{j}$ .

The derivation of the Tukey confidence interval of Equation 3.36 for equal sample sizes is straightforward. For the studentized range statistic q, we have

$$
P \left(\frac {\max (\overline {{y}} _ {i .} - \mu_ {i}) - \min (\overline {{y}} _ {i .} - \mu_ {i})}{\sqrt {M S _ {E} / n}} \leq q _ {\alpha} (a, f)\right) = 1 - \alpha
$$

If $\max (\overline{y}_{i.} - \mu_i) - \min (\overline{y}_{i.} - \mu_i)$ is less than or equal to $q_{\alpha}(a,f)\sqrt{MS_E / n}$ , it must be true that $|\overline{y}_{i.} - \mu_i) - \overline{y}_{j.} - \mu_j)|\leq q_{\alpha}(a,f)\sqrt{MS_E / n}$ for every pair of means. Therefore

$$
P \left(- q _ {\alpha} (a, f) \sqrt {\frac {M S _ {E}}{n}} \leq \overline {{{y}}} _ {i.} - \overline {{{y}}} _ {j.} - (\mu_ {i} - \mu_ {j}) \leq q _ {\alpha} (a, f) \sqrt {\frac {M S _ {E}}{n}}\right) = 1 - \alpha
$$

Rearranging this expression to isolate $\mu_{i}-\mu_{j}$ between the inequalities will lead to the set of $100(1-\alpha)$ percent simultaneous confidence intervals given in Equation 3.38.

The Fisher Least Significant Difference (LSD) Method. The Fisher method for comparing all pairs of means controls the error rate $\alpha$ for each individual pairwise comparison but does not control the experimentwise or family error rate. This procedure uses the t statistic for testing $H_{0} : \mu_{i} = \mu_{j}$

$$
t _ {0} = \frac {\overline {{y}} _ {i .} - \overline {{y}} _ {j .}}{\sqrt {M S _ {E} \left(\frac {1}{n _ {i}} + \frac {1}{n _ {j}}\right)}}\tag{3.39}
$$

Assuming a two-sided alternative, the pair of means $\mu_{i}$ and $\mu_{j}$ would be declared significantly different if $|\overline{y}_{i.} - \overline{y}_{j.}| > t_{\alpha /2,N - a}\sqrt{MS_E(1 / n_i + 1 / n_j)}$ . The quantity

$$
\mathrm{LSD} = t _ {\alpha / 2, N - \alpha} \sqrt {M S _ {E} \left(\frac {1}{n _ {i}} + \frac {1}{n _ {j}}\right)}\tag{3.40}
$$

is called the least significant difference. If the design is balanced, $n_{1}=n_{2}=\cdots=n_{a}=n$ , and

$$
\mathrm{LSD} = t _ {\alpha / 2, N - a} \sqrt {\frac {2 M S _ {E}}{n}}\tag{3.41}
$$

To use the Fisher LSD procedure, we simply compare the observed difference between each pair of averages to the corresponding LSD. If $|\overline{y}_{i.}-\overline{y}_{j.}|>$ LSD, we conclude that the population means $\mu_{i}$ and $\mu_{j}$ differ. The t statistic in Equation 3.39 could also be used.

## EXAMPLE 3.8

To illustrate the procedure, if we use the data from the experiment in Example 3.1, the LSD at $\alpha = 0.05$ is

$$
\overline {{{y}}} _ {1.} - \overline {{{y}}} _ {3.} = 5 5 1. 2 - 6 2 5. 4 = - 7 4. 2 ^ {*}
$$

$$
\overline {{{y}}} _ {1.} - \overline {{{y}}} _ {4.} = 5 5 1. 2 - 7 0 7. 0 = - 1 5 5. 8 ^ {*}
$$

$$
\mathrm{LSD} = t _ {. 0 2 5, 1 6} \sqrt {\frac {2 M S _ {E}}{n}} = 2. 1 2 0 \sqrt {\frac {2 (3 3 3 . 7 0)}{5}} = 2 4. 4 9
$$

$$
\overline {{{y}}} _ {2.} - \overline {{{y}}} _ {3.} = 5 8 7. 4 - 6 2 5. 4 = - 3 8. 0 ^ {*}
$$

Thus, any pair of treatment averages that differ in absolute value by more than 24.49 would imply that the corresponding pair of population means are significantly different. The differences in averages are

$$
\overline {{{y}}} _ {2.} - \overline {{{y}}} _ {4.} = 5 8 7. 4 - 7 0 7. 0 = - 1 1 9. 6 ^ {*}
$$

$$
\overline {{{y}}} _ {3.} - \overline {{{y}}} _ {4.} = 6 2 5. 4 - 7 0 7. 0 = - 8 1. 6 ^ {*}
$$

$$
\overline {{{y}}} _ {1.} - \overline {{{y}}} _ {2.} = 5 5 1. 2 - 5 8 7. 4 = - 3 6. 2 ^ {*}
$$

The starred values indicate pairs of means that are significantly different. Clearly, all pairs of means differ significantly.

Note that the overall $\alpha$ risk may be considerably inflated using this method. Specifically, as the number of treatments a gets larger, the experimentwise or family type I error rate (the ratio of the number of experiments in which at least one type I error is made to the total number of experiments) becomes large.

Which Pairwise Comparison Method Do I Use? Certainly, a logical question at this point is as follows: Which one of these procedures should I use? Unfortunately, there is no clear-cut answer to this question, and professional statisticians often disagree over the utility of the various procedures. Carmer and Swanson (1973) have conducted Monte Carlo simulation studies of a number of multiple comparison procedures, including others not discussed here. They report that the least significant difference method is a very effective test for detecting true differences in means if it is applied only after the F-test in the ANOVA is significant at 5 percent. However, this method does not contain the experimentwise error rate. Because the Tukey method does control the overall error rate, many statisticians prefer to use it.

As indicated above, there are several other multiple comparison procedures. For articles describing these methods, see O'Neill and Wetherill (1971), Miller (1977), and Nelson (1989). The books by Miller (1991) and Hsu (1996) are also recommended.

## 3.5.8 Comparing Treatment Means with a Control

In many experiments, one of the treatments is a control, and the analyst is interested in comparing each of the other $a - 1$ treatment means with the control. Thus, only $a - 1$ comparisons are to be made. A procedure for making these comparisons has been developed by Dunnett (1964). Suppose that treatment $a$ is the control and we wish to test the hypotheses

$$
\begin{array}{l} {H _ {0}: \mu_ {i} = \mu_ {a}} \\ {H _ {1}: \mu_ {i} \neq \mu_ {a}} \end{array}
$$

for $i = 1,2,\dots ,a - 1$ . Dunnett's procedure is a modification of the usual $t$ -test. For each hypothesis, we compute the observed differences in the sample means

$$
| \overline {{{y}}} _ {i.} - \overline {{{y}}} _ {a.} | \quad i = 1, 2, \ldots , a - 1
$$

The null hypothesis $H_0: \mu_i = \mu_a$ is rejected using a type I error rate $\alpha$ if

$$
\left| \bar {y} _ {i.} - \bar {y} _ {a.} \right| > d _ {\alpha} (a - 1, f) \sqrt {M S _ {E} \left(\frac {1}{n _ {i}} + \frac {1}{n _ {a}}\right)}\tag{3.42}
$$

where the constant $d_{\alpha}(a-1,f)$ is given in Appendix Table VI. (Both two- and one-sided tests are possible.) Note that $\alpha$ is the joint significance level associated with all a-1 tests.

## EXAMPLE 3.9

To illustrate Dunnett's test, consider the experiment from Example 3.1 with treatment 4 considered as the control. In this example, a = 4, a - 1 = 3, f = 16, and $n_{i} = n = 5$ . At the five percent level, we find from Appendix Table VI that $d_{0.05}(3, 16) = 2.59$ . Thus, the critical difference becomes

$$
d _ {0. 0 5} (3, 1 6) \sqrt {\frac {2 M S _ {E}}{n}} = 2. 5 9 \sqrt {\frac {2 (3 3 3 . 7 0)}{5}} = 2 9. 9 2
$$

(Note that this is a simplification of Equation 3.42 resulting from a balanced design.) Thus, any treatment mean that differs in absolute value from the control by more than 29.92 would be declared significantly different. The observed differences are

$$
1 \text { vs. } 4: \overline {{y}} _ {1.} - \overline {{y}} _ {4.} = 5 5 1. 2 - 7 0 7. 0 = - 1 5 5. 8
$$

$$
2 \text { vs. } 4: \overline {{y}} _ {2.} - \overline {{y}} _ {4.} = 5 8 7. 4 - 7 0 7. 0 = - 1 1 9. 6
$$

$$
3 \mathrm{vs.} 4: \overline {{{y}}} _ {3.} - \overline {{{y}}} _ {4.} = 6 2 5. 4 - 7 0 7. 0 = - 8 1. 6
$$

Note that all differences are significant. Thus, we would conclude that all power settings are different from the control.

When comparing treatments with a control, it is a good idea to use more observations for the control treatment (say $n_{a}$ ) than for the other treatments (say n), assuming equal numbers of observations for the remaining a - 1 treatments. The ratio $n_{a}/n$ should be chosen to be approximately equal to the square root of the total number of treatments. That is, choose $n_{a}/n = \sqrt{a}$ .

## 3.6 Sample Computer Output

Computer programs for supporting experimental design and performing the analysis of variance are widely available. The output from one such program, Design-Expert, is shown in Figure 3.12, using the data from the plasma etching experiment in Example 3.1. The sum of squares corresponding to the “Model” is the usual $SS_{Treatments}$ for a single-factor design. That source is further identified as “A.” When there is more than one factor in the experiment, the model sum of squares will be decomposed into several sources (A, B, etc.). Notice that the analysis of variance summary at the top of the computer output contains the usual sums of squares, degrees of freedom, mean squares, and test statistic $F_{0}$ . The column “Prob > F” is the P-value (actually, the upper bound on the P-value because probabilities less than 0.0001 are defaulted to 0.0001).

In addition to the basic analysis of variance, the program displays some other useful information. The quantity "R-squared" is defined as

$$
R ^ {2} = \frac {S S _ {\text { Model }}}{S S _ {\text { Total }}} = \frac {6 6 , 8 7 0 . 5 5}{7 2 , 2 0 9 . 7 5} = 0. 9 2 6 1
$$

and is loosely interpreted as the proportion of the variability in the data “explained” by the ANOVA model. Thus, in the plasma etching experiment, the factor “power” explains about 92.61 percent of the variability in etch rate. Clearly, we must have $0 \leq R^{2} \leq 1$ , with larger values being more desirable. There are also some other $R^{2}$ -like statistics displayed in the output. The “adjusted” $R^{2}$ is a variation of the ordinary $R^{2}$ statistic that reflects the number of factors in the model. It can be a useful statistic for more complex experiments with several design factors when we wish to evaluate the impact of increasing or decreasing the number of model terms. “Std. Dev.” is the square root of the error mean square, $\sqrt{333.70} = 18.27$ , and “C.V.” is the coefficient of variation, defined as $(\sqrt{MS_{E}}/\overline{y})100$ . The coefficient of variation measures the unexplained or residual variability in the data as a percentage of the mean of the response variable. “PRESS” stands for “prediction error sum of squares,” and it is a measure of how well the model for the experiment is likely to predict the responses in a new experiment. Small values of PRESS are desirable. Alternatively, one can calculate an $R^{2}$ for prediction based on PRESS (we will show how to do this later). This $R_{Pred}^{2}$ in our problem is 0.8845, which is not unreasonable, considering that the model accounts for about 93 percent of the variability in the current experiment. The “adequate precision” statistic is computed by dividing the difference between the maximum predicted response and the minimum predicted response by the average standard deviation of all predicted responses. Large values of this quantity are desirable, and values that exceed four usually indicate that the model will give reasonable performance in prediction.

Treatment means are estimated, and the standard error (or sample standard deviation of each treatment mean, $\sqrt{MS_{E}/n}$ ) is displayed. Differences between pairs of treatment means are investigated by using a hypothesis testing version of the Fisher LSD method described in Section 3.5.7.

The computer program also calculates and displays the residuals, as defined in Equation 3.16. The program will also produce all of the residual plots that we discussed in Section 3.4. There are also several other residual diagnostics displayed in the output. Some of these will be discussed later. Design-Expert also displays the studentized residual (called “Student Residual” in the output) calculated as

$$
r _ {i j} = \frac {e _ {i j}}{\sqrt {M S _ {E} (1 - \mathrm{Leverage} _ {i j})}}
$$

where Leverage $_{ij}$ is a measure of the influence of the ij $^{th}$ observation on the model. We will discuss leverage in more detail and show how it is calculated in Chapter 10. Studentized residuals are considered to be more effective in identifying potential outliers rather than either the ordinary residuals or standardized residuals.

Finally, notice that the computer program also has some interpretative guidance embedded in the output. This “advisory” information is fairly standard in many PC-based statistics packages. Remember in reading such guidance that it is written in very general terms and may not exactly suit the report writing requirements of any specific experimenter. This advisory output may be hidden upon request by the user.

■ FIGURE 3.12 Design-Expert computer output for Example 3.1

## Response: Etch Rate

ANOVA for Selected Factorial Model Analysis of variance table [Partial sum of squares]

<table><tr><td>Source</td><td>Sum of Squares</td><td>DF</td><td>Mean Square</td><td>F Value</td><td>Prob &gt; F</td></tr><tr><td>Model</td><td>66870.55</td><td>3</td><td>22290.18</td><td>66.80</td><td>&lt;0.0001 significant</td></tr><tr><td>A</td><td>66870.55</td><td>3</td><td>22290.18</td><td>66.80</td><td>&lt;0.0001</td></tr><tr><td>Pure Error</td><td>5338.20</td><td>16</td><td>333.70</td><td></td><td></td></tr><tr><td>Cor Total</td><td>72209.75</td><td>19</td><td></td><td></td><td></td></tr></table>

The Model F-value of 66.80 implies that the model is significant. There is only a 0.01% chance that a "Model F-Value" this large could occur due to noise.

Values of "Prob > F" less than 0.0500 indicate that model terms are significant. In this case, A are significant model terms.

Values greater than 0.1000 indicate that the model terms are not significant.
If there are many insignificant model terms (not counting those required to support hierarchy), model reduction may improve your model.

<table><tr><td>Std. Dev.</td><td>18.27</td><td>R-Squared</td><td>0.9261</td></tr><tr><td>Mean</td><td>617.75</td><td>Adj R-Squared</td><td>0.9122</td></tr><tr><td>C.V.</td><td>2.96</td><td>Pred R-Squared</td><td>0.8846</td></tr><tr><td>PRESS</td><td>8342.50</td><td>Adeq Precision</td><td>19.071</td></tr></table>

The "Pred R-Squared" of 0.8845 is in reasonable agreement with the "Adj R-Squared" of 0.9122.

"Adeq Precision" measures the signal-to-noise ratio. A ratio greater than four is desirable. Your ratio of 19.071 indicates an adequate signal. This model can be used to navigate the design space.

Treatment Means (Adjusted, If Necessary)

<table><tr><td></td><td colspan="2">Estimated Mean</td><td colspan="3">Standard Error</td></tr><tr><td>1-160</td><td colspan="2">551.20</td><td colspan="3">8.17</td></tr><tr><td>2-180</td><td colspan="2">587.40</td><td colspan="3">8.17</td></tr><tr><td>3-200</td><td colspan="2">625.40</td><td colspan="3">8.17</td></tr><tr><td>4-220</td><td colspan="2">707.00</td><td colspan="3">8.17</td></tr><tr><td></td><td colspan="2">Mean</td><td rowspan="2">Standard Error</td><td rowspan="2">t for  $H_0$ Coeff=0</td><td rowspan="2">Prob &gt; |t|</td></tr><tr><td>Treatment</td><td>Difference</td><td>DF</td></tr><tr><td>1 vs 2</td><td>-36.20</td><td>1</td><td>11.55</td><td>-3.13</td><td>0.0064</td></tr><tr><td>1 vs 3</td><td>-74.20</td><td>1</td><td>11.55</td><td>-6.42</td><td>&lt;0.0001</td></tr><tr><td>1 vs 4</td><td>-155.80</td><td>1</td><td>11.55</td><td>-13.49</td><td>&lt;0.0001</td></tr><tr><td>2 vs 3</td><td>-38.00</td><td>1</td><td>11.55</td><td>-3.29</td><td>0.0046</td></tr><tr><td>2 vs 4</td><td>-119.60</td><td>1</td><td>11.55</td><td>-10.35</td><td>&lt;0.0001</td></tr><tr><td>3 vs 4</td><td>-81.60</td><td>1</td><td>11.55</td><td>-7.06</td><td>&lt;0.0001</td></tr></table>

Values of "Prob > |t|" less than 0.0500 indicate that the difference in the treatment means is significant.  
Values of "Prob > |t|" greater than 0.1000 indicate that the difference in the two treatment means is not significant.

## Diagnostics Case Statistics

<table><tr><td>Standard Order</td><td>Actual Value</td><td>Predicted Value</td><td>Residual</td><td>Leverage</td><td>Student Residual</td><td>Cook&#x27;s Distance</td><td>Outlier t</td><td>Run Order</td></tr><tr><td>1</td><td>575.00</td><td>551.20</td><td>23.80</td><td>0.200</td><td>1.457</td><td>0.133</td><td>1.514</td><td>13</td></tr><tr><td>2</td><td>542.00</td><td>551.20</td><td>-9.20</td><td>0.200</td><td>-0.563</td><td>0.020</td><td>-0.551</td><td>14</td></tr><tr><td>3</td><td>530.00</td><td>551.20</td><td>-21.20</td><td>0.200</td><td>-1.298</td><td>0.105</td><td>-1.328</td><td>8</td></tr><tr><td>4</td><td>539.00</td><td>551.20</td><td>-12.20</td><td>0.200</td><td>-0.747</td><td>0.035</td><td>-0.736</td><td>5</td></tr><tr><td>5</td><td>570.00</td><td>551.20</td><td>18.80</td><td>0.200</td><td>1.151</td><td>0.083</td><td>1.163</td><td>4</td></tr><tr><td>6</td><td>565.00</td><td>587.40</td><td>-22.40</td><td>0.200</td><td>-1.371</td><td>0.117</td><td>-1.413</td><td>18</td></tr><tr><td>7</td><td>593.00</td><td>587.40</td><td>5.60</td><td>0.200</td><td>0.343</td><td>0.007</td><td>0.333</td><td>9</td></tr><tr><td>8</td><td>590.00</td><td>587.40</td><td>2.60</td><td>0.200</td><td>0.159</td><td>0.002</td><td>0.154</td><td>6</td></tr><tr><td>9</td><td>579.00</td><td>587.40</td><td>-8.40</td><td>0.200</td><td>-0.514</td><td>0.017</td><td>-0.502</td><td>16</td></tr><tr><td>10</td><td>610.00</td><td>587.40</td><td>22.60</td><td>0.200</td><td>1.383</td><td>0.120</td><td>1.427</td><td>17</td></tr><tr><td>11</td><td>600.00</td><td>625.40</td><td>-25.40</td><td>0.200</td><td>-1.555</td><td>0.151</td><td>-1.634</td><td>7</td></tr><tr><td>12</td><td>651.00</td><td>625.40</td><td>25.60</td><td>0.200</td><td>1.567</td><td>0.153</td><td>1.649</td><td>19</td></tr><tr><td>13</td><td>610.00</td><td>625.40</td><td>-15.40</td><td>0.200</td><td>-0.943</td><td>0.056</td><td>-0.939</td><td>10</td></tr><tr><td>14</td><td>637.00</td><td>625.40</td><td>11.60</td><td>0.200</td><td>0.710</td><td>0.032</td><td>0.699</td><td>20</td></tr><tr><td>15</td><td>629.00</td><td>625.40</td><td>3.60</td><td>0.200</td><td>0.220</td><td>0.003</td><td>0.214</td><td>1</td></tr><tr><td>16</td><td>725.00</td><td>707.00</td><td>18.00</td><td>0.200</td><td>1.102</td><td>0.076</td><td>1.110</td><td>2</td></tr><tr><td>17</td><td>700.00</td><td>707.00</td><td>-7.00</td><td>0.200</td><td>-0.428</td><td>0.011</td><td>-0.417</td><td>3</td></tr><tr><td>18</td><td>715.00</td><td>707.00</td><td>8.00</td><td>0.200</td><td>0.490</td><td>0.015</td><td>0.478</td><td>15</td></tr><tr><td>19</td><td>685.00</td><td>707.00</td><td>-22.00</td><td>0.200</td><td>-1.346</td><td>0.113</td><td>-1.385</td><td>11</td></tr><tr><td>20</td><td>710.00</td><td>707.00</td><td>3.00</td><td>0.200</td><td>0.184</td><td>0.002</td><td>0.178</td><td>12</td></tr></table>

Proceed to Diagnostic Plots (the next icon in progression). Be sure to look at the  
(1) Normal probability plot of the studentized residuals to check for normality of residuals.  
(2) Studentized residuals versus predicted values to check for constant error.  
(3) Outlier t versus run order to look for outliers, i.e., influential values.  
(4) Box-Cox plot for power transformations.

If all the model statistics and diagnostic plots are OK, finish up with the Model Graphs icon.

Power = 200 subtracted from

<table><tr><td>Power</td><td>Lower</td><td>Center</td><td>Upper</td></tr><tr><td>180</td><td>3.11</td><td>36.20</td><td>69.29</td></tr><tr><td>200</td><td>41.11</td><td>74.20</td><td>107.29</td></tr><tr><td>220</td><td>122.71</td><td>155.80</td><td>188.89</td></tr></table>

![](images/figure3.12.jpg)

Figure 3.13 presents the output from Minitab for the plasma etching experiment. The output is very similar to the Design-Expert output in Figure 3.12. Note that confidence intervals on each individual treatment mean are provided and that the pairs of means are compared using Tukey's method. However, the Tukey method is presented using the confidence interval format instead of the hypothesis-testing format that we used in Section 3.5.7. None of the Tukey confidence intervals includes zero, so we would conclude that all of the means are different.

Figure 3.14 is the output from JMP for the plasma etch experiment in Example 3.1. The output information is very similar to that from Design-Expert and Minitab. The plots of actual observations versus the predicted values and residuals versus the predicted values are default output. There is an option in JMP to provide the Fisher LSD procedure or Tukey's method to compare all pairs of means.

## One-way ANOVA: Etch Rate versus Power

<table><tr><td>Source</td><td>DF</td><td>SS</td><td>MS</td><td>F</td><td>P</td></tr><tr><td>Power</td><td>3</td><td>66871</td><td>22290</td><td>66.80</td><td>0.000</td></tr><tr><td>Error</td><td>16</td><td>5339</td><td>334</td><td></td><td></td></tr><tr><td>Total</td><td>19</td><td>72210</td><td></td><td></td><td></td></tr><tr><td colspan="6">S = 18.27 R-Sq = 92.61% R-Sq (adj) = 91.22%</td></tr></table>

Pooled Std. Dev. = 18.27  
Turkey 95% Simultaneous Confidence Intervals All Pairwise Comparisons among Levels of Power

Individual confidence level = 98.87%

Power = 160 subtracted from  
Power = 180 subtracted from  
![](images/figure3.13.jpg)  
■ FIGURE 3.13 Minitab computer output for Example 3.1

Response Etch rate
Whole Model

Actual by Predicted Plot  
![](images/40d05bb3f3c18adfd7ff394a6d302d5c1317577e0216f59765875c8e69e4c551.jpg)

Summary of Fit

<table><tr><td>RSquare</td><td>0.92606</td></tr><tr><td>RSquare Adj</td><td>0.912196</td></tr><tr><td>Root Mean Square Error</td><td>18.26746</td></tr><tr><td>Mean of Response</td><td>617.75</td></tr><tr><td>Observations (or Sum Wgts)</td><td>20</td></tr></table>

Analysis of Variance

<table><tr><td>Source</td><td>DF</td><td>Sum of Squares</td><td>Mean Square</td><td>F Ratio</td></tr><tr><td>Model</td><td>3</td><td>66870.550</td><td>22290.2</td><td>66.7971</td></tr><tr><td>Error</td><td>16</td><td>5339.200</td><td>333.7</td><td>Prob &gt; F</td></tr><tr><td>C.Total</td><td>19</td><td>72209.750</td><td></td><td>&lt;.0001</td></tr></table>

<table><tr><td colspan="6">Effect Tests</td></tr><tr><td>Source</td><td>Nparm</td><td>DF</td><td>Sum of Squares</td><td>F Ratio</td><td>Prob &gt; F</td></tr><tr><td>RF power</td><td>3</td><td>3</td><td>66870.550</td><td>66.7971</td><td>&lt;.0001</td></tr></table>

Residual by Predicted Plot  
![](images/6999cf7eec4de2df828002740afa5073d74673684fc628ad82ade8b027f6540b.jpg)  
RF power

Least Squares Means Table

<table><tr><td>Level</td><td>Least Sq Mean</td><td>Std Error</td><td>Mean</td></tr><tr><td>160</td><td>551.20000</td><td>8.1694553</td><td>551.200</td></tr><tr><td>180</td><td>587.40000</td><td>8.1694553</td><td>587.400</td></tr><tr><td>200</td><td>625.40000</td><td>8.1694553</td><td>625.400</td></tr><tr><td>220</td><td>707.00000</td><td>8.1694553</td><td>707.000</td></tr></table>

■ FIGURE 3.14 JMP output from Example 3.1

## 3.7 Determining Sample Size

In any experimental design problem, a critical decision is the choice of sample size—that is, determining the number of replicates to run. Generally, if the experimenter is interested in detecting small effects, more replicates are required than if the experimenter is interested in detecting large effects. In this section, we discuss several approaches to determining sample size. Although our discussion focuses on a single-factor design, most of the methods can be used in more complex experimental situations.

## 3.7.1 Operating Characteristic and Power Curves

Recall that an operating characteristic (OC) curve is a plot of the type II error probability $\beta$ of a statistical test for a particular sample size versus a parameter that reflects the extent to which the null hypothesis is false. Alternatively, a Power Curve plots power or $1 - \beta$ versus this parameter. Power and/or OC curves can be constructed from software and are useful in guiding the experimenter in selecting the number of replicates so that the design will be sensitive to important potential differences in the treatments.

We consider the probability of type II error of the fixed effects model for the case of equal sample sizes per treatment, say

$$
\begin{array}{r l} \beta & = 1 - P \{\text { Reject } H _ {0} | H _ {0} \text { is   false } \} \\ & = 1 - P \{F _ {0} > F _ {\alpha , a - 1, N - a} | H _ {0} \text { is   false } \} \end{array}\tag{3.43}
$$

To evaluate the probability statement in Equation 3.43, we need to know the distribution of the test statistic $F_0$ if the null hypothesis is false. It can be shown that, if $H_0$ is false, the statistic $F_0 = MS_{\text{Treatments}} / MS_E$ is distributed as a noncentral $F$ random variable with $a - 1$ and $N - a$ degrees of freedom and the noncentrality parameter $\delta$ . If $\delta = 0$ , the noncentral $F$ distribution becomes the usual (central) $F$ distribution.

We will illustrate the sample size determination method implemented in JMP. Consider the plasma etching experiment described in Example 3.1. Suppose that the experimenter is interested in rejecting the null hypothesis with a probability of at least 0.9 (power = 0.9) if the true treatment means are

$$
\mu_ {1} = 5 7 5, \mu_ {2} = 6 0 0, \mu_ {3} = 6 5 0, \mathrm{and} \mu_ {1} = 6 7 5
$$

The experimenter feels that the standard deviation of etch rate will be no larger than $\sigma = 25$ Å/min. The input and output from the JMP power and sample size platform for comparing several means are shown in the following display:

![](images/c03uf001.jpg)

The graph on the right is a plot of power versus the total sample size. This plot indicates that at least four replicates are required to obtain a power that exceeds 0.90.

A potential problem with this approach to determining sample size is that it can be difficult to select a set of treatment means on which the sample size decision should be based. An alternate approach is to select a sample size such that if the difference between any two treatment means exceeds a specified value, the null hypothesis should be rejected.

Minitab uses this approach to perform power calculations and find sample sizes for single-factor ANOVAs. Consider the following display:

```txt
Power and Sample Size
One-way ANOVA
Alpha = 0.01 Assumed standard deviation = 25
Number of Levels = 4
Sample Maximum
SS Means Size Power Difference
2812.5 5 0.804838 75
The sample size is for each level.
Power and Sample Size
One-way ANOVA
Alpha = 0.01 Assumed standard deviation = 25
Number of Levels 5 4
Sample Target Maximum
SS Means Size Power Actual Power Difference
2812.5 6 0.9 0.915384 75
The sample size is for each level.
```

In the upper portion of the display, we asked Minitab to calculate the power for n = 5 replicates when the maximum difference in treatment means is 75. The bottom portion of the display is the output when the experimenter requests the sample size to obtain a target power of at least 0.90.

## 3.7.2 Confidence Interval Estimation Method

This approach assumes that the experimenter wishes to express the final results in terms of confidence intervals and is willing to specify in advance how wide he or she wants these confidence intervals to be. For example, suppose that in the plasma etching experiment from Example 3.1, we wanted a 95 percent confidence interval on the difference in mean etch rate for any two power settings to be $\pm30\ \AA/min$ and a prior estimate of $\sigma$ is 25. Then, using Equation 3.13, we find that the accuracy of the confidence interval is

$$
\pm t _ {\alpha / 2, N - a} \sqrt {\frac {2 M S _ {E}}{n}}
$$

Suppose that we try n = 5 replicates. Then, using $\sigma^{2} = (25)^{2} = 625$ as an estimate of $MS_{E}$ , the accuracy of the confidence interval becomes

$$
\pm 2. 1 2 0 \sqrt {\frac {2 (6 2 5)}{5}} = \pm 3 3. 5 2
$$

which does not meet the requirement. Trying n = 6 gives

$$
\pm 2. 0 8 6 \sqrt {\frac {2 (6 2 5)}{6}} = \pm 3 0. 1 1
$$

Trying $n = 7$ gives

$$
\pm 2. 0 6 4 \sqrt {\frac {2 (6 2 5)}{7}} = \pm 2 7. 5 8
$$

Clearly, $n = 7$ is the smallest sample size that will lead to the desired accuracy.

The quoted level of significance in the above illustration applies only to one confidence interval. However, the same general approach can be used if the experimenter wishes to prespecify a set of confidence intervals about which a joint or simultaneous confidence statement is made (see the comments about simultaneous confidence intervals in Section 3.3.3). Furthermore, the confidence intervals could be constructed about more general contrasts in the treatment means than the pairwise comparison illustrated above.

## 3.8 Other Examples of Single-Factor Experiments

## 3.8.1 Chocolate and Cardiovascular Health

An article in Nature describes an experiment to investigate the effect of consuming chocolate on cardiovascular health (“Plasma Antioxidants from Chocolate,” Nature, Vol. 424, 2003, pp. 1013). The experiment consisted of using three different types of chocolates: 100 g of dark chocolate, 100 g of dark chocolate with 200 mL of full-fat milk, and 200 g of milk chocolate. A total of 12 subjects were used, seven women and five men, with an average age range of $32.2 \pm 1$ years, an average weight of $65.8 \pm 3.1$ kg, and body-mass index of $21.9 \pm 0.4$ kg m $^{-2}$ . On different days a subject consumed one of the chocolate-factor levels and 1 hour later the total antioxidant capacity of their blood plasma was measured in an assay. Data similar to that summarized in the article are shown in Table 3.12.

Figure 3.15 presents box plots for the data from this experiment. The result is an indication that the blood antioxidant capacity one hour after eating the dark chocolate is higher than for the other two treatments. The variability in the sample data from all three treatments seems very similar. Table 3.13 is the Minitab ANOVA output. The test statistic is highly significant (Minitab reports a P-value of 0.000, which is clearly wrong because P-values cannot be zero; this means that the P-value is less than 0.001), indicating that some of the treatment means are different. The output also contains the Fisher LSD analysis for this experiment. This indicates that the mean antioxidant capacity after consuming dark chocolate is higher than after consuming dark chocolate plus milk or milk chocolate alone. Furthermore, the mean antioxidant capacity after consuming dark chocolate plus milk or milk chocolate alone is equal. Figure 3.16 is the normal probability plot of the residual and Figure 3.17 is the plot of residuals versus predicted values. These plots do not suggest any problems with model assumptions. We conclude that consuming dark chocolate results in higher mean blood antioxidant capacity after one hour than consuming either dark chocolate plus milk or milk chocolate alone.

## TABLE 3.12

Blood Plasma Levels One Hour Following Chocolate Consumption

<table><tr><td rowspan="2">Factor</td><td colspan="12">Subjects (Observations)</td></tr><tr><td>1</td><td>2</td><td>3</td><td>4</td><td>5</td><td>6</td><td>7</td><td>8</td><td>9</td><td>10</td><td>11</td><td>12</td></tr><tr><td>DC</td><td>118.8</td><td>122.6</td><td>115.6</td><td>113.6</td><td>119.5</td><td>115.9</td><td>115.8</td><td>115.1</td><td>116.9</td><td>115.4</td><td>115.6</td><td>107.9</td></tr><tr><td>DC + MK</td><td>105.4</td><td>101.1</td><td>102.7</td><td>97.1</td><td>101.9</td><td>98.9</td><td>100.0</td><td>99.8</td><td>102.6</td><td>100.9</td><td>104.5</td><td>93.5</td></tr><tr><td>MC</td><td>102.1</td><td>105.8</td><td>99.6</td><td>102.7</td><td>98.8</td><td>100.9</td><td>102.8</td><td>98.7</td><td>94.7</td><td>97.8</td><td>99.7</td><td>98.6</td></tr></table>

■ FIGURE 3.15 Box plots of the blood antioxidant capacity data from the chocolate consumption experiment

![](images/figure3.15.jpg)

```txt
One-way ANOVA: DC, DC+MK, MC
Source DF SS MS F P
Factor 2 1952.6 976.3 93.58 0.000
Error 33 344.3 10.4
Total 35 2296.9
S = 3.230 R-Sq = 85.01% R-Sq(adj) = 84.10%
Individual 95% CIs For Mean Based on
Pooled StDev
Level N Mean StDev ----+----+----+----+----
DC 12 116.06 3.53 (----*----)
DC+MK 12 100.70 3.24 (--*----)
MC 12 100.18 2.89 (--*----)
100.0 105.0 110.0 115.0
Pooled StDev = 3.23
Fisher 95% Individual Confidence Intervals
All Pairwise Comparisons
Simultaneous confidence level = 88.02
DC subtracted from:
Lower Center Upper -+----+----+----+----
DC+MK -18.041 -15.358 -12.675 (--*----)
MC -18.558 -15.875 -13.192 (--*----)
-+----+----+----+----
-18.0 -12.0 -6.0 0.0
DC+MK subtracted from:
Lower Center Upper -+----+----+----+----
MC -3.200 -0.517 2.166 (----*----)
-+----+----+----+----+----
-18.0 -12.0 -6.0 0.0
```

![](images/figure3.16.jpg)  
■ FIGURE 3.16 Normal probability plot of the residuals from the chocolate consumption experiment

![](images/figure3.17.jpg)  
■ FIGURE 3.17 Plot of residuals versus the predicted values from the chocolate consumption experiment

## 3.8.2 A Real Economy Application of a Designed Experiment

Designed experiments have had tremendous impact on manufacturing industries, including the design of new products and the improvement of existing ones, development of new manufacturing processes, and process improvement. In the last 15 years, designed experiments have begun to be widely used outside of this traditional environment. These applications are in financial services, telecommunications, health care, e-commerce, legal services, marketing, logistics and transportation, and many of the nonmanufacturing components of manufacturing businesses. These types of businesses are sometimes referred to as the real economy. It has been estimated that manufacturing accounts for only about 20 percent of the total US economy, so applications of experimental design in the real economy are of growing importance. In this section, we present an example of a designed experiment in marketing.

A soft drink distributor knows that end-aisle displays are an effective way to increase sales of the product. However, there are several ways to design these displays: by varying the text displayed, the colors used, and the visual images. The marketing group has designed three new end-aisle displays and wants to test their effectiveness. They have identified 15 stores of similar size and type to participate in the study. Each store will test one of the displays for a period of one month. The displays are assigned at random to the stores, and each display is tested in five stores. The response variable is the percentage increase in sales activity over the typical sales for that store when the end-aisle display is not in use. The data from this experiment are shown in Table 3.14.

Table 3.15 shows the analysis of the end-aisle display experiment. This analysis was conducted using JMP. The P-value for the model F-statistic in the ANOVA indicates that there is a difference in the mean percentage increase in sales between the three display types. In this application, we had JMP use the Fisher LSD procedure to compare the

TABLE 3.14  
The End-Aisle Display Experimental Design

<table><tr><td>Display Design</td><td></td><td colspan="4">Sample Observations, Percent Increase in Sales</td></tr><tr><td>1</td><td>5.43</td><td>5.71</td><td>6.22</td><td>6.01</td><td>5.29</td></tr><tr><td>2</td><td>6.24</td><td>6.71</td><td>5.98</td><td>5.66</td><td>6.60</td></tr><tr><td>3</td><td>8.79</td><td>9.20</td><td>7.90</td><td>8.15</td><td>7.55</td></tr></table>

Response Sales Increase

Whole Model

Actual by Predicted Plot  
![](images/c03uf002a.jpg)

Summary of Fit

<table><tr><td>RSquare</td><td>0.856364</td></tr><tr><td>RSquare Adj</td><td>0.832425</td></tr><tr><td>Root Mean Square Error</td><td>0.512383</td></tr><tr><td>Mean of Response</td><td>6.762667</td></tr><tr><td>Observations (or Sum Wgts)</td><td>15</td></tr></table>

Analysis of Variance

<table><tr><td>Source</td><td>DF</td><td>Sum of Squares</td><td>Mean Square</td><td>F Ratio</td></tr><tr><td>Model</td><td>2</td><td>18.783053</td><td>9.39153</td><td>35.7722</td></tr><tr><td>Error</td><td>12</td><td>3.150440</td><td>0.26254</td><td>Prob &gt; F</td></tr><tr><td>C.Total</td><td>14</td><td>21.933493</td><td></td><td>&lt; .0001</td></tr></table>

Effect Tests

<table><tr><td>Source</td><td>Nparm</td><td>DF</td><td>Sum of Squares</td><td>F Ratio</td><td>Prob &gt; F</td></tr><tr><td>Display</td><td>2</td><td>2</td><td>18.783053</td><td>35.7722</td><td>&lt; .001</td></tr></table>

Residual by Predicted Plot  
![](images/c03uf002b.jpg)

■ TABLE 3.15 (Continued)

<table><tr><td colspan="4">Least Squares Means Table</td></tr><tr><td>Level</td><td>Least Sq Mean</td><td>Std Error</td><td>Mean</td></tr><tr><td>1</td><td>5.7320000</td><td>0.22914479</td><td>5.73200</td></tr><tr><td>2</td><td>6.2380000</td><td>0.22914479</td><td>6.23800</td></tr><tr><td>3</td><td>8.3180000</td><td>0.22914479</td><td>8.31800</td></tr></table>

LSMeans Differences Student's t  
$a = 0.050t = 2.17881$  
LSMean[i] By LSMean [i]

<table><tr><td>Mean[i]-Mean [i]Std Err DifLower CL DifUpper CL Dif</td><td>1</td><td>2</td><td>3</td></tr><tr><td rowspan="4">1</td><td>0</td><td>-0.506</td><td>-2.586</td></tr><tr><td>0</td><td>0.32406</td><td>-2.586</td></tr><tr><td>0</td><td>-1.2121</td><td>-3.2921</td></tr><tr><td>0</td><td>0.20007</td><td>-1.8799</td></tr><tr><td rowspan="4">2</td><td>0.506</td><td>0</td><td>-2.08</td></tr><tr><td>0.32406</td><td>0</td><td>0.32406</td></tr><tr><td>-0.2001</td><td>0</td><td>-2.7861</td></tr><tr><td>1.21207</td><td>0</td><td>-1.3739</td></tr><tr><td rowspan="4">3</td><td>2.586</td><td>2.08</td><td>0</td></tr><tr><td>0.32406</td><td>0.32406</td><td>0</td></tr><tr><td>1.87993</td><td>1.37393</td><td>0</td></tr><tr><td>3.29207</td><td>2.78607</td><td>0</td></tr></table>

<table><tr><td>Level</td><td></td><td>Least Sq Mean</td></tr><tr><td>3</td><td>A</td><td>8.3180000</td></tr><tr><td>2</td><td>B</td><td>6.2380000</td></tr><tr><td>1</td><td>B</td><td>5.7320000</td></tr></table>

Levels not connected by same letter are significantly different.

pairs of treatment means (JMP labels these as the least squares means). The results of this comparison are presented as confidence intervals on the difference in pairs of means. For pairs of means where the confidence interval includes zero, we would not declare that the pairs of means are different. The JMP output indicates that display designs 1 and 2 are similar in that they result in the same mean increase in sales, but that display design 3 is different from both designs 1 and 2 and that the mean increase in sales for display 3 exceeds that of both designs 1 and 2. Notice that JMP automatically includes some useful graphics in the output, a plot of the actual observations versus the predicted values from the model, and a plot of the residuals versus the predicted values. There is some mild indication that display design 3 may exhibit more variability in sales increase than the other two designs.

## 3.8.3 Discovering Dispersion Effects

We have focused on using the analysis of variance and related methods to determine which factor levels result in differences among treatment or factor level means. It is customary to refer to these effects as location effects. If there was inequality of variance at the different factor levels, we used transformations to stabilize the variance to improve our inference on the location effects. In some problems, however, we are interested in discovering whether the different factor levels affect variability; that is, we are interested in discovering potential dispersion effects. This will occur whenever the standard deviation, variance, or some other measure of variability is used as a response variable.

TABLE 3.16  
Data for the Smelting Experiment

<table><tr><td rowspan="2">Ratio Control Algorithm</td><td colspan="6">Observations</td></tr><tr><td>1</td><td>2</td><td>3</td><td>4</td><td>5</td><td>6</td></tr><tr><td>1</td><td>4.93(0.05)</td><td>4.86(0.04)</td><td>4.75(0.05)</td><td>4.95(0.06)</td><td>4.79(0.03)</td><td>4.88(0.05)</td></tr><tr><td>2</td><td>4.85(0.04)</td><td>4.91(0.02)</td><td>4.79(0.03)</td><td>4.85(0.05)</td><td>4.75(0.03)</td><td>4.85(0.02)</td></tr><tr><td>3</td><td>4.83(0.09)</td><td>4.88(0.13)</td><td>4.90(0.11)</td><td>4.75(0.15)</td><td>4.82(0.08)</td><td>4.90(0.12)</td></tr><tr><td>4</td><td>4.89(0.03)</td><td>4.77(0.04)</td><td>4.94(0.05)</td><td>4.86(0.05)</td><td>4.79(0.03)</td><td>4.76(0.02)</td></tr></table>

To illustrate these ideas, consider the data in Table 3.16, which resulted from a designed experiment in an aluminum smelter. Aluminum is produced by combining alumina with other ingredients in a reaction cell and applying heat by passing electric current through the cell. Alumina is added continuously to the cell to maintain the proper ratio of alumina to other ingredients. Four different ratio control algorithms were investigated in this experiment. The response variables studied were related to cell voltage. Specifically, a sensor scans cell voltage several times each second, producing thousands of voltage measurements during each run of the experiment. The process engineers decided to use the average voltage and the standard deviation of cell voltage (shown in parentheses) over the run as the response variables. The average voltage is important because it affects cell temperature, and the standard deviation of voltage (called “pot noise” by the process engineers) is important because it affects the overall cell efficiency.

An analysis of variance was performed to determine whether the different ratio control algorithms affect average cell voltage. This revealed that the ratio control algorithm had no location effect; that is, changing the ratio control algorithms does not change the average cell voltage. (Refer to Problem 3.33.)

To investigate dispersion effects, it is usually best to use

$$
\log (s) \quad \text { or } \quad \log (s ^ {2})
$$

as a response variable since the log transformation is effective in stabilizing variability in the distribution of the sample standard deviation. Because all sample standard deviations of pot voltage are less than unity, we will use

$$
y = - \ln (s)
$$

as the response variable. Table 3.17 presents the analysis of variance for this response, the natural logarithm of “pot noise.” Notice that the choice of a ratio control algorithm affects pot noise; that is, the ratio control algorithm has a dispersion effect. Standard tests of model adequacy, including normal probability plots of the residuals, indicate that there are no problems with experimental validity. (Refer to Problem 3.34.)

## TABLE 3.17

Analysis of Variance for the Natural Logarithm of Pot Noise

<table><tr><td>Source of Variation</td><td>Sum of Squares</td><td>Degrees of Freedom</td><td>Mean Square</td><td> $F_0$ </td><td>P-Value</td></tr><tr><td>Ratio control algorithm</td><td>6.166</td><td>3</td><td>2.055</td><td>21.96</td><td>&lt; 0.001</td></tr><tr><td>Error</td><td>1.872</td><td>20</td><td>0.094</td><td></td><td></td></tr><tr><td>Total</td><td>8.038</td><td>23</td><td></td><td></td><td></td></tr></table>

![](images/figure3.18.jpg)

$$
\sqrt {M S _ {E} / n} =
$$

Figure 3.18 plots the average log pot noise for each ratio control algorithm and also presents a scaled t distribution for use as a reference distribution in discriminating between ratio control algorithms. This plot clearly reveals that ratio control algorithm 3 produces greater pot noise or greater cell voltage standard deviation than the other algorithms. There does not seem to be much difference between algorithms 1, 2, and 4.

## 3.9 The Random Effects Model

## 3.9.1 A Single Random Factor

An experimenter is frequently interested in a factor that has a large number of possible levels. If the experimenter randomly selects $a$ of these levels from the population of factor levels, then we say that the factor is random. Because the levels of the factor actually used in the experiment were chosen randomly, inferences are made about the entire population of factor levels. We assume that the population of factor levels is either of infinite size or is large enough to be considered infinite. Situations in which the population of factor levels is small enough to employ a finite population approach are not encountered frequently. Refer to Bennett and Franklin (1954) and Searle and Fawcett (1970) for a discussion of the finite population case.

The linear statistical model is

$$
y _ {i j} = \mu + \tau_ {i} + \epsilon_ {i j} \left\{ \begin{array}{l} i = 1, 2, \ldots , a \\ j = 1, 2, \ldots , n \end{array} \right.\tag{3.44}
$$

where both the treatment effects $\tau_{i}$ and $\epsilon_{ij}$ are random variables. We will assume that the treatment effects $\tau_{i}$ are NID $(0,\sigma_{\tau}^{2})$ , random variables $^{3}$ and that the errors are NID $(0,\sigma^{2})$ , random variables, and that the $\tau_{i}$ and $\epsilon_{ij}$ are independent. Because $\tau_{i}$ is independent of $\epsilon_{ij}$ , the variance of any observation is

$$
V (y _ {i j}) = \sigma_ {\tau} ^ {2} + \sigma^ {2}
$$

The variances $\sigma_{\tau}^{2}$ and $\sigma^{2}$ are called variance components, and the model (Equation 3.44) is called the components of variance or random effects model. The observations in the random effects model are normally distributed because they are linear combinations of the two normally and independently distributed random variables $\tau_{i}$ and $\epsilon_{ij}$ . However, unlike the fixed effects case in which all of the observations $y_{ij}$ are independent, in the random model the observations $y_{ij}$ are only independent if they come from different factor levels. Specifically, we can show that the covariance of any two observations is

$$
\begin{array}{l l} {C o v (y _ {i j}, y _ {i j ^ {\prime}}) = \sigma_ {\tau} ^ {2}} & {j \neq j ^ {\prime}} \\ {C o v (y _ {i j}, y _ {i ^ {\prime} j ^ {\prime}}) = 0} & {i \neq i ^ {\prime}} \end{array}
$$

Note that the observations within a specific factor level all have the same covariance, because before the experiment is conducted, we expect the observations at that factor level to be similar because they all have the same random component. Once the experiment has been conducted, we can assume that all observations can be assumed to be independent, because the parameter $\tau_{i}$ has been determined and the observations in that treatment differ only because of random error.

We can express the covariance structure of the observations in the single-factor random effects model through the covariance matrix of the observations. To illustrate, suppose that we have a = 3 treatments and n = 2 replicates. There are N = 6 observations, which we can write as a vector

$$
\mathbf {y} = \left[ \begin{array}{c} y _ {1 1} \\ y _ {1 2} \\ y _ {2 1} \\ y _ {2 2} \\ y _ {3 1} \\ y _ {3 2} \end{array} \right]
$$

and the $6 \times 6$ covariance matrix of these observations is

$$
C o v (\mathbf {y}) = \left[ \begin{array}{c c c c c c} \sigma_ {\tau} ^ {2} + \sigma^ {2} & \sigma_ {\tau} ^ {2} & 0 & 0 & 0 & 0 \\ \sigma_ {\tau} ^ {2} & \sigma_ {\tau} ^ {2} + \sigma^ {2} & 0 & 0 & 0 & 0 \\ 0 & 0 & \sigma_ {\tau} ^ {2} + \sigma^ {2} & \sigma_ {\tau} ^ {2} & 0 & 0 \\ 0 & 0 & \sigma_ {\tau} ^ {2} & \sigma_ {\tau} ^ {2} + \sigma^ {2} & 0 & 0 \\ 0 & 0 & 0 & 0 & \sigma_ {\tau} ^ {2} + \sigma^ {2} & \sigma^ {2} \\ 0 & 0 & 0 & 0 & \sigma_ {\tau} ^ {2} & \sigma_ {\tau} ^ {2} + \sigma^ {2} \end{array} \right]
$$

The main diagonals of this matrix are the variances of each individual observation and every off-diagonal element is the covariance of a pair of observations.

## 3.9.2 Analysis of Variance for the Random Model

The basic ANOVA sum of squares identity

$$
S S _ {T} = S S _ {\text { Treatments }} + S S _ {E}\tag{3.45}
$$

is still valid. That is, we partition the total variability in the observations into a component that measures the variation between treatments ( $SS_{Treatments}$ ) and a component that measures the variation within treatments ( $SS_{E}$ ). Testing hypotheses about individual treatment effects is not very meaningful because they were selected randomly, we are more interested in the population of treatments, so we test hypotheses about the variance component $\sigma_{\tau}^{2}$ .

$$
\begin{array}{l} H _ {0}: \sigma_ {\tau} ^ {2} = 0 \\ H _ {1}: \sigma_ {\tau} ^ {2} > 0 \end{array}\tag{3.46}
$$

If $\sigma_{\tau}^{2} = 0$ , all treatments are identical; but if $\sigma_{\tau}^{2} = 0$ , variability exists between treatments. As before, $SS_{E} / \sigma^{2}$ is distributed as chi-square with $N - a$ degrees of freedom and, under the null hypothesis, $SS_{\text{Treatments}} / \sigma^{2}$ is distributed as chi-square with $a - 1$ degrees of freedom. Both random variables are independent. Thus, under the null hypothesis $\sigma_{\tau}^{2} = 0$ , the ratio

$$
F _ {0} = \frac {\frac {S S _ {\text { Treatments }}}{a - 1}}{\frac {S S _ {E}}{N - a}} = \frac {M S _ {\text { Treatments }}}{M S _ {E}}\tag{3.47}
$$

is distributed as $F$ with $a - 1$ and $N - a$ degrees of freedom. However, we need to examine the expected mean squares to fully describe the test procedure.

Consider

$$
\begin{array}{l} E (M S _ {\text {Treatments}}) = \frac {1}{a - 1} E (S S _ {\text {Treatments}}) = \frac {1}{a - 1} E \left[ \sum_ {i = 1} ^ {a} \frac {y _ {i .} ^ {2}}{n} - \frac {y _ {. .} ^ {2}}{N} \right] \\ = \frac {1}{a - 1} E \left[ \frac {1}{n} \sum_ {i = 1} ^ {a} \left(\sum_ {j = 1} ^ {n} \mu + \tau_ {i} + \epsilon_ {i j}\right) ^ {2} - \frac {1}{N} \left(\sum_ {i = 1} ^ {a} \sum_ {j = 1} ^ {n} \mu + \tau_ {i} + \epsilon_ {i j}\right) ^ {2} \right] \end{array}
$$

When squaring and taking expectation of the quantities in brackets, we see that terms involving $\tau_{i}^{2}$ are replaced by $\sigma_{\tau}^{2}$ as $E(\tau_{i}) = 0$ . Also, terms involving $\epsilon_{i.}^{2}, \epsilon_{..}^{2}$ , and $\sum_{i=1}^{a} \sum_{j=1}^{n} \tau_{i}^{2}$ are replaced by $n\sigma^{2}, an\sigma^{2}$ , and $an^{2}$ , respectively. Furthermore, all cross-product terms involving $\tau_{i}$ and $\epsilon_{ij}$ have zero expectation. This leads to

$$
E (M S _ {\text { Treatments }}) = \frac {1}{a - 1} [ N \mu^ {2} + N \sigma_ {\tau} ^ {2} + a \sigma^ {2} - N \mu^ {2} - n \sigma_ {\tau} ^ {2} - \sigma^ {2} ]
$$

or

$$
E (M S _ {\text { Treatments }}) = \sigma^ {2} + n \sigma_ {\tau} ^ {2}\tag{3.48}
$$

Similarly, we may show that

$$
E (M S _ {E}) = \sigma^ {2}\tag{3.49}
$$

From the expected mean squares, we see that under $H_{0}$ both the numerator and denominator of the test statistic (Equation 3.47) are unbiased estimators of $\sigma^{2}$ , whereas under $H_{1}$ the expected value of the numerator is greater than the expected value of the denominator. Therefore, we should reject $H_{0}$ for values of $F_{0}$ that are too large. This implies an upper-tail, one-tail critical region, so we reject $H_{0}$ if $F_{0} > F_{\alpha,a-1,N-a}$ .

The computational procedure and ANOVA for the random effects model are identical to those for the fixed effects case. The conclusions, however, are quite different because they apply to the entire population of treatments.

## 3.9.3 Estimating the Model Parameters

We are usually interested in estimating the variance components $(\sigma^{2}$ and $\sigma_{\tau}^{2})$ in the model. One very simple procedure that we can use to estimate $\sigma^{2}$ and $\sigma_{\tau}^{2}$ is called the analysis of variance method because it makes use of the lines in the analysis of variance table. The procedure consists of equating the expected mean squares to their observed values in the ANOVA table and solving for the variance components. In equating observed and expected mean squares in the single-factor random effects model, we obtain

$$
M S _ {\text { Treatments }} = \sigma^ {2} + n \sigma_ {\tau} ^ {2}
$$

and

$$
M S _ {E} = \sigma^ {2}
$$

Therefore, the estimators of the variance components are

$$
\hat {\sigma} ^ {2} = M S _ {E}\tag{3.50}
$$

and

$$
\hat {\sigma} _ {\tau} ^ {2} = \frac {M S _ {\text { Treatments }} - M S _ {E}}{n}\tag{3.51}
$$

For unequal sample sizes, replace $n$ in Equation 3.51 by

$$
n _ {0} = \frac {1}{a - 1} \left[ \sum_ {i = 1} ^ {a} n _ {i} - \frac {\sum_ {i = 1} ^ {a} n _ {i} ^ {2}}{\sum_ {i = 1} ^ {a} n _ {i}} \right]\tag{3.52}
$$

The analysis of variance method of variance component estimation is a method of moments procedure. It does not require the normality assumption. It does yield estimators of $\sigma^{2}$ and $\sigma_{\tau}^{2}$ that are best quadratic unbiased (i.e., of all unbiased quadratic functions of the observations, these estimators have minimum variance). There is a different method based on maximum likelihood that can be used to estimate the variance components that will be introduced later.

Occasionally, the analysis of variance method produces a negative estimate of a variance component. Clearly, variance components are by definition nonnegative, so a negative estimate of a variance component is viewed with some concern. One course of action is to accept the estimate and use it as evidence that the true value of the variance component is zero, assuming that sampling variation led to the negative estimate. This has intuitive appeal, but it suffers from some theoretical difficulties. For instance, using zero in place of the negative estimate can disturb the statistical properties of other estimates. Another alternative is to reestimate the negative variance component using a method that always yields nonnegative estimates. Still another alternative is to consider the negative estimate as evidence that the assumed linear model is incorrect and reexamine the problem. Comprehensive treatment of variance component estimation is given by Searle (1971a, 1971b), Searle, Casella, and McCullogh (1992), and Burdick and Graybill (1992).

## EXAMPLE 3.10

A textile company weaves a fabric on a large number of looms. It would like the looms to be homogeneous so that it obtains a fabric of uniform strength. The process engineer suspects that, in addition to the usual variation in strength within samples of fabric from the same loom, there may also be significant variations in strength between looms. To investigate this, she selects four looms at random and makes four strength determinations on the fabric manufactured on each loom. This experiment is run in random order, and the data obtained are shown in Table 3.18. The ANOVA is conducted and is shown in Table 3.19. From the ANOVA, we conclude that the looms in the plant differ significantly.

The variance components are estimated by $\hat{\sigma}^2 = 1.90$ and

$$
\hat {\sigma} _ {\tau} ^ {2} = \frac {2 9 . 7 3 - 1 . 9 0}{4} = 6. 9 6
$$

Therefore, the variance of any observation on strength is estimated by

$$
\hat {\sigma} _ {y} = \hat {\sigma} ^ {2} + \hat {\sigma} _ {\tau} ^ {2} = 1. 9 0 + 6. 9 6 = 8. 8 6.
$$

Most of this variability is attributable to differences between looms.

TABLE 3.18  
Strength Data for Example 3.10

<table><tr><td rowspan="2">Looms</td><td colspan="4">Observations</td><td rowspan="2"> $y_{i.}$ </td></tr><tr><td>1</td><td>2</td><td>3</td><td>4</td></tr><tr><td>1</td><td>98</td><td>97</td><td>99</td><td>96</td><td>390</td></tr><tr><td>2</td><td>91</td><td>90</td><td>93</td><td>92</td><td>366</td></tr><tr><td>3</td><td>96</td><td>95</td><td>97</td><td>95</td><td>383</td></tr><tr><td>4</td><td>95</td><td>96</td><td>99</td><td>98</td><td>388</td></tr></table>

$$
1 5 2 7 = y _ {\cdot \cdot}
$$

TABLE 3.19  
Analysis of Variance for the Strength Data

<table><tr><td>Source of Variation</td><td>Sum of Squares</td><td>Degrees of Freedom</td><td>Mean Square</td><td> $F_0$ </td><td>P-Value</td></tr><tr><td>Looms</td><td>89.19</td><td>3</td><td>29.73</td><td>15.68</td><td>&lt;0.001</td></tr><tr><td>Error</td><td>22.75</td><td>12</td><td>1.90</td><td></td><td></td></tr><tr><td>Total</td><td>111.94</td><td>15</td><td></td><td></td><td></td></tr></table>

![](images/figure3.19.jpg)  

■ FIGURE 3.19 Process output in the fiber strength problem  
(b) Variability of process output if $\sigma_{\tau}^{2} = 0$

This example illustrates an important use of variance components—isolating different sources of variability that affect a product or system. The problem of product variability frequently arises in quality assurance, and it is often difficult to isolate the sources of variability. For example, this study may have been motivated by an observation that there is too much variability in the strength of the fabric, as illustrated in Figure 3.19a. This graph displays the process output (fiber strength) modeled as a normal distribution with variance $\hat{\sigma}_{y}^{2}=8.86$ . (This is the estimate of the variance of any observation on strength from Example 3.10.) Upper and lower specifications on strength are also shown in Figure 3.19a, and it is relatively easy to see that a fairly large proportion of the process output is outside the specifications (the shaded tail areas in Figure 3.19a). The process engineer has asked why so much fabric is defective and must be scrapped, reworked, or downgraded to a lower quality product. The answer is that most of the product strength variability is the result of differences between looms. Different loom performance could be the result of faulty setup, poor maintenance, ineffective supervision, poorly trained operators, defective input fiber, and so forth.

The process engineer must now try to isolate the specific causes of the differences in loom performance. If she could identify and eliminate these sources of between-loom variability, the variance of the process output could be reduced considerably, perhaps to as low as $\hat{\sigma}_{y}^{2}=1.90$ , the estimate of the within-loom (error) variance component in Example 3.10. Figure 3.19b shows a normal distribution of fiber strength with $\hat{\sigma}_{y}^{2}=1.90$ . Note that the proportion of defective product in the output has been dramatically reduced. Although it is unlikely that all of the between-loom variability can be eliminated, it is clear that a significant reduction in this variance component would greatly increase the quality of the fiber produced.

We may easily find a confidence interval for the variance component $\sigma^{2}$ . If the observations are normally and independently distributed, then $(N-a)MS_{E}/\sigma^{2}$ is distributed as $\chi_{N-a}^{2}$ . Thus,

$$
P \left[ \chi_ {1 - (\alpha / 2), N - a} ^ {2} \leq \frac {(N - a) M S _ {E}}{\sigma^ {2}} \leq \chi_ {\alpha / 2, N - a} ^ {2} \right] = 1 - \alpha
$$

and a $100(1 - \alpha)$ percent confidence interval for $\sigma^2$ is

$$
\frac {(N - a) M S _ {E}}{\chi_ {\alpha / 2 , N - a} ^ {2}} \leq \sigma^ {2} \leq \frac {(N - a) M S _ {E}}{\chi_ {1 - (\alpha / 2) , N - a} ^ {2}}\tag{3.53}
$$

Since $MS_{E} = 190$ , $N = 16$ , $a = 4$ , $\chi^{2}_{0.025,12} = 23,3367$ and $\chi^{2}_{0.975,12} = 4.4038$ , the $95\%$ CI on $\sigma^2$ is $0.9770 \leq \sigma^2 \leq 5.1775$ .

Now consider the variance component $\sigma_{\tau}^{2}$ . The point estimator of $\sigma_{\tau}^{2}$ is

$$
\hat {\sigma} _ {\tau} ^ {2} = \frac {M S _ {\text { Treatments }} - M S _ {E}}{n}
$$

The random variable $(a-1)MS_{\text{Treatments}}/(\sigma^{2}+n\sigma_{\tau}^{2})$ is distributed as $\chi_{a-1}^{2}$ , and $(N-a)MS_{E}/\sigma^{2}$ is distributed as $\chi_{N-a}^{2}$ . Thus, the probability distribution of $\hat{\sigma}_{\tau}^{2}$ is a linear combination of two chi-square random variables, say

$$
u _ {1} \chi_ {a - 1} ^ {2} - u _ {2} \chi_ {N - a} ^ {2}
$$

where

$$
u _ {1} = \frac {\sigma^ {2} + n \sigma_ {\tau} ^ {2}}{n (a - 1)} \quad \mathrm{and} \quad u _ {2} = \frac {\sigma^ {2}}{n (N - a)}
$$

Unfortunately, a closed-form expression for the distribution of this linear combination of chi-square random variables cannot be obtained. Thus, an exact confidence interval for $\sigma_{\tau}^{2}$ cannot be constructed. Approximate procedures are given in Graybill (1961) and Searle (1971a). Also see Section 13.6 of Chapter 13.

It is easy to find an exact expression for a confidence interval on the ratio $\sigma_{\tau}^{2}/(\sigma_{\tau}^{2}+\sigma^{2})$ . This ratio is called the intraclass correlation coefficient, and it reflects the proportion of the variance of an observation [recall that $V(y_{ij})=\sigma_{\tau}^{2}+\sigma^{2}]$ that is the result of differences between treatments. To develop this confidence interval for the case of a balanced design, note that $MS_{Treatments}$ and $MS_{E}$ are independent random variables and, furthermore, it can be shown that

$$
\frac {M S _ {\mathrm{Treatments}} / (n \sigma_ {\tau} ^ {2} + \sigma^ {2})}{M S _ {E} / \sigma^ {2}} \sim F _ {a - 1, N - a}
$$

Thus,

$$
\mathrm{P} \left(F _ {1 - \alpha / 2, a - 1, N - a} \leq \frac {M S _ {\text {Treatments}}}{M S _ {E}} \frac {\sigma^ {2}}{n \sigma_ {\tau} ^ {2} + \sigma^ {2}} \leq F _ {\alpha / 2, a - 1, N - a}\right) = 1 - \alpha\tag{3.54}
$$

By rearranging Equation 3.54, we may obtain the following:

$$
P \left(L \leq \frac {\sigma_ {\tau} ^ {2}}{\sigma^ {2}} \leq U\right) = 1 - \alpha\tag{3.55}
$$

where

$$
L = \frac {1}{n} \left(\frac {M S _ {\text { Treatments }}}{M S _ {E}} \frac {1}{F _ {\alpha / 2 , a - 1 , N - a}} - 1\right)\tag{3.56a}
$$

and

$$
U = \frac {1}{n} \left(\frac {M S _ {\mathrm{Treatments}}}{M S _ {E}} \frac {1}{F _ {1 - \alpha / 2 , a - 1 , N - a}} - 1\right)\tag{3.56b}
$$

Note that L and U are $100(1 - \alpha)$ percent lower and upper confidence limits, respectively, for the ratio $\sigma_{\tau}^{2}/\sigma^{2}$ . Therefore, a $100(1 - \alpha)$ percent confidence interval for $\sigma_{\tau}^{2}/(\sigma_{\tau}^{2} + \sigma^{2})$ is

$$
\frac {L}{1 + L} \leq \frac {\sigma_ {\tau} ^ {2}}{\sigma_ {\tau} ^ {2} + \sigma^ {2}} \leq \frac {U}{1 + U}\tag{3.57}
$$

To illustrate this procedure, we find a 95 percent confidence interval on $\sigma_{\tau}^{2}/(\sigma_{\tau}^{2}+\sigma^{2})$ for the strength data in Example 3.10. Recall that $MS_{\text{Treatments}} = 29.73, MS_{E} = 1.90, a = 4, n = 4, F_{0.025,3,12} = 4.47$ , and $F_{0.975,3,12} = 1/F_{0.025,12,3} = 1/14.34 = 0.070$ . Therefore, from Equation 3.56a and b,

$$
L = \frac {1}{4} \left[ \left(\frac {2 9 . 7 3}{1 . 9 0}\right) \left(\frac {1}{4 . 4 7}\right) - 1 \right] = 0. 6 2 5
$$

$$
U = \frac {1}{4} \left[ \left(\frac {2 9 . 7 3}{1 . 9 0}\right) \left(\frac {1}{0 . 0 7 0}\right) - 1 \right] = 5 5. 6 3 3
$$

and from Equation 3.57, the 95 percent confidence interval on $\sigma_{\tau}^{2}/(\sigma_{\tau}^{2}+\sigma^{2})$ is

$$
\frac {0 . 6 2 5}{1 . 6 2 5} \leq \frac {\sigma^ {2}}{\sigma_ {\tau} ^ {2} + \sigma^ {2}} \leq \frac {5 5 . 6 3 3}{5 6 . 6 3 3}
$$

or

$$
0. 3 8 \leq \frac {\sigma^ {2}}{\sigma_ {\tau} ^ {2} + \sigma^ {2}} \leq 0. 9 8
$$

We conclude that variability between looms accounts for between 38 and 98 percent of the variability in the observed strength of the fabric produced. This confidence interval is relatively wide because of the small number of looms used in the experiment. Clearly, however, the variability between looms ( $\sigma_{\tau}^{2}$ ) is not negligible.

Estimation of the Overall Mean $\mu$ . In many random effects experiments, the experimenter is interested in estimating the overall mean $\mu$ . From the basic model assumptions, it is easy to see that the expected value of any observation is just the overall mean. Consequently, an unbiased estimator of the overall mean is

$$
\hat {\mu} = \overline {{y}} _ {..}
$$

So for Example 3.10 the estimate of the overall mean strength is

$$
\hat {\mu} = \overline {{{y}}} _ {\cdot \cdot} = \frac {y _ {\cdot \cdot}}{N} = \frac {1 5 2 7}{1 6} = 9 5. 4 4
$$

It is also possible to find a $100(1-\alpha)\%$ confidence interval on the overall mean. The variance of $\bar{y}$ is

$$
V (\overline {{{y}}} _ {\cdot \cdot}) = V \left(\frac {\sum_ {i = 1} ^ {a} \sum_ {j = 1} ^ {n} y _ {i j}}{a n}\right) = \frac {n \sigma_ {\tau} ^ {2} + \sigma^ {2}}{a n}
$$

The numerator of this ratio is estimated by the treatment mean square, so an unbiased estimator of $V(\overline{y})$ is

$$
\hat {V} (\overline {{y}} _ {\cdot \cdot}) = \frac {M S _ {\text { Treatments }}}{a n}
$$

Therefore, the $100(1 - \alpha)\%$ CI on the overall mean is

$$
\overline {{{\mathrm{y}}}} _ {\cdot \cdot} - t _ {\alpha / 2, a (n - 1)} \sqrt {\frac {M S _ {\text { Treatments }}}{a n}} \leq \mu \leq \overline {{{\mathrm{y}}}} _ {\cdot \cdot} + t _ {\alpha / 2, a (n - 1)} \sqrt {\frac {M S _ {\text { Treatments }}}{a n}}\tag{3.58}
$$

To find a 95% CI on the overall mean in the fabric strength experiment from Example 3.10, we need $MS_{Treatments} = 29.73$ and $t_{0.025,12} = 2.18$ . The CI is computed from Equation 3.58 as follows:

$$
\overline {{{{\mathrm{y}}}}} _ {\cdot \cdot} - t _ {\alpha / 2, a (n - 1)} \sqrt {\frac {M S _ {\text { Treatments }}}{a n}} \leq \mu \leq \overline {{{{\mathrm{y}}}}} _ {\cdot \cdot} + t _ {\alpha / 2, a (n - 1)} \sqrt {\frac {M S _ {\text { Treatments }}}{a n}}
$$

$$
9 5. 4 4 - 2. 1 8 \sqrt {\frac {2 9 . 7 3}{1 6}} \leq \mu \leq 9 5. 4 4 + 2. 1 8 \sqrt {\frac {2 9 . 7 3}{1 6}}
$$

$$
9 2. 4 7 \leq \mu \leq 9 8. 4 1
$$

So, at 95 percent confidence the mean strength of the fabric produced by the looms in this facility is between 92.47 and 98.41. This is a relatively wide confidence interval because a small number of looms were sampled and there is a large difference between looms as reflected by the large portion of total variability that is accounted for by the differences between looms.

Maximum Likelihood Estimation of the Variance Components. Earlier in this section we presented the analysis of variance method of variance component estimation. This method is relatively straightforward to apply and makes use of familiar quantities—the mean squares in the analysis of variance table. However, the method has some disadvantages. As we pointed out previously, it is a method of moments estimator, a technique that mathematical statisticians generally do not prefer to use for parameter estimation because it often results in parameter estimates that do not have good statistical properties. One obvious problem is that it does not always lead to an easy way to construct confidence intervals on the variance components of interest. For example, in the single-factor random model, there is not a simple way to construct confidence intervals on $\sigma_{\tau}^{2}$ , which is certainly a parameter of primary interest to the experimenter. The preferred parameter estimation technique is called the method of maximum likelihood. The implementation of this method can be somewhat involved, particularly for an experimental design model, but it has been incorporated in some modern computer software packages that support designed experiments, including JMP.

A complete presentation of the method of maximum likelihood is beyond the scope of this book, but the general idea can be illustrated very easily. Suppose that x is a random variable with probability distribution $f(x,\theta)$ , where $\theta$ is an unknown parameter. Let $x_{1}, x_{2}, \ldots, x_{n}$ be a random sample of n observations. The joint probability distribution of the sample is $\prod_{i=1}^{n} f(x_{i}, \theta)$ . The likelihood function is just this joint probability distribution with the sample observations consider fixed and the parameter $\theta$ unknown. Note that the likelihood function, say

$$
L (x _ {1}, x _ {2}, \dots , x _ {n}; \theta) = \prod_ {i = 1} ^ {n} f (x _ {i}, \theta)
$$

is now a function of only the unknown parameter $\theta$ . The maximum likelihood estimator of $\theta$ is the value of $\theta$ that maximizes the likelihood function $L(x_{1}, x_{2}, \ldots, x_{n}; \theta)$ . To illustrate how this applies to an experimental design model with random effects, let y be the $an \times 1$ vector of observations for a single-factor random effects model with a treatments and n replicates and let $\sum$ be the $an \times an$ covariance matrix of the observations. Refer to Section 3.9.1 where we developed this covariance matrix for the special case where a = 3 and n = 2. The likelihood function is

$$
L \left(x _ {1 1}, x _ {1 2}, \dots , x _ {a, n}; \mu , \sigma_ {\tau} ^ {2}, \sigma^ {2}\right) = \frac {1}{(2 \pi) ^ {N / 2} [ \sum ] ^ {1 / 2}} \exp \left[ - \frac {1}{2} (\mathbf {y} - \mathbf {j} _ {N} \mu) ^ {\prime} \sum^ {- 1} (\mathbf {y} - \mathbf {j} _ {N} \mu) \right]
$$

where N = an is the total number of observations, $j_{N}$ is an $N \times 1$ vector of 1s, and $\mu$ is the overall mean in the model. The maximum likelihood estimates of the parameters $\mu$ , $\sigma_{\tau}^{2}$ , and $\sigma^{2}$ are the values of these quantities that maximize the likelihood function.

Maximum likelihood estimators (MLEs) have some very useful properties. For large samples, they are unbiased, and they have a normal distribution. Furthermore, the inverse of the matrix of second derivatives of the likelihood function (multiplied by -1) is the covariance matrix of the MLEs. This makes it relatively easy to obtain approximate confidence intervals on the MLEs.

The standard variant of maximum likelihood estimation that is used for estimating variance components is known as the residual maximum likelihood (REML) method. It is popular because it produces unbiased estimators and like all MLEs, it is easy to find CIs. The basic characteristic of REML is that it takes the location parameters in the model into account when estimating the random effects. As a simple example, suppose that we want to estimate the mean and variance of a normal distribution using the method of maximum likelihood. It is easy to show that the MLEs are

$$
\hat {\mu} = \frac {\sum_ {i = 1} ^ {n} y _ {i}}{n} = \overline {{y}}
$$

$$
\hat {\sigma} ^ {2} = \frac {\sum_ {i = 1} ^ {n} (y _ {i} - \overline {{{y}}}) ^ {2}}{n}
$$

Notice that the MLE $\hat{\sigma}^{2}$ is not the familiar sample standard deviation. It does not take the estimation of the location parameter $\mu$ into account. The REML estimator would be

$$
S ^ {2} = \frac {\sum_ {i = 1} ^ {n} (y _ {i} - \overline {{{{y}}}}) ^ {2}}{n - 1}
$$

The REML estimator is unbiased.

■ TABLE 3.20
JMP Output for the Loom Experiment in Example 3.10

<table><tr><td colspan="7">Response Y</td></tr><tr><td colspan="7">Summary of Fit</td></tr><tr><td>RSquare</td><td></td><td>0.793521</td><td></td><td></td><td></td><td></td></tr><tr><td>RSquare Adj</td><td></td><td>0.793521</td><td></td><td></td><td></td><td></td></tr><tr><td>Root Mean Square Error</td><td></td><td>1.376893</td><td></td><td></td><td></td><td></td></tr><tr><td>Mean of Response</td><td></td><td>95.4375</td><td></td><td></td><td></td><td></td></tr><tr><td>Observations (or Sum Wgts)</td><td></td><td>16</td><td></td><td></td><td></td><td></td></tr><tr><td colspan="7">Parameter Estimates</td></tr><tr><td>Term</td><td>Estimate</td><td>Std Error</td><td>DFDen</td><td>t Ratio</td><td colspan="2">Prob &gt; |t|</td></tr><tr><td>Intercept</td><td>95.4375</td><td>1.363111</td><td>3</td><td>70.01</td><td colspan="2">&lt; .0001*</td></tr><tr><td colspan="7">REML Variance Component Estimates</td></tr><tr><td>Random Effect</td><td>Var Ratio</td><td>Var Component</td><td>Std Error</td><td>95% Lower</td><td>95% Upper</td><td>Pct of Total</td></tr><tr><td>X1</td><td>3.6703297</td><td>6.9583333</td><td>6.0715247</td><td>-4.941636</td><td>18.858303</td><td>78.588</td></tr><tr><td>Residual</td><td></td><td>1.8958333</td><td>0.7739707</td><td>0.9748608</td><td>5.1660065</td><td>21.412</td></tr><tr><td>Total</td><td></td><td>8.8541667</td><td></td><td></td><td></td><td>100.000</td></tr><tr><td colspan="7">Covariance Matrix of Variance Component Estimates</td></tr><tr><td>Random Effect</td><td></td><td>X1</td><td>Residual</td><td></td><td></td><td></td></tr><tr><td>X1</td><td></td><td>36.863412</td><td>-0.149758</td><td></td><td></td><td></td></tr><tr><td>Residual</td><td></td><td>-0.149758</td><td>0.5990307</td><td></td><td></td><td></td></tr></table>

To illustrate the REML method, Table 3.20 presents the JMP output for the loom experiment in Example 3.10. The REML estimates of the model parameters $\mu$ , $\sigma_{\tau}^{2}$ , and $\sigma^{2}$ are shown in the output. Note that the REML estimates of the variance components are identical to those found earlier by the ANOVA method. These two methods will agree for balanced designs. However, the REML output also contains the covariance matrix of the variance components. The square roots of the main diagonal elements of this matrix are the standard errors of the variance components. If $\hat{\theta}$ is the MLE of $\theta$ and $\hat{\sigma}(\hat{\theta})$ is its estimated standard error, then the approximate $100(1 - \alpha)$ percent confidence interval on $\theta$ is

$$
\hat {\theta} - Z _ {\alpha / 2} \hat {\sigma} (\hat {\theta}) \leq \theta \leq \hat {\theta} + Z _ {\alpha / 2} \hat {\sigma} (\hat {\theta})
$$

JMP uses this approach to find the approximate CIs of $\sigma_{\tau}^{2}$ and $\sigma^{2}$ shown in the output. The 95 percent CI from REML for $\sigma^{2}$ is very similar to the chi-square-based interval computed earlier in Section 3.9.

## 3.10 The Regression Approach to the Analysis of Variance

We have given an intuitive or heuristic development of the analysis of variance. However, it is possible to give a more formal development. The method will be useful later in understanding the basis for the statistical analysis of more complex designs. Called the general regression significance test, the procedure essentially consists of finding the reduction in the total sum of squares for fitting the model with all parameters included and the reduction in sum of squares when the model is restricted to the null hypotheses. The difference between these two sums of squares is the treatment sum of squares with which a test of the null hypothesis can be conducted. The procedure requires the least squares estimators of the parameters in the analysis of variance model. We have given these parameter estimates previously (in Section 3.3.3); however, we now give a formal development.

## 3.10.1 Least Squares Estimation of the Model Parameters

We now develop estimators for the parameter in the single-factor ANOVA fixed-effects model

$$
y _ {i j} = \mu + \tau_ {i} + \epsilon_ {i j}
$$

using the method of least squares. To find the least squares estimators of $\mu$ and $\tau_{i}$ , we first form the sum of squares of the errors

$$
L = \sum_ {i = 1} ^ {a} \sum_ {j = 1} ^ {n} \epsilon_ {i j} ^ {2} = \sum_ {i = 1} ^ {a} \sum_ {j = 1} ^ {n} (y _ {i j} - \mu - \tau_ {i}) ^ {2}\tag{3.59}
$$

and then choose values of $\mu$ and $\tau_{i}$ , say $\hat{\mu}$ and $\hat{\tau}_{i}$ , that minimize L. The appropriate values would be the solutions to the $a + 1$ simultaneous equations

$$
\begin{array}{l} \frac {\partial L}{\partial \mu} \bigg | _ {\hat {\mu}, \hat {\tau} _ {i}} = 0 \\ \frac {\partial L}{\partial \tau_ {i}} \bigg | _ {\hat {\mu}, \hat {\tau} _ {i}} = 0 \quad i = 1, 2, \ldots , a \end{array}
$$

Differentiating Equation 3.59 with respect to $\mu$ and $\tau_{i}$ and equating to zero, we obtain

$$
- 2 \sum_ {i = 1} ^ {a} \sum_ {j = 1} ^ {n} (y _ {i j} - \hat {\mu} - \hat {\tau} _ {i}) = 0
$$

and

$$
- 2 \sum_ {j = 1} ^ {n} (y _ {i j} + \hat {\mu} - \hat {\tau} _ {i}) = 0 \quad i = 1, 2, \dots , a
$$

which, after simplification, yield

$$
\begin{array}{c c c} N \hat {\mu} + n \hat {\tau} _ {1} + n \hat {\tau} _ {2} + \dots + n \hat {\tau} _ {a} = y _ {..} \\ n \hat {\mu} + n \hat {\tau} _ {1} & = y _ {1}. \\ n \hat {\mu} + n \hat {\tau} _ {2} & = y _ {2}. \\ \vdots & \vdots \\ n \hat {\mu} & + n \hat {\tau} _ {a} = y _ {a}. \end{array}\tag{3.60}
$$

The $a + 1$ equations (Equation 3.60) in $a + 1$ unknowns are called the least squares normal equations. Notice that if we add the last a normal equations, we obtain the first normal equation. Therefore, the normal equations are not linearly independent, and no unique solution for $\mu, \tau_{i}, \ldots, \tau_{a}$ exists. This has happened because the effects model is overparameterized. This difficulty can be overcome by several methods. Because we have defined the treatment effects as deviations from the overall mean, it seems reasonable to apply the constraint

$$
\sum_ {i = 1} ^ {a} \hat {\tau} _ {i} = 0\tag{3.61}
$$

Using this constraint, we obtain as the solution to the normal equations

$$
\begin{array}{l} \hat {\mu} = \overline {{y}} _ {..} \\ \hat {\tau} _ {i} = \overline {{y}} _ {i.} - \overline {{y}} _ {..} \quad i = 1, 2, \ldots , a \end{array}\tag{3.62}
$$

This solution is obviously not unique and depends on the constraint (Equation 3.61) that we have chosen. At first this may seem unfortunate because two different experimenters could analyze the same data and obtain different results if they apply different constraints. However, certain functions of the model parameters are uniquely estimated, regardless of the constraint. Some examples are $\tau_{i} - \tau_{j}$ , which would be estimated by $\hat{\tau}_i - \hat{\tau}_j = \overline{y}_{i.} - \overline{y}_{j.}$ , and the $i$ th treatment mean $\mu_i = \mu + \tau_i$ , which would be estimated by $\hat{\mu}_i = \hat{\mu} + \hat{\tau}_i = \overline{y}_{i.}$ .

Because we are usually interested in differences among the treatment effects rather than their actual values, it causes no concern that the $\tau_{i}$ cannot be uniquely estimated. In general, any function of the model parameters that is a linear combination of the left-hand side of the normal equations (Equations 3.60) can be uniquely estimated. Functions that are uniquely estimated regardless of which constraint is used are called estimable functions. For more information, see the supplemental material for this chapter. We are now ready to use these parameter estimates in a general development of the analysis of variance.

## 3.10.2 The General Regression Significance Test

A fundamental part of this procedure is writing the normal equations for the model. These equations may always be obtained by forming the least squares function and differentiating it with respect to each unknown parameter, as we did in Section 3.9.1. However, an easier method is available. The following rules allow the normal equations for any experimental design model to be written directly:

RULE 1. There is one normal equation for each parameter in the model to be estimated.

RULE 2. The right-hand side of any normal equation is just the sum of all observations that contain the parameter associated with that particular normal equation.

To illustrate this rule, consider the single-factor model. The first normal equation is for the parameter $\mu$ ; therefore, the right-hand side is $y_{\infty}$ because all observations contain $\mu$ .

RULE 3. The left-hand side of any normal equation is the sum of all model parameters, where each parameter is multiplied by the number of times it appears in the total on the right-hand side. The parameters are written with a circumflex ( $^{\wedge}$ ) to indicate that they are estimators and not the true parameter values.

For example, consider the first normal equation in a single-factor experiment. According to the aforementioned rules, it would be

$$
N \hat {\mu} + n \hat {\tau} _ {1} + n \hat {\tau} _ {2} + \dots + n \hat {\tau} _ {a} = y _ {..}
$$

because $\mu$ appears in all N observations, $\tau_{1}$ appears only in the n observations taken under the first treatment, $\tau_{2}$ appears only in the n observations taken under the second treatment, and so on. From Equation 3.60, we verify that the equation shown above is correct. The second normal equation would correspond to $\tau_{1}$ and is

$$
n \hat {\mu} + n \hat {\tau} _ {1} = y _ {1}.
$$

because only the observations in the first treatment contain $\tau_{1}$ (this gives $y_{1}$ , as the right-hand side), $\mu$ and $\tau_{1}$ appear exactly n times in $y_{1}$ , and all other $\tau_{i}$ appear zero times. In general, the left-hand side of any normal equation is the expected value of the right-hand side.

Now, consider finding the reduction in the sum of squares by fitting a particular model to the data. By fitting a model to the data, we “explain” some of the variability; that is, we reduce the unexplained variability by some amount. The reduction in the unexplained variability is always the sum of the parameter estimates, each multiplied by the right-hand side of the normal equation that corresponds to that parameter. For example, in a single-factor experiment, the reduction due to fitting the full model $y_{ij} = \mu + \tau_i + \epsilon_{ij}$ is

$$
\begin{array}{r l} & {R (\mu , \tau) = \hat {\mu} y _ {..} + \hat {\tau} _ {1} y _ {1.} + \hat {\tau} _ {2} y _ {2.} + \dots + \hat {\tau} _ {a} y _ {a.}} \\ & {\qquad = \hat {\mu} y _ {..} + \sum_ {i = 1} ^ {a} \hat {\tau} _ {i} y _ {i.}} \end{array}\tag{3.63}
$$

The notation $R(\mu,\tau)$ means that reduction in the sum of squares from fitting the model containing $\mu$ and $\{\tau_{i}\}$ . $R(\mu,\tau)$ is also sometimes called the “regression” sum of squares for the full model $y_{ij} = \mu + \tau_{i} + \epsilon_{ij}$ . The number of degrees of freedom associated with a reduction in the sum of squares, such as $R(\mu, \tau)$ , is always equal to the number of linearly independent normal equations. The remaining variability unaccounted for by the model is found from

$$
S S _ {E} = \sum_ {i = 1} ^ {a} \sum_ {j = 1} ^ {n} y _ {i j} ^ {2} - R (\mu , \tau)\tag{3.64}
$$

This quantity is used in the denominator of the test statistic for $H_{0}:\tau_{1}=\tau_{2}=\ldots=\tau_{a}=0$ .

We now illustrate the general regression significance test for a single-factor experiment and show that it yields the usual one-way analysis of variance. The model is $y_{ij} = \mu + \tau_i + \epsilon_{ij}$ , and the normal equations are found from the above rules as

$$
\begin{array}{c c c} {N \hat {\mu} + n \hat {\tau} _ {1} + n \hat {\tau} _ {2} + \dots + n \hat {\tau} _ {a} = y _ {..}} \\ {n \hat {\mu} + n \hat {\tau} _ {1}} & & {= y _ {1}.} \\ {n \hat {\mu}} & {+ n \hat {\tau} _ {2}} & {= y _ {2}.} \\ & \vdots & \vdots \\ {n \hat {\mu}} & & {+ n \hat {\tau} _ {a} = y _ {a}.} \end{array}
$$

Compare these normal equations with those obtained in Equation 3.60.

Applying the constraint $\sum_{i=1}^{a} \hat{\tau}_i = 0$ , we find that the estimators for $\mu$ and $\tau_i$ are

$$
\hat {\mu} = \overline {{{y}}} _ {..} \quad \hat {\tau} _ {i} = \overline {{{y}}} _ {i.} - \overline {{{y}}} _ {..} i = 1, 2, \dots , a
$$

The reduction in the sum of squares due to fitting this full model is found from Equation 3.48 as

$$
\begin{array}{r l} R (\mu , \tau) & = \hat {\mu} y _ {..} + \sum_ {i = 1} ^ {a} \hat {\tau} _ {i} y _ {i.} \\ & = (\overline {{y}} _ {..}) y _ {..} + \sum_ {i = 1} ^ {a} (\overline {{y}} _ {i.} - \overline {{y}} _ {..}) y _ {i.} \\ & = \frac {y _ {. .} ^ {2}}{N} + \sum_ {i = 1} ^ {a} \overline {{y}} _ {i.} y _ {i.} - \overline {{y}} _ {..} \sum_ {i = 1} ^ {a} y _ {i.} \\ & = \sum_ {i = 1} ^ {a} \frac {y _ {i .} ^ {2}}{n} \end{array}
$$

which has $a$ degrees of freedom because there are $a$ linearly independent normal equations. The error sum of squares is, from Equation 3.64,

$$
\begin{array}{c} S S _ {E} = \sum_ {i = 1} ^ {a} \sum_ {j = 1} ^ {n} y _ {i j} ^ {2} - R (\mu , \tau) \\ = \sum_ {i = 1} ^ {a} \sum_ {j = 1} ^ {n} y _ {i j} ^ {2} - \sum_ {i = 1} ^ {a} \frac {y _ {i .} ^ {2}}{n} \end{array}
$$

and has $N - a$ degrees of freedom.

To find the sum of squares resulting from the treatment effects (the $\{\tau_{i}\}$ ), we consider a reduced model; that is, the model to be restricted to the null hypothesis ( $\tau_{i}=0$ for all i). The reduced model is $y_{ij}=\mu+\epsilon_{ij}$ . There is only one normal equation for this model:

$$
N \hat {\mu} = y _ {..}
$$

and the estimator of $\mu$ is $\hat{\mu} = \overline{y}_{\cdot}$ . Thus, the reduction in the sum of squares that results from fitting the reduced model containing only $\mu$ is

$$
R (\mu) = (\overline {{{y}}} _ {\cdot \cdot}) (y _ {\cdot \cdot}) = \frac {y _ {\cdot \cdot} ^ {2}}{N}
$$

Because there is only one normal equation for this reduced model, $R(\mu)$ has one degree of freedom. The sum of squares due to the $\{\tau_{i}\}$ , given that $\mu$ is already in the model, is the difference between $R(\mu,\tau)$ and $R(\mu)$ , which is

$$
\begin{array}{r l} R (\tau | \mu) & = R (\mu , \tau) - R (\mu) \\ & = R (\text { Full   Model }) - R (\text { Reduced   Model }) \\ & = \frac {1}{n} \sum_ {i = 1} ^ {a} y _ {i.} ^ {2} - \frac {y _ {. .} ^ {2}}{N} \end{array}
$$

with a-1 degrees of freedom, which we recognize from Equation 3.9 as $SS_{Treatments}$ . Making the usual normality assumption, we obtain the appropriate statistic for testing $H_{0}:\tau_{1}=\tau_{2}=\cdots=\tau_{a}=0$

$$
F _ {0} = \frac {R (\tau | \mu) (/ (a - 1)}{\left[ \sum_ {i = 1} ^ {a} \sum_ {j = 1} ^ {n} y _ {i j} ^ {2} - R (\mu , \tau) \right] / (N - a)}
$$

which is distributed as $F_{a-1,N-a}$ under the null hypothesis. This is, of course, the test statistic for the single-factor analysis of variance.

## 3.11 Nonparametric Methods in the Analysis of Variance

## 3.11.1 The Kruskal-Wallis Test

In situations where the normality assumption is unjustified, the experimenter may wish to use an alternative procedure to the F-test analysis of variance that does not depend on this assumption. Such a procedure has been developed by Kruskal and Wallis (1952). The Kruskal–Wallis test is used to test the null hypothesis that the a treatments are identical against the alternative hypothesis that some of the treatments generate observations that are larger than others. Because the procedure is designed to be sensitive for testing differences in means, it is sometimes convenient to think of the Kruskal–Wallis test as a test for equality of treatment means. The Kruskal–Wallis test is a nonparametric alternative to the usual analysis of variance.

To perform a Kruskal–Wallis test, first rank the observations $y_{ij}$ in ascending order and replace each observation by its rank, say $R_{ij}$ , with the smallest observation having rank 1. In the case of ties (observations having the same value), assign the average rank to each of the tied observations. Let $R_{i}$ be the sum of the ranks in the ith treatment. The test statistic is

$$
H = \frac {1}{S ^ {2}} \left[ \sum_ {i = 1} ^ {a} \frac {R _ {i .} ^ {2}}{n _ {i}} - \frac {N (N + 1) ^ {2}}{4} \right]\tag{3.65}
$$

where $n_{i}$ is the number of observations in the ith treatment, N is the total number of observations, and

$$
S ^ {2} = \frac {1}{N - 1} \left[ \sum_ {i = 1} ^ {a} \sum_ {j = 1} ^ {n _ {i}} R _ {i j} ^ {2} - \frac {N (N + 1) ^ {2}}{4} \right]\tag{3.66}
$$

Note that $S^2$ is just the variance of the ranks. If there are no ties, $S^2 = N(N + 1) / 12$ and the test statistic simplifies to

$$
H = \frac {1 2}{N (N + 1)} \sum_ {i = 1} ^ {a} \frac {R _ {i .} ^ {2}}{n _ {i}} - 3 (N + 1)\tag{3.67}
$$

When the number of ties is moderate, there will be little difference between Equations 3.66 and 3.67, and the simpler form (Equation 3.67) may be used. If the $n_i$ are reasonably large, say $n_i \geq 5, H$ is distributed approximately as $\chi_{a-1}^2$ under the null hypothesis. Therefore, if

$$
H > \chi_ {\alpha , a - 1} ^ {2}
$$

the null hypothesis is rejected. The P-value approach could also be used.

## EXAMPLE 3.11

The data from Example 3.1 and their corresponding ranks are shown in Table 3.21. There are ties, so we use Equation 3.65 as the test statistic. From Equation 3.65

$$
S ^ {2} = \frac {1}{1 9} \left[ 2 8 6 9. 5 0 - \frac {2 0 (2 1) ^ {2}}{4} \right] = 3 4. 9 7
$$

and the test statistic is

$$
\begin{array}{r l} H & = \frac {1}{S ^ {2}} \left[ \sum_ {i = 1} ^ {a} \frac {R _ {i .} ^ {2}}{n _ {i}} - \frac {N (N + 1) ^ {2}}{4} \right] \\ & = \frac {1}{3 4 . 9 7} [ 2 7 9 6. 3 0 - 2 2 0 5 ] \\ & = 1 6. 9 1 \end{array}
$$

## TABLE 3.21

Data and Ranks for the Plasma Etching Experiment in Example 3.1

<table><tr><td colspan="8">Power</td></tr><tr><td colspan="2">160</td><td colspan="2">180</td><td colspan="2">200</td><td colspan="2">220</td></tr><tr><td> $y_{1j}$ </td><td> $R_{1j}$ </td><td> $y_{2j}$ </td><td> $R_{2j}$ </td><td> $y_{3j}$ </td><td> $R_{3j}$ </td><td> $y_{4j}$ </td><td> $R_{4j}$ </td></tr><tr><td>575</td><td>6</td><td>565</td><td>4</td><td>600</td><td>10</td><td>725</td><td>20</td></tr><tr><td>542</td><td>3</td><td>593</td><td>9</td><td>651</td><td>15</td><td>700</td><td>17</td></tr><tr><td>530</td><td>1</td><td>590</td><td>8</td><td>610</td><td>11.5</td><td>715</td><td>19</td></tr><tr><td>539</td><td>2</td><td>579</td><td>7</td><td>637</td><td>14</td><td>685</td><td>16</td></tr><tr><td>570</td><td>5</td><td>610</td><td>11.5</td><td>629</td><td>13</td><td>710</td><td>18</td></tr><tr><td> $R_{i.}$ </td><td>17</td><td></td><td>39.5</td><td></td><td>63.5</td><td></td><td>90</td></tr></table>

Because $H > \chi_{0.01,3}^{2} = 11.34$ , we would reject the null hypothesis and conclude that the treatments differ. (The

P-value for H = 16.91 is $P = 7.38 \times 10^{-4}$ . This is the same conclusion as given by the usual analysis of variance F-test.

## 3.11.2 General Comments on the Rank Transformation

The procedure used in the previous section of replacing the observations by their ranks is called the rank transformation. It is a very powerful and widely useful technique. If we were to apply the ordinary F-test to the ranks rather than to the original data, we would obtain

$$
F _ {0} = \frac {H / (a - 1)}{(N - 1 - H) / (N - a)}\tag{3.68}
$$

as the test statistic [see Conover (1980), p. 337]. Note that as the Kruskal–Wallis statistic H increases or decreases, $F_{0}$ also increases or decreases, so the Kruskal–Wallis test is equivalent to applying the usual analysis of variance to the ranks.

The rank transformation has wide applicability in experimental design problems for which no nonparametric alternative to the analysis of variance exists. This includes many of the designs in subsequent chapters of this book. If the data are ranked and the ordinary F-test is applied, an approximate procedure that has good statistical properties results [see Conover and Iman (1976, 1981)]. When we are concerned about the normality assumption or the effect of outliers or “wild” values, we recommend that the usual analysis of variance be performed on both the original data and the ranks. When both procedures give similar results, the analysis of variance assumptions are probably satisfied reasonably well, and the standard analysis is satisfactory. When the two procedures differ, the rank transformation should be preferred because it is less likely to be distorted by nonnormality and unusual observations. In such cases, the experimenter may want to investigate the use of transformations for nonnormality and examine the data and the experimental procedure to determine whether outliers are present and why they have occurred.