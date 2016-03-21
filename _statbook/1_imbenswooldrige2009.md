---
layout: mathpost
title: "Reading Imbens and Wooldrige 2009"
author: "Anh Le"
description: "This is my reading of it"
header-img: "img/about-bg.jpg"
comments: true
---

Reference: Imbens & Wooldrige (2009). Recent developments in the econometrics of program evaluation.[Open-access link](http://dash.harvard.edu/bitstream/handle/1/3043416/imbens_recent.pdf?sequence=2)

## Section 5.3: Using regression method to estimate average treatment effect.

> To describe the general approach to regression methods for estimating average treatment effects, define $$\mu_0(x)$$ and $$\mu_1(x)$$ to be the two regression functions for the potential outcomes.
> $$\mu_0(x) = E[Y_i(0) \vert X_i=x] = E[Y_i(0) \vert W_i=0, X_i=x] = \alpha_0 + \beta_0 ' (x - \psi_X)$$
> $$\mu_1(x) = E[Y_i(1) \vert X_i=x] = E[Y_i(1) \vert W_i=1, X_i=x] =\alpha_1 + \beta_1 ' (x - \psi_X)$$

**Comment:** The first equality is the definition. The second equality is true due to unconfoundedness. The third equality is because we're assuming the simplest case in which the conditional mean is linear in the parameters.

Note that $$x, \psi_X$$ refer to the covariates of the entire population, not just of the treated or untreated. This is because $$\mu_0, \mu_1$$ refer to the two potential outcomes, and all observations have both potential outcomes.

In other words, for each observation $$i$$, we can plug in $$i$$'s covariates $$x$$ into $$\mu_0(x)$$ and $$\mu_1(x)$$ in order to get the untreated and treated potential outcomes for observation $$i$$.

> Then $$\hat\tau_{reg}$$ is simply
> $$\hat\tau_{reg} = \hat{\alpha}_1 - \hat{\alpha}_0 \qquad (13)$$

**Comment**: The authors skip some steps here, which I'll elaborate as follows.

$$
\begin{align}
\tau_{reg} &= E[Y_i(1) - Y_i(0)] & \text{The average treatment effect} \\
&= E[E[Y_i(1) - Y_i(0) \vert X_i=x]] & \text{Law of Iterated Expectation} \\
&= E[E[Y_i(1) \vert X_i=x] - E[Y_i(0) \vert X_i=x]] \\
&= E[E[Y_i(1) \vert W_i=1, X_i=x] - E[Y_i(0) \vert W_i=0,X_i=x]] & \text{Due to unconfoundedness, as explained on p 23} \\
&= E[\mu_1(x) - \mu_0(x)] 
\end{align}
$$

Notice that we can't estimate $$E[Y_i(1) \vert X_i=x]$$ since $$Y_i(1)$$ is not available for the untreated $$i$$'s. However, $$E[Y_i(1) \vert W_i=1, X_i=x]$$ *can* be estimated, so $$\hat\mu_1(x)$$ exists. In order to estimate $$\hat\mu_0(x), \hat\mu_1(x)$$, we give them a parametric model like the authors do (p. 24):

$$
\begin{align}
\mu_0(x) &= \alpha_0 + \beta_0 ' (x - \psi_X) \\
\mu_1(x) &= \alpha_1 + \beta_1 ' (x - \psi_X)
\end{align}
$$

and just replacing all the parameters with its estimates (we'll talk about what these estimates actually are):

$$
\begin{align}
\hat\mu_0(x) &= \hat\alpha_0 + \hat\beta_0 ' (x - \bar X) \\
\hat\mu_1(x) &= \hat\alpha_1 + \hat\beta_1 ' (x - \bar X)
\end{align}
$$

But what is $$\hat\mu_0(x)$$? Recall that $$\mu_0(x) = E[Y_i(0) \vert W_i=0, X_i=x]$$. Since we observe $$Y_i(0)$$ for all the untreated observations, $$\hat\mu_0(x_i)$$ is simply $$Y_i$$!

To derive Eq. 13, just plug in the formula for $$\hat\mu_0(x), \hat\mu_1(x)$$.

$$
\begin{align}
\hat\tau_{reg} &= \frac{1}{N} \sum_{i=1}^{N} \left( \hat\mu_1(X_i) - \hat\mu_0(X_i) \right) \\
&= \frac{1}{N} \sum_{i=1}^{N} \left( \hat\alpha_1 + \hat\beta_1 ' (x - \bar X) - \hat\alpha_0 - \hat\beta_0 ' (x - \bar X) \right) \\
&= \hat\alpha_1 - \hat\alpha_0 + \hat\beta_1' \frac{1}{N} \sum_{i=1}^{N} (x - \bar X) - \hat\beta_0 ' \frac{1}{N} \sum_{i=1}^{N} (x - \bar X) \\
&= \hat\alpha_1 - \hat\alpha_0 + 0 + 0 \\
&= \hat\alpha_1 - \hat\alpha_0 & \text{Arrived at their Eq. 13}
\end{align}
$$

> A different representation of $$\hat\tau_{reg}$$ is useful in order to illustrate some of the concerns with regression estimators in this setting. Suppose we do use the linear model in (12). It can be shown that
> 
> $$\hat\tau_{reg} = \bar Y_1 - \bar Y_0 - \left( \frac{N_0}{N_0 + N_1} \hat\beta_1 + \frac{N_1}{N_0 + N_1} \hat\beta_0 \right)' (\bar X_1 - \bar X_0) \qquad (14)$$

**Comment**: To derive Eq. 14, we have to substitute out $$\hat\alpha_0, \hat\alpha_1$$ with the other terms. Notice that:

$$\hat\alpha_0 = \hat\mu_0(x) - \hat\beta_0'(x - \bar X)$$

But what is $$\hat\mu_0(x)$$? Recall that $$\mu_0(x) = E[Y_i(0) \vert W_i=0, X_i=x]$$. Since we observe $$Y_i(0)$$ for all the untreated observations, $$\hat\mu_0(x_i)$$ is simply $$Y_i$$!

Therefore,

$$
\begin{align}
\hat\mu_0(x) &= \hat\alpha_0 + \hat\beta_0 ' (x - \bar X) \\
\frac{1}{N_0} \sum_{i=1}^{N_0} \hat\mu_0(x_i) &= \hat\alpha_0 + \hat\beta_0 ' \frac{1}{N_0} \sum_{i=1}^{N_0} (x_i - \bar X) & \text{summing over untreated units}\\
\frac{1}{N_0} \sum_{i=1}^{N_0} Y_i &= \hat\alpha_0 + \hat\beta_0 ' \frac{1}{N_0} \sum_{i=1}^{N_0} (x_i - \bar X) & \text{Using the fact above, i.e.} \hat\mu_0(x_i) = Y_i\\
\bar Y_0 &= \hat\alpha_0 + \hat\beta_0 '(\bar X_0 - \bar X) \\
\hat\alpha_0 &= \bar Y_0 - \hat\beta_0 '(\bar X_0 - \bar X)
\end{align}
$$

We plug the above formula for $$\hat\alpha_0, \hat\alpha_1$$ into $$\hat\tau_{reg}$$ and do some transformation to get Eq. 14. Specifically,

$$
\begin{align}
\hat\alpha_0 &= \bar Y_0 - \hat\beta_0 '(\bar X_0 - \bar X) \\
&= \bar Y_0 - \hat\beta_0 '(\bar X_0 - \frac{\bar X_0 N_0 + \bar X_1 N_1}{N_0 + N_1}) \\
&= \bar Y_0 - \hat\beta_0 ' \frac{N_1}{N_0+N_1} \left( \bar X_0 \frac{N_0 + N_1}{N_1} - \bar X_0 \frac{N_0}{N_1} - \bar X_1 \right) \\
&= \bar Y_0 - \hat\beta_0 ' \frac{N_1}{N_0+N_1} (\bar X_0 - \bar X_1) \\
\hat\alpha_1 &= \bar Y_1 - \hat\beta_1 ' \frac{N_0}{N_0+N_1} (\bar X_1 - \bar X_0) & \text{similar to above}
\end{align}
$$

Plugging $$\hat\alpha_0, \hat\alpha_1$$ into $$\hat\tau_{reg}$$, we now have Eq. 14