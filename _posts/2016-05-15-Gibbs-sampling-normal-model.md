---
title: "Gibbs sampling with Rcpp"
author: "Anh Le"
layout: post
toc: yes
comments: true
---

When Gibbs sampling in R is slow, we can use Rcpp to speed it up by integrating C++ code into R. However, we cannot run Rcpp code line-by-line because C++ code needs to be compiled before being run. Rcpp code is thus harder to debug and it's paramount that we ensure the Rcpp code is doing what we wants it to.

The solution is to start with small Rcpp programs and make sure that they produce the same results as their pure R counterpart. Following this advice, I'm writing a series of increasingly complex Markov chain Monte Carlo (MCMC) code in Rcpp, eventually building towards my final goal of a custom Metropolis-Hastings algorithm for my work on modeling FDI flow.

The first installment will be a Gibbs sampler of a semi-conjugate normal model, arguably the most common introductory model to Bayesian statistics and MCMC. (Both [Hoff's First Course in Bayesian](http://www.stat.washington.edu/people/pdhoff/book.php) and [Gelman's Bayesian Data Analysis](http://www.stat.columbia.edu/~gelman/book/) cover it.)

I'll follow Hoff's example in Chapter 6 so that I have a ground truth to check whether my code is working correctly. For Rcpp basics, I read [Hadley Wickham's introduction to Rcpp](http://adv-r.had.co.nz/Rcpp.html). I'll note additional pitfalls along the way of converting R to Rcpp.

# Normal model with semi-conjugate priors

A normal model with semi-conjugate priors has the following specifications. For more details and intuition behind the model, see Chapter 5 & 6 in [Hoff's book](http://www.stat.washington.edu/people/pdhoff/book.php).

## Likelihood

$$
Data = Y_1, \dots, Y_n \sim i.i.d. N(\theta, \sigma^2)
$$

## Priors

$$
\begin{align}
p(\theta) &\sim N(\mu_0, \tau^2_0) \\
p(\sigma^2) &\sim Inverse-Gamma(\nu_0 / 2, \nu_0 \sigma^2_0 / 2)
\end{align}
$$

## Full conditionals

$$
\begin{align}
p(\theta | \tilde \sigma^2, Data) &= N(\mu_n, \tau^2_n) \\
p(\tilde \sigma^2 | \theta, Data) &= Gamma(\frac{\nu_n}{2}, \frac{\nu_n \sigma_n^2(\theta)}{2})
\end{align}
$$

with

$$
\begin{align}
\tau_n^2 = \frac{1}{\frac{1}{\tau_0^2} + \frac{n}{\sigma^2}} \qquad &\text{and} \qquad
\mu_n = \frac{\frac{\mu_0}{\tau_0^2} + \frac{n\bar y}{\sigma^2}}{\frac{1}{\tau_0^2} + \frac{n}{\sigma^2}} \\
\nu_n = \nu_0 + n \qquad &\text{and} \qquad \sigma^2_n(\theta) = \frac{1}{\nu_n} \left[ \nu_0\sigma_0^2 + n \frac{\sum (y_i - \theta)^2}{n}\right]
\end{align}
$$

There is a quick intuitive interpretation of the full conditional of $$\theta$$,
$$p(\theta | \tilde\sigma^2, Data)$$. Its posterior mean $$\mu_n$$ is the weighted sum of the prior mean $$\mu_0$$ and the sample mean $$\bar y$$. This makes sense if we think of the posterior as combining the information from the prior and the sample.

But how should we weigh the information from the prior and the sample? Intuitively, if there's more information in the prior, we should weigh it more. The amount of information is captured in the prior precision, $$\frac{1}{\tau_0^2}$$. The precision is the inverse of variance, because the smaller the variance around our prior belief about $$\theta$$, the more information we have in the prior. 

In sum, because precision captures the amount of information available, we should weigh the prior mean $$\mu_0$$ by the prior precision $$\frac{1}{\tau_0^2}$$. That is indeed the weight of the prior mean in the weighted sum.

Similarly, we weight the sample mean $$\bar y$$ by the sample precision $$\frac{n}{\sigma^2}$$. Note that we multiply the precision by $$n$$ to capture the fact that the bigger the sample size, the more information there is in the sample.

# Pure R Gibbs sampler

Below is a direct translation of the full conditionals below into a Gibbs sampler. Inside the big loop, for the iteration $$s$$ we sample:

1. $$p(\theta^{(s)} | \tilde \sigma^{2(s - 1)}, Data)$$
2. $$p(\tilde \sigma^{2(s)} | \theta^{(s)}, Data)$$

and repeats.


{% highlight r %}
f_gibbsR <- function(y, mu_0, tau2_0, sigma2_0, nu_0, S = 1000) {
  # Sample statistics
  n <- length(y) ; ybar <- mean(y)

  # Initialize storage
  res_theta <- rep(NA, S)
  res_sigma2 <- rep(NA, S)

  # Starting value as the sample variance and mean
  res_sigma2[1] <- var(y)
  res_theta[1] <- ybar
  
  # Gibbs sampler -- big loop
  for (s in 2:S) {
    # Full conditional for theta
    sigma2 <- res_sigma2[s - 1] # Get the current value of sigma2

    tau2_n <- 1 / (1 / tau2_0 + n / sigma2)
    mu_n <- (mu_0 / tau2_0 + n * ybar / sigma2) / (1 / tau2_0 + n / sigma2)

    # Generate and save theta. Note that rnorm takes SD, not variance
    res_theta[s] <- rnorm(1, mu_n, sqrt(tau2_n))

    # Full conditional for sigma
    theta <- res_theta[s] # Get the current value of theta

    nu_n <- nu_0 + n
    nu_sigma2_n <- nu_0 * sigma2_0 + sum((y - theta) ** 2)
    res_sigma2[s] <- 1 / rgamma(1, nu_n / 2, nu_sigma2_n / 2)
  }

  return(list(theta = res_theta, sigma2 = res_sigma2))
}
{% endhighlight %}

We replicate Hoff's example with the following


{% highlight r %}
y <- c(1.64,1.70,1.72,1.74,1.82,1.82,1.82,1.90,2.08)
set.seed(1)
gibbsR <- f_gibbsR(y = y, S = 1000,
                   mu_0 = 1.9, tau2_0 = 0.95 ** 2,
                   sigma2_0 = 0.01, nu_0 = 1)
quantile(gibbsR$theta, c(0.025, 0.5, 0.975))
{% endhighlight %}



{% highlight text %}
##     2.5%      50%    97.5% 
## 1.707282 1.804348 1.901129
{% endhighlight %}


{% highlight r %}
# set.seed(1)
# gibbsR <- f_gibbsR(y = y, S = 1000,
#                    mu_0 = 1.9, tau2_0 = 0.95 ** 2,
#                    sigma2_0 = 0.01, nu_0 = 1)
# set.seed(1)
# Rcpp::sourceCpp("_source/gibbsC.cpp")
# gibbsC <- f_gibbsC(y = y, S = 1000,
#                    mu_0 = 1.9, tau2_0 = 0.95 ** 2,
#                    sigma2_0 = 0.01, nu_0 = 1)
# 
# 
# quantile(gibbsR$theta, c(0.025, 0.5, 0.975))
# quantile(gibbsC$theta, c(0.025, 0.5, 0.975))
# quantile(sqrt(gibbsR$sigma2), c(0.025, 0.5, 0.975))
# quantile(sqrt(gibbsC$sigma2), c(0.025, 0.5, 0.975))
# 
# library(microbenchmark)
# 
# microbenchmark(f_gibbsR(y = y, S = 1000,
#                         mu_0 = 1.9, tau2_0 = 0.95 ** 2,
#                         sigma2_0 = 0.01, nu_0 = 1))
# 
{% endhighlight %}

