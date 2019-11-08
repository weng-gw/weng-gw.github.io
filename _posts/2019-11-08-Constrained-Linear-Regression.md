---
layout: post
title:  "On general linear hypothesis test"
date:  2019-11-08
categories: statistics 
---

I am teaching an undergraduate course on statistical computing this
semester. When I was preparing for the course slides on simulation study for
general linear hypothesis test, I found that some of the optimization problems
in this topic were talked in very detail in most statistics (I didn't think
about the optimization perspect of the problem when I took the course in my
first year of PhD study). So here, I want to write a thorough note on this problem
to provide a detailed derivation  on things related to this test.

## Set-up 

Consider the following learning regression model:

$$Y=X\beta +\epsilon$$,


where $$Y=(Y_1,\ldots,Y_n)^T$$ is the vector of responses, $$X$$ is a
$$n\times p$$ design matrix and $$\epsilon=(\epsilon_1,\ldots, \epsilon_n)^T$$
is the vector containing the random errors. Moreover, we assume the random
errors are iid $$N(0,\sigma^2)$$. In the problem of general liner hypothesis
test, we are interested in testing the hypotheses $$H_0: A\beta=b$$ against
$$H_a:A\beta\ne b$$. 

## Least-squares estimator and constrained least-squares estimator

The test statistic is constructed based on the residual sum
of squares using the unconstrained least-squares estimator and the residual sum
of squares using the constrained least-squares. The unconstrained least-squares
estimator can be calculated as $$\hat{\beta}=(X^TX)^{-1}X^TY$$ provided $$X^TX$$
is invertible. The corresponding residual sum of squares is
$$RSS_f=(Y-X\hat{\beta})^T(Y-X\hat{\beta})$$.

For the constrained problem, the least-squares estimator $$\tilde{\beta}$$ is
the solution to the following problem

$$\begin{eqnarray}\text{minimize}& (Y-X\beta)^T(Y-X\beta)\\
s.t& A\beta=b
\end{eqnarray}
$$

This is an quadratic optimization problem with linear equality constraints. This
is usually not talked much in a typical statistical course. And my origin motivation
for this post is to provide the solution to this problem in detail. We will use
KKT condition to solve this problem and the solution could be found in section
10.1.1 of [*Convex Optimization*](https://web.stanford.edu/~boyd/cvxbook/) by
Stephen Boyd  and Lieven Vandenberghe.

The Lagrangian functio, with Lagrange parameter $$\lambda$$ is 

$$L(\beta,\lambda)=(Y-X\beta)^T(Y-X\beta)+\lambda^T(A\beta-b)$$

By KKT condition, the optimality conditions are  $$\frac{\partial L}{\partial
\beta}(\tilde{\beta},\lambda)=0$$ and $$\frac{\partial L}{\partial
\lambda}(\tilde{\beta},\lambda)=0$$. So we will have the linear equation system
$$2X^TX\tilde{\beta}-2X^TY+A^T\lambda=0$$ and $$A\tilde{\beta}-b=0$$. Or in
matrix form, this could be written as

$$
\left[\begin{array}{cc}
2X^TX&A^T\\
A&0
\end{array}\right]
\left[\begin{array}{c}
\tilde{\beta}\\
\lambda
\end{array}\right]=
\left[\begin{array}{c}
2X^TY\\
b
\end{array}\right]
$$

The coefficient matrix is called the KKT matrix and is invertible if and only if
$$A$$ has independent rows and $$[X^T,A^T]^T$$ has independent columns. Then by
block matrix inverse formula we have the inverse of KKT matrix can be written as 

$$
\begin{eqnarray}
\left[\begin{array}{cc}
\frac{1}{2}(X^TX)^{-1}+\frac{1}{2}(X^TX)^{-1}A^T(A(X^TX)^{-1}A^T)^{-1}A(X^TX)^{-1}&-(X^TX)^{-1}A^T(A(X^TX)^{-1}A^T)^{-1}\\
-(A(X^TX)^{-1}A^T)^{-1}A(X^TX)^{-1}&(A(X^TX)^{-1}A^T)^{-1}
\end{array}\right]
\end{eqnarray}
$$

Then by some simple calculation we can show 

$$\tilde{\beta}=\hat{\beta}-(X^TX)^{-1}A^T(A(X^TX)^{-1}A^T)^{_1}(A\hat{\beta}-b)$$

And the corresponding residual sum of squares is
$$RSS_0=(Y-X\tilde{\beta})^T(Y-X\tilde{\beta})$$.

## Test statistic
The test statistic for the test is 

$$
F=\frac{(RSS_0-RSS_f)/d}{RSS_f/(n-p)}
$$,

Basically, it measures the difference of goodness of fit between the full model
and the null model and the intuition is to reject the null when $$F$$ is large. 

By plugging in the expressions for $\hat{\beta}$ and $\tilde{\beta}$, we have

$$RSS_0-RSS_f=(A\hat{\beta}-b)^T(A(X^TX)^{-1}A^T)^{-1}(A\hat{\beta}-b)^T$$


So $$F$$ is simplied to 

$$F=\frac{n-p}{d(Y-X\hat{\beta})^T(Y-X\hat{\beta})}(A\hat{\beta}-b)^T(A(X^TX)^{-1}A^T)^{-1}(A\hat{\beta}-b)^T$$.

We know under normality assumption and the null hypothesis

$$A\hat{\beta}-b\sim N(0,\sigma^2A(X^TX)^{-1}A^T)$$

Thus $$RSS_0-RSS_f$$ is $$\sigma^2\chi^2_d$$, while $$RSS_f$$ is
$$\sigma^2\chi^2_{(n-p)}$$. Then by the independence between $\hat{\beta}$ and
$$RSS_f$$, we have $$F\sim F_{d,n-p}$$ under the null.
