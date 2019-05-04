---
title: Small-Sample Odds Ratios in C
author: Yeng Miller-Chang
date: "2019-05-04"
draft: true
---

*This was the final project that I did for STAT 580: Statistical Computing in the Fall 2019 semester, taught by Dr. Karin Dorman at Iowa State University. Within one semester, the goal of the course was to learn how to perform statistical computation from C, starting from only knowing a bit about R, and nothing about C. If you ever get the opportunity to take this course, I strongly recommend it. It has been one of the best learning experiences I've ever had.*

# Mathematical Exposition

Mathematically, my project can be described as follows:

<center>
Given a ``$2 \times 2$`` contingency table of counts with ``$\alpha \in (0, 1)$`` fixed, generate a ``$100(1-\alpha)\%$`` exact confidence interval for the odds ratio for the table.
</center>

This problem, although seemingly simple at its surface, requires a substantial amount of computational work.

Suppose we have a ``$2 \times 2$`` contingency table of proportions which constitute a probability distribution:
``\begin{equation*}
\begin{array}{|c|c|}
\hline
\pi_{11} & \pi_{12} \\
\hline
\pi_{21} & \pi_{22} \\
\hline
\end{array}
\end{equation*}``

The odds ratio is defined by ``$\theta = \dfrac{\pi_{11}\pi_{22}}{\pi_{12}\pi_{21}}$``. ^[Agresti, A. (2013), *Categorical Data Analysis*  (3rd ed.), Hoboken, NJ: John Wiley \& Sons, p. 606.]

Consider a ``$2 \times 2$`` contingency table of counts derived from a sample:
``\begin{equation*}
\begin{array}{|c|c|}
\hline
n_{11} & n_{12} \\
\hline
n_{21} & n_{22} \\
\hline
\end{array}
\end{equation*}``

Let ``$+$`` denote summation over indices; for example, ``$n_{+1} = \sum_{i=1}^{2}n_{i1}$``, and ``$n_{++} = \sum_{i=1}^{2}\sum_{j=1}^{2}n_{ij}$``. The sample odds ratio is defined by^[Agresti, A. (2013), *Categorical Data Analysis*  (3rd ed.), Hoboken, NJ: John Wiley \& Sons, p. 69.]
``\begin{equation*}
\hat{\theta} = \dfrac{n_{11}n_{22}}{n_{21}n_{22}}\text{.}
\end{equation*}``
It can be shown that, using a multinomial assumption, the Central Limit Theorem, and the Delta Method, that an approximate standard error for ``$\log(\hat{\theta})$`` is^[Agresti, A. (2013), *Categorical Data Analysis*  (3rd ed.), Hoboken, NJ: John Wiley \& Sons, pp. 70-75.]
``\begin{equation*}
\hat{\sigma}_{\log(\hat{\theta})} = \sqrt{\dfrac{1}{n_{11}} + \dfrac{1}{n_{12}} + \dfrac{1}{n_{21}} + \dfrac{1}{n_{22}}}\text{.}
\end{equation*}``
One could easily use the approximation above to generate ``$100(1-\alpha)\%$`` confidence intervals for ``$\log(\theta)$`` based on a normal approximation; however, this is only appropriate for large samples. Suppose that the above table has been stratified by a variable (say, for example, age), for which each value of said variable has its own ``$2 \times 2$`` contingency table. As an example, suppose we have the following ``$2 \times 2$`` contingency table:
``\begin{equation*}
\begin{array}{|c|c|}
\hline
30 & 20\\
\hline
40 & 30 \\
\hline
\end{array}
\end{equation*}``
and that we decide to stratify these data based on another factor:
``\begin{align*}
\text{Group A} \qquad& \begin{array}{|c|c|}
\hline
20 & 15\\
\hline
10 & 15 \\
\hline
\end{array} \\
\text{Group B} \qquad& \begin{array}{|c|c|}
\hline
10 & 5\\
\hline
30 & 15 \\
\hline
\end{array}
\end{align*}``
In cases such as the one above, it does not seem reasonable to use a normal approximation. Set ``$n = n_{++}$``. It turns out that with the marginal totals given, it can be shown that ``$n_{11}$`` has probability mass function
``\begin{equation}
f(t \mid n_{1+}, n_{+1}, n) = \dfrac{\displaystyle\binom{n_{1+}}{t} \binom{n-n_{1+}}{n_{+1} - t} \theta^t}{\displaystyle\sum_{u=m_{-}}^{m_{+}}\binom{n_{1+}}{u} \binom{n-n_{1+}}{n_{+1} - u} \theta^u}\label{hypergeom-pmf}
\end{equation}``
for ``$m_{-} \leq t \leq m_{+}$`` with ``$m_{-} = \max(0, n_{1+} + n_{+1} - n)$`` and ``$m_{+} = \min(n_{1+}, n_{+1})$``, which is of the "noncentral hypergeometric distribution." Cornfield (1956)^[Cornfield, J. (1956), "A Statistical Problem Arising from Retrospective Studies," *Proceedings of the Third Berkeley Symposium on Mathematical Statistics and Probability*, 4, 135--148,  Berkeley, CA: University of California Press. Available at https://projecteuclid.org/euclid.bsmsp/1200502552.] provides one method to find a ``$100(1-\alpha)\%$`` confidence interval for ``$\theta$``: one could solve for ``$\theta_0$`` and ``$\theta_1$`` in the equations
``\begin{equation}
\dfrac{\alpha}{2} = \sum\limits_{t \geq n_{11}}f(t \mid n_{1+}, n_{+1}, n, \theta_1) = \sum\limits_{t \leq n_{11}}f(t \mid n_{1+}, n_{+1}, n, \theta_0)\text{.}\label{hypergeom-CI}
\end{equation}``

# Computational Exposition

The crux of this problem is to find the roots of the following functions of ``$\theta$``:
``\begin{align}
\sum\limits_{t \geq n_{11}}f(t \mid n_{1+}, n_{+1}, n, \theta) - \dfrac{\alpha}{2} &= \dfrac{\displaystyle\sum\limits_{t \geq n_{11}}\binom{n_{1+}}{t} \binom{n-n_{1+}}{n_{+1} - t} \theta^t}{\displaystyle\sum_{u=m_{-}}^{m_{+}}\binom{n_{1+}}{u} \binom{n-n_{1+}}{n_{+1} - u} \theta^u} - \dfrac{\alpha}{2} \label{upper-poly}\\ 
\sum\limits_{t \leq n_{11}}f(t \mid n_{1+}, n_{+1}, n, \theta) - \dfrac{\alpha}{2} &= \dfrac{\displaystyle\sum\limits_{t \leq n_{11}}\binom{n_{1+}}{t} \binom{n-n_{1+}}{n_{+1} - t} \theta^t}{\displaystyle\sum_{u=m_{-}}^{m_{+}}\binom{n_{1+}}{u} \binom{n-n_{1+}}{n_{+1} - u} \theta^u} - \dfrac{\alpha}{2}\text{.}
\end{align}``
These two problems can be described more succinctly as follows: let ``$\{a_t\}$`` and ``$\{b_u\}$`` be non-negative finite (with starting and stopping points) sequences. Then we wish to find the root of
``\begin{equation}
\dfrac{\displaystyle\sum_{t} a_t \theta^t}{\displaystyle\sum_{u}b_u\theta^u} - \dfrac{\alpha}{2}\text{.} 
\end{equation}``
The above function of ``$\theta$`` is problematic. First of all, it is not continuous at ``$\theta = 0$``, and who knows where else it could be discontinuous? Note the denominator of the first fraction implies that the function above is discontinuous for all ``$\theta$`` satisfying ``$\sum_{u}b_u\theta^u = 0$``. Therefore, using Newton-Raphson on this function of ``$\theta$`` above is quite risky. Furthermore, since we do not know anything about the powers of ``$\theta$`` beforehand, we can't determine where the function of ``$\theta$`` above will be positive or negative, so the bisection method is off limits.

Instead, let's choose to ignore the problem of continuity. Then what happens is as follows: set the aforementioned function of ``$\theta$`` equal to ``$0$`` to obtain
``\begin{equation*}
\dfrac{\displaystyle\sum_{t} a_t \theta^t}{\displaystyle\sum_{u}b_u\theta^u} - \dfrac{\alpha}{2} = 0 \implies \sum_{t} a_t \theta^t = \dfrac{\alpha}{2}\sum_{u}b_u\theta^u \implies \sum_{t} a_t \theta^t - \dfrac{\alpha}{2}\sum_{u}b_u\theta^u = 0\text{.}
\end{equation*}``
We define
``\begin{equation}
\text{poly}(\theta) = \sum_{t} a_t \theta^t - \dfrac{\alpha}{2}\sum_{u}b_u\theta^u\text{.}
\end{equation}``
This function is a difference of two polynomials, so that it is not only continuous in ``$\mathbb{R}$``, but differentiable in ``$\mathbb{R}$`` as well. ``$\text{poly}$`` also has the same roots as our previous problematic function, as well as a few additional roots. 

But this is a chance we are willing to take. Why? First of all, note that in the original equations  

``\begin{align}
\sum\limits_{t \geq n_{11}}f(t \mid n_{1+}, n_{+1}, n, \theta) - \dfrac{\alpha}{2} &= \dfrac{\displaystyle\sum\limits_{t \geq n_{11}}\binom{n_{1+}}{t} \binom{n-n_{1+}}{n_{+1} - t} \theta^t}{\displaystyle\sum_{u=m_{-}}^{m_{+}}\binom{n_{1+}}{u} \binom{n-n_{1+}}{n_{+1} - u} \theta^u} - \dfrac{\alpha}{2} \\ 
\sum\limits_{t \leq n_{11}}f(t \mid n_{1+}, n_{+1}, n, \theta) - \dfrac{\alpha}{2} &= \dfrac{\displaystyle\sum\limits_{t \leq n_{11}}\binom{n_{1+}}{t} \binom{n-n_{1+}}{n_{+1} - t} \theta^t}{\displaystyle\sum_{u=m_{-}}^{m_{+}}\binom{n_{1+}}{u} \binom{n-n_{1+}}{n_{+1} - u} \theta^u} - \dfrac{\alpha}{2}\text{.}
\end{align}``

by the fact that ``$m_{-} \leq t \leq m_{+}$`` that - if we ignore the ``$\alpha/2$`` term for the moment - the remaining fraction has a denominator which must contain each term of the numerator. Therefore, in our more general formulation of the problem, it follows that ``$\sum_{u}b_u\theta^u$`` must contain each term of ``$\sum_{t} a_t \theta^t$``. 

Looking back at ``$\text{poly}(\theta)$``, what does this imply? Since ``$\{a_t\}$`` and ``$\{b_u\}$`` are non-negative, this implies that where ``$t = u$``, the coefficient of ``$\theta^t$`` must be ``$a_t\left(1 - \dfrac{\alpha}{2}\right)$`` in ``$\text{poly}(\theta)$``, which is a non-negative coefficient. Every other case where there is a coefficient of ``$\theta^u$`` which is not in the ``$\theta^t$`` summation yields a negative coefficient.

In other words, we may express ``$\text{poly}(\theta)$`` as follows:
``\begin{equation}
\text{poly}(\theta) = \sum_{t} a_t \left(1-\dfrac{\alpha}{2}\right) \theta^t - \dfrac{\alpha}{2}\sum_{v}c_v\theta^v
\end{equation}``
where ``$\{c_v\}$`` is a non-negative finite sequence. We also know that the ``$t$``-index set is either the set such that ``$t \geq n_{11}$`` or ``$t \leq n_{11}$``. Therefore, the ``$v$``-index set must be the complement of the ``$t$``-index set, and we have our main result:

**Theorem**. There is only one pair of sign changes in the coefficients of     ``$\text{poly}(\theta)$`` when the terms are arranged in ascending power of ``$\theta$``.

By Descartes' Rule of Signs,^[Weisstein, Eric W. "Descartes' Sign Rule." *Mathworld* \[online\]. Available at http://mathworld.wolfram.com/DescartesSignRule.html.] the previous theorem immediately implies that

**Corollary**. ``$\text{poly}(\theta)$`` has at most one positive root.

Therefore, we are safe in using Newton-Raphson on ``$\text{poly}$`` to find its positive roots as ``$\text{poly}$`` has at most one positive root (the odds ratio we desire), and as long as we choose a reasonable guess to begin with. The derivative of ``$\text{poly}$`` is
``\begin{equation}
\text{poly_deriv}(\theta) = \sum_{t} ta_t \theta^{t-1} - \dfrac{\alpha}{2}\sum_{u}ub_u\theta^{u-1}\text{.}\label{poly-deriv}
\end{equation}``
Therefore, the Newton-Raphson procedure for finding the root of ``$\text{poly}$`` is as follows:

1. Choose an initial ``$\theta^{(0)}$``.

2. Until convergence is met, or the number of iterations exceeds the maximum number of iterations desired, set for ``$k = 1, 2, \dots$``:
``\begin{equation*}
\theta^{(k)} \leftarrow \theta^{(k-1)} - \dfrac{\text{poly}\left(\theta^{(k-1)}\right)}{\text{poly_deriv}\left(\theta^{(k-1)}\right)}\text{.}
\end{equation*}``

This computation is performed by using two of three arrays as inputs:

- the array of coefficients and powers of the "upper polynomial"     ``$\displaystyle\sum\limits_{t \geq n_{11}}\binom{n_{1+}}{t} \binom{n-n_{1+}}{n_{+1} - t} \theta^t$``,

- the array of coefficients and powers of the "lower polynomial" ``$\displaystyle\sum\limits_{t \leq n_{11}}\binom{n_{1+}}{t} \binom{n-n_{1+}}{n_{+1} - t} \theta^t$``, and 

- the array of coefficients and powers of the "denominator polynomial" ``$\displaystyle\sum_{u=m_{-}}^{m_{+}}\binom{n_{1+}}{u} \binom{n-n_{1+}}{n_{+1} - u} \theta^u$``.

Two of these three arrays (either the "upper" and "denominator" **or** the "lower" and "denominator") are then used as inputs to determine ``$\text{poly}$`` and ``$\text{poly_deriv()}$``.