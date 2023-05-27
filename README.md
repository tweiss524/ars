# Adaptive Rejection Sampling Algorithm

This project was completed in collaboration with Jeffrey Kuo and Hsiang-Chuan Sha at the University of California, Berkeley.

The following paragraphs give a brief run-through of the contents of this repository. A more detailed version of the steps is written in [ars_report.pdf](report/ars_report.pdf).

# Package Installation

In order to install the package in RStudio, refer to the following instructions. Note that required packages are `assertthat`, `numDeriv`, `stats`, `rmutil`, and `testthat`. Also remember to clean your environment before installing the package to ensure there are no errors with installation.
1. clone the repository to your local environment using git clone "https://github.berkeley.edu/tweiss/ars.git"
2. set your working directory in RStudio to the `ars` directory you just cloned
3. type `devtools::load_all()` into your console
4. `library(ars)`

# Introduction

This project aims to create an `R` package called ars that simulates the adaptive rejection sampling algorithm described in the paper Adaptive Rejection Sampling for Gibbs Sampling by W.R. Gilks and P Wild. In this paper, Gilks and Wild propose a method for rejection sampling which is valid for any univariate, log-concave probability density function. As part of this project, we created a main file called `ars.R` which contained the main ars() function used to calculate and produce the samples, with helper functions to complete each step in the algorithm defined in `ars_functions.R`.

# Methods

The main function to calculate and produce the samples is called `ars()`. For purposes of organization and readability, we decided to use a functional programming approach where we defined a number of auxiliary functions to complete subtasks of the `ars()` function, and wrote all of these functions in a separate `R` script entitled `ars_functions.R`. The main `ars()` function takes in the following inputs from the user:
- `f`: The density function from which the user wants to sample from. Requirement of the argument is that it must be a univariate, log-concave function.
- `n`: The number of samples the user wants to generate from the function `f`. Requirement of the argument is that it must be a positive integer value.
- `bounds`: denotes the left and right bounds of the domain of the function `f`, which is the set of points $x$
for which $f(x) > 0$. Requirement of the argument is that it be a vector of length 2, where the first element is the left
(lower) bound and the second element is the right (upper) bound.
- `x_init`: The initial point within the bounds used as one of the initial abscissae points for future
sampling. Requirement of the argument is that it is a single point defined within the bounds.

The first task the function carries out is checking the validity of the user inputs. After passing these initial argument checks, the `ars()` function moves on to defining two simple functions: `h()` and `hprime()`, where both functions take in a point $x$. `h()` returns the value $\log f(x)$ and `hprime()` returns the derivative of this value; that is, $\frac{d}{dx}\log f(x)$, where we used the external package numDeriv to compute the derivative. Next, the `ars()` function is broken down into the three steps Gilks and Wild suggest: (1) initialization, (2) sampling, and (3) updating.

# Initialization Step

The basic idea of the initialization step is to initialize the abscissae points, which we define in the function as $T_k$. In order to carry out this step, we define an auxiliary function `initialize_abscissae()`, which uses the `x_init` point as a reference in order to output a vector of $k = 20$ initial abscissae points. After defining these 20 initial abscissae points and computing their derivatives, we initialize the functions $u_k(x)$, $s_k(x)$, $l_k(x)$ and calculate the intersection of the tangent lines, $z_j$, for each $x_j$ and $x_{j+1}$, where $j = 1,\dots, k − 1$. As described in the paper by Gilks and Wild, we derive the rejection envelope on $T_k$ and calculate the value of every $x_k$ on the rejection envelope, which was done utilizing the function $u_k(x)$. This function first identifies the segment where the inputted value $x$ lies, and calculates the $u_k(x)$ for the corresponding slope and intercept of the piecewise function. We similarly calculate the squeezing function $l_k(x)$ with our function $l_k(x)$. Therefore, we derived both the rejection sampling value and the squeezing value for each point within $T_k$ from $u_k(k)$, and $l_k(x)$, respectively. For $s_k(x)$, we approximate the density function by calculating the area under the piecewise $e^{u_j(x)}$ for each segment $(z_j, z_{j+1})$, including the bounds the function is defined on as $z_0$ and $z_k$, in order to find the unnormalized probability $q_j$ , which takes the following form:

$$ q_j=\int_{z_j}^{z_j+1}e^{u_j(x)}dx$$

As defined in the paper, $u_j(x) = h(x_j) + (x-x_j)h'(x_j)$, which simplifies the integral to two cases based on the value of $h'(x_j)$:

$$ 
q_j=
 \begin{cases}
   e^{u_j(z_j)}(z_{j+1}-z_j) & \text{if } h'(x_j)=0 \\
   \frac{e^{u_j(z_{j+1})}-e^{u_j(z_j)}}{h'(x_j)} & \text{if } h'(x_j)\neq 0 
  \end{cases}
$$

The normalized probability of sampling each segment is derived by dividing each $q_j$ with the summation of all of the unnormalized probabilities $\sum\limits_{j=1}^{k+1}q_j$. The probability to sample each interval $(z_j,z_{j+1})$ is therefore equal to $$p_j = \frac{q_j}{\sum\limits_{j=1}^{k+1}q_j}$$ 

# Sampling Step

We followed the sampling method that Gilks and Wild derived. That is, we first sample the segment using the probabilities calculated above. Then we sample an $x^\*$ with our `sample_sk()` function by applying the inverse CDF of that segment to a sampled draw of $w$ from a $Unif(0, 1)$ distribution which represents the probability that $X$ is less than or equal to $x^\*$:

$$w = \mathbb{P}(X\leq x^\*) = \frac{\int^{x^\*}_{z_j}e^{u_j(x)}dx}{\int^{z_{j+1}}_{z_j}e^{u_j(x)}dx}$$

where the denominator is the area under the curve of the sampled interval calculated above. Let $\int^{z_{j+1}}_{z_j}e^{u_j(x)}dx = A_j$ for simplicity. For $h'(x_j)\neq0$, we solve for $x^\*$ in the equation:
$$
\begin{aligned}
w &= \frac{\int^{x^\*}_{z_j}e^{u_j(x)}dx}{A_j} \\
\implies x^\* &= \begin{cases}
                  \frac{\log\left(wAh'(x_j)+e^{u_j(z_j)}-h(x_j)\right)}{h'(x_j)} +x_j & \text{if } h'(x_j)\neq 0\\
                  z_j+w(z_{j+1}-z_j) & \text{if } h'(x_j)=0 
                \end{cases}
\end{aligned}
$$
where the $h'(x_j) = 0$ case is akin to scaling the uniform by the length of the interval.\
\
We then perform the following tests:\
\
First is the squeezing test: if $w\leq e^{l_k(x^\*) - u_k(x^\*)}$, we add the value $x^\*$ to our sample vector. If this is not the case, we move on to the rejection test: if $w\leq e^{h(x^\*)-u_k(x^\*)}$, we add the value $x^\*$ to our sample vector. Otherwise, $x^\*$ is rejected. In the event a value of $x^\*$ leads us to perform the rejection test, we also complete the updating step, defined in the next section.


# Updating Step

As just mentioned, if $x^\*$ fails the squeezing test, we perform the updating step, which works as follows: we append the value of $x^\*$ to our vector $T_k$ such that we now have $T_{k+1}$. We then sort the values of $T_{k+1}$ so they are in ascending order, and compute $h'(x^\*)$, placing it in its analogous location to its position in $T_{k+1}$. We then calculate the new intersection points of the tangent lines in of the points in $T_{k+1}$ and return to the beginning of the sampling step. We complete this process until we have obtained a total of $n$ samples, and the vector of $n$ samples are returned from the `ars()` function to the user.

# Tests
In order to run our specified set of tests, refer to the following instructions.
1. In order to test the main ars() function: testthat::test_file("tests/testthat/test-ars.R")
2. In order to test the auxiliary functions: testthat::test_file("tests/testthat/test-ars_functions.R")

In order to run the tests for our function, we created an R script which utilized the R package testthat. We designed a number of tests in order to ensure that our function works as desired, including checks for log-concavity of the function f, checks for invalid function inputs, and checks that the function returned a vector of samples that resembled a vector of samples of the same length from its respective distribution. In order to do this, we implemented the Kolmogorov-Smirnov test.

In order to run this test in R, we used the function ks.test() and ensured that the obtained p-value was greater than $\alpha = 0.05$. We also wrote an R script to test that our auxiliary functions were providing reasonable results; namely `check_log_convave()`, `initialize_abscissae()`, `calc_z()`, `u()`, `l()`, `calc_probs()`, and `sample_sk()`. The tests for `check_log_convave()` aimed at targeting vectors with non-decreasing values, and we also included a test case for a vector where every value was the same (as this would still be log-concave). The tests for `initialize_abscissae()` aimed to try each case where the left bound was -Inf, the right bound was Inf, or both bounds were finite, ensuring that a vector of length $k=20$ was properly returned for each test case. For testing `calc_z()`, `u()`, `l()`, `calc_probs()`, and `sample_sk()`, we provided sample values for the functions to perform computations on, and then compared the result of the function with a manual computation carrying out the same task for a variety of test cases. Given that each of these functions work under a "black-box" methodology, the sample values provided were random and contained no real meaning, but were still effective in testing the basic tasks of each function.

To view examples of the final product, refer to the last section of [ars_report.pdf](report/ars_report.pdf).

