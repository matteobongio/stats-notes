#import "@preview/mannot:0.2.3": *
#import "@preview/codly:1.3.0": *
#import "@preview/codly-languages:0.1.1": *
#show: codly-init.with()
#codly(languages: codly-languages)


#let title = [Stats]
#set page(
   paper: "a4",
   header: align(left, title),
   numbering: "1",
)

#align(center, text(17pt)[
   *#title*
])

$ 
P(A|B) = frac(P(A inter B), p(B)) \
"assume" P(B) != 0 \
P(epsilon_1 union epsilon_2 | H) = P(epsilon_1 | H) + P(epsilon_2 | H) 
- P(epsilon_1 inter epsilon_2 | H)
$

= Law of total Probablilities
$ P(epsilon) = P(epsilon | H) P(H) + P(epsilon | H^c) P(H^c) $

= Bayes' Formula
$ P(H|E) = 
frac(
P(epsilon | H) P(H),
P(epsilon)
) \
= frac(
P(epsilon | H) P(H),
P(epsilon|H)P(H) + P(epsilon|H^c)P(H^c)
) $

= Discrete Random Variable

*Probablility Mass Function*

$ f_x (x) = P(X = x) "forall" x $

= Continuous Random Variable

*Probablility Mass Function*

A random vvariable $X$ is continuous if its cumulative distrobution
function $F_x(x)$ is continuous

$
P(X <= x) = F_x (x)
= integral^x_(- infinity) f_x (y) d y "forall" x in RR \
"for all real numbers" a < b \
P(a <= X <= b) = integral^a_b f_x (y) d y
$

= Expectation
$ 
E[x] = cases(
Sigma_x x f_x (x) "(discrete)",
integral^infinity_(-infinity) x f_x (x) d x "(pdf "f_x")"
) 
$

= Variation
$
"var"(x) = E[x^2] - (E[x])^2
$

= Summary Methods
$
"mean"(x) = frac(1, N) Sigma^N_(i=1) x_i \
"median"(x) = "middle value" \
"Mode of " x: "most frequent value" \
"quantile"(x, k) = "the value such that " k% " of the data is " <= \
"IQR" = "quantile"(x, 75) - "quantile"(x,25)
$

= Discrete Random Variation

== Uniform distrobution

$
D(X=x) = frac(1, N) \
E[x] = frac(1, N) Sigma^N_(i = 1) x_i = mu \
"Var"(X) = frac(1, N) Sigma^N_(i=1) (x_i - mu)^2
$

== Bernoulli Distrobution

Used to model events with 2 mutually exclusive outcomes

$
P(X = 1) = p \
E[X] = p \
"Var"(X) = p (1-p)
$

== Binomial Distrobution
finite number of independent Bernoulli trials

$
P(X = i) = binom(N, i) p^i (1-p)^(N-i) \
E[X] = N p \
"Var"(x) = N p (1-p)
$

== Geometric Distrobution
Is used to model an infinite suite of independent events with the
same probability of success. $P(X = i)$ - the probability that the
event is a success for the 1st time at the i-th trial.

$
P(X = i) = (1 - p)^(i - 1) p \
E[x] = frac(1, p) \
"Var"(X) = frac(1-p,p^2)
$

== Poisson Distrobution

Used to model the number of occurences of an event in a time interval or space
that occur independent from the last time it occured. $P(X=i)$ represents
the probablility that the even $X$ takes place $i$ times per unit (space/time)

$
P(X=i) = frac(e^(-lambda) lambda^i, i!) \
E[X] = lambda
"Var"(X) = lambda
$

*Poisson Convergence*: Binomial distrobution convergences to Poisson

$
P(N p) ~ B (N, p) "(for large "N" and small "N p")"
$

*DeMoivre - Laplace*
$ B(N, p) ~ N(N p, sqrt(N p q)) "For large "N \ q = 1-p $


= Distrobutions in R

#align(center)[
#show table.cell: it => {
   set text(18pt)
      set align(left)
   if it.x == 0 or it.y == 0 {
      set text(orange)
      strong(it)
   } else {
      set text(blue)
         it
   }
}
#table(
stroke: none,
gutter: 0.2em,
columns: 4,
table.header(
   [Distrobution], [Sample], [PMF], [CMF]
   ),
   [ Uniform ] , [ runif ], [ dunif ], [ punif ],
   [ Binomial ] , [ rbinom ], [ dbinom ], [ pbinom ],
   [ Geometric ] , [ rgeom ], [ dgeom ], [ pgeom ],
   [ Poisson ] , [ rpois ], [ dpois ], [ ppois ],
)
]

$
("PMF") P(Y=x) \
("CMF") P(Y<=x) \
(1-"CMF") P(Y>x)
$

= Uniform Random Variables (Uniform Distribution)

We don't know anything except they have upper and lower bound

$
f(x) = cases(frac(1, b-a)", " a <= x < b, 0 )\
E[x] = frac(a + b, 2) \
"Var"(X) = frac((b-a)^2, 12)
$

= Bayesian Inference

$
mark(P(H|E), tag: #<post>, color: #blue) = mark(P(H), tag:#<prior>) mark(frac(P(E|H),P(E)), tag: #<marg>, color: #green)

#annot(<post>, pos:bottom+left, yshift: 2em)[
#set text(8pt)
Posterior Probablility
]
#annot(<prior>, pos:bottom, yshift: 2em)[
#set text(8pt)
Prior Probablility
]
#annot(<marg>, pos:bottom+right, yshift: 2em)[
#set text(8pt)
Marginal liklihood of $epsilon$ given $H$
]
$
#v(10em)


#align(center)[
#show table.cell: it => {
   set text(18pt)
      set align(left)
   if it.y == 0 {
      set text(orange)
      strong(it)
   } else {
      set text(blue)
         it
   }
}
#table(
stroke: none,
gutter: 0.2em,
columns: 3,
table.header(
   [liklihood], [Prior], [Posterior] ),
   [Binom], [Beta], [$alpha'=alpha+k, beta'=beta+n-k$],
   [Poisson], [Gamma], [$alpha'=alpha+Sigma^n_(i=1) x_i, beta'=beta + n$],
   [Normal], [Normal], 
   [$mu'=(
   frac(mu_0, sigma^2_0)
   + 
   frac(Sigma^n_(i=1) x_i, sigma^2)
   )
   (
   frac(1, sigma^2_0)
   +
   frac(n, sigma^2)
   )^(-1)
   \
   sigma^2'=(frac(1, sigma^2_0) + frac(n, sigma^2))^(-1)
   $],
)
]

= Beta

$
f (x) = frac(Gamma (alpha + beta), Gamma(alpha) Gamma(beta)) x^(alpha - 1) (1-x)^(beta - 1) \
E[X] = frac(alpha, alpha + beta) \
"Var"(X) = frac(alpha beta, (alpha + beta)^2 (alpha + beta + 1))
$

= Gamma

$
f (x) = frac(beta^alpha, Gamma(alpha)) x^(alpha - 1) e^(- beta x) \
E[X] = frac(alpha, beta) \
"Var"(X) = frac(alpha, beta^2)
$

= Exponential Distrobution
time or space between 2 events

$
f (x) = cases(lambda e^(-lambda x)  x >= 0, 0) \
E[X] = frac(1, lambda) \
"Var"(X) = frac(1, lambda^2)
$

= Normal Distrobution

$
f (x) = frac(1, sigma sqrt(2 pi)) e^(- frac((x - mu)^2, 2 sigma^2)) \
E[X] = mu \
"Var"(X) = mu^2 \ 
N(0, 1) "(is the standard normal distrobution)"
$

= Central Limit Theorem

let $X_1, ..., X_N$ be a sequencee of independent and identically distributed random
variables, each with the same mean and variance

$
frac(
Sigma^N_(i=1) X_i - N_mu,
sqrt(N sigma^2)
)
~ N(0,1)
$

= Chebyshev's Inequality

$X$ be a random variable

$
forall t > 0, t = b - frac(a + b, 2) \
P(|X - E[X]| >= t) <= frac("Var"(x), t^2)
$

= Markov's Inequality

How fast $P(X >= t)$ goes to $0$ as $t -> infinity$

$X$ be a random non-negative with variance


$
forall t > 0, P(X >= t) <= frac(E[X], t)
$

= Hypothesis Testing

$
H_0 : "null hypothesis, assumed true" \
H_a : "alternative hypothesis " not H_0
$

*Significance level*

$ alpha : "probability of rejecting " H_0 " when it is true" $

*p-value*: lowest level $alpha$ to which $H_0$ could have been rejected

*Outcomes*: 
- reject $H_0 => H_a$
- do not reject $H_0 =>$ don't claim anything

$ p < a => "reject " H_0", else nothing" $

*type 1 error*: rejecting $H_0$ when its true

*type 2 error*: failing to reject when false

+ Formulate $H_0$ and $H_a$
+ specify significance level $alpha$ usually ${0.1, 0.05, 0.01}$
+ Given $H_0$ look up critical value of p-value
+ Do not reject $H_0$ if found value $>$ significance level


= Test Formulate

*Mean of Normal Population with _known_ variance*

$ z = frac( overline(x) - mu_0, sigma / sqrt(n)) ~ N(0,1) $

*Mean of Normal Population with _unknown_ variance*

$ T = frac(overline(x) - mu_0, S / sqrt(n)) ~ t_(n - 1) $

*Proportion of a Population*

$ z = frac(accent(p, hat) - p_0, sqrt(p_0 (1 - p_0)/n)) ~ N(0,1) $

*Difference of means of Normal Population with _unknown_ variance*

$
T = frac(
(overline(x_1) - overline(x_2)) - (mu_1 - mu_2),
S_p (1/n_1 + 1/n_2)
) ~ t_(n_1 + n_2 - 2)
$

#let hp = $accent(p, hat)$
*Difference of 2 proportions*
$
z = frac( hp_1 - hp_2, sqrt( hp (1 - hp) ( 1/n_1 + 1/n_2))) ! N(0,1) \
#hp = frac(n_1 hp_1 + n_2 hp_2, n_1 + n_2)
$

*$chi^2$ Test for goodness of fit*

$ chi^2 = Sigma^k_(i=1) frac((O_i - E_i)^2, E_i) ~ chi^2_(k-1) $

*$chi^2$ Test for independence*

$
chi^2 = Sigma^k_(i = 1) Sigma^r_(j = 1) frac((O_(i j) - E_(i j))^2, E_(i j)) ~ chi^2_((k - 1)(r - 1))
$

= R Code
== Libraries
```R
library(tidyverse)
# cool stuff
data({dataset})
summarise({dataset})
mean({dataset})
median({dataset})
sd({dataset})
```

== filtering
```R
data$Model <- rownames(data)
data %>% filter(str_detect(Model,'Merc')) %>% summarise(mean(mpg))
```

== GGplot

*General Structure*
```R
ggplot(data = <DATA>, aes(x = <X>, y = <Y>)) + <GEOM_FUNCTION>()
```

- `geom_point()` – Scatter plot
- `geom_line()` – Line chart
- `geom_bar()` – Bar chart
- `geom_histogram()` – Histogram
- `geom_boxplot()` – Boxplot

== Distrobution Simulation
```R
#define the outcome space
outcomeSpace <- c(1,2,3,4,5,6)
#define the function to simulate the rolls
getOne <- function(){
   trial <- 0
   while (TRUE){
      outcome <- sample(outcomeSpace,1,prob = rep(0.166666667,6))
      trial <- trial +1
      if (outcome==1){
         break
      }
   }
   return(trial)
}
#replicate the experiment
avgNumberOfTrials <- c()
for (i in c(10,50,100,500,1000,5000,10000,50000,100000)) {
   set.seed(999)
   trials <- replicate(i,getOne())
   avgNumberOfTrials <- append(avgNumberOfTrials,mean(trials))
}
```

#pagebreak()

== Bayesian Inference: Calculating the Posterior Distrobution

=== Normal Distrobution

```R
data(ChickWeight)
weights <- ChickWeight$weight[ChickWeight$Diet == "1"]
prior_mean <- 200
prior_var <- 20ˆ2
# Data summary
n <- length(weights)
sample_mean <- mean(weights)
sample_var <- var(weights)
# Posterior parameters
posterior_sd <- sqrt(1 / (1 / prior_var + n / sample_var))
posterior_mean <- posterior_sd**2 * (prior_mean / prior_var + n * sample_mean / sample_var)
posterior_mean <- round(posterior_mean, 2)
glue("Posterior_Mean = {posterior_mean}")
```

== Credible Interval
```R
z_alpha2 <- -qnorm(0.05/2)
lowerbound <- posterior_mean - z_alpha2 * posterior_sd
upperbound <- posterior_mean + z_alpha2 * posterior_sd
c(lowerbound, upperbound)
```

== Confidence Interval
```R
# two sided
supp_oj <- ToothGrowth$len[ToothGrowth$supp == 'OJ']
supp_vc <- ToothGrowth$len[ToothGrowth$supp == 'VC']
t_test_result <- t.test(supp_oj, supp_vc)

```

== Linear Regression Model
```R
data <- lm(iris$Sepal.Length - iris$Petal.Width)
```
