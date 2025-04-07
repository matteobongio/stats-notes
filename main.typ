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
P(Epsilon_1 union Epsilon_2 | H) = P(Epsilon_1 | H) + P(Epsilon_2 | H) 
- P(Epsilon_1 inter Epsilon_2 | H)
$

= Law of total Probablilities
$ P(Epsilon) = P(Epsilon | H) P(H) + P(Epsilon | H^c) P(H^c) $

= Bayes' Formula
$ P(H|E) = 
frac(
P(Epsilon | H) P(H),
P(Epsilon)
) \
= frac(
P(Epsilon | H) P(H),
P(Epsilon|H)P(H) + P(Epsilon|H^c)P(H^c)
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

$
P(X = i) = (1 - p)^(i - 1) p \
E[x] = frac(1, p) \
"Var"(X) = frac(1-p,p^2)
$

== Poisson Distrobution
