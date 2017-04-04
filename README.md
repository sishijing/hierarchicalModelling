# hierarchicalModelling
several ways to perform hierarchical modelling

Here I will provide several algorithm to fit hierarchical
models

First level: Y_{i}\sim{N}(mu_{i}, sigma^{2});
Second level: mu_{i}\sim{N}(gamma, tau^{2});
with Data=(Y_{1},..., Y_{n}), Mu=(mu_{1}, ..., mu_{n}),
sigma the common standard deviation for the observed data.
Mu is the true values, in real cases unknown.

We want to use Bayesian hierarchical models to obtain samples for parameters mu_{1}, ..., mu_{n}
