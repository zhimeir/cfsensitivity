# cfsensitivity
An R package for the conformalized counterfactual estimation with the presence of unobserved confounders.
The methodology is based on the paper: [Sensitivity Analysis of Individual Treatment Effects: A Robust Conformal Inference Approach](https://arxiv.org/abs/2111.12161).

## Overview
This R package implements the method for sensitivity analysis of individual treatment effects 
proposed in our paper: Sensitivity Analysis of Individual Treatment Effects: A Robust Conformal 
Inference Approach.

- `cfsens_mgn()` conducts sensitivity analysis of individual treatment effects with the marginal guarantee.
- `cfsens_pac()` conducts sensitivity analysis of individual treatment effects with the PAC-type guarantee.
- `cfsens_cf_mgn()` produces robust predictice intervals for the counterfactuals with the marginal guarantee.
- `cfsens_cf_pac()` produces robust predictice intervals for the counterfactuals with the PAC-type guarantee.


## Installation
To install the package, run the following commands in R:
```{r}
if (!require("devtools")){
    install.packages("devtools")
}
devtools::install_github("zhimeir/cfsensitivity")
```
The [grf](https://github.com/grf-labs/grf) package is required.

## Usage examples
We illustrate the usage of the package with synthetic dataset.
```{r}
## Data generating model
data_model <- function(n, p, Gamma, seed){
  set.seed(seed)
  beta <- matrix(c(-0.531,0.126,-0.312,0.018,rep(0,p-4)), nrow=p)
  X <- matrix(runif(n * p), n, p)
  U <- rnorm(n) * abs(1 + 0.5 * sin(2.5 * X[,1]))
  Y1 <-  X %*% beta + U
  Y0 <- X %*% beta - U
  prop.x <- exp(X %*% beta) / (1 + exp(X %*% beta))
  p.x <- 1-(1/(prop.x + (1-prop.x)/Gamma ) -1)/( 1/(prop.x + (1-prop.x)/Gamma) - 1/(prop.x + Gamma*(1-prop.x)))
  t.x <- qnorm(1-p.x/2) * abs(1+0.5*sin(2.5*X[,1]))
  prop.xu <- (prop.x/(prop.x+Gamma*(1-prop.x)))*(abs(U)<=t.x) + (prop.x/(prop.x+ (1-prop.x)/Gamma))*(abs(U)>t.x)
  TT <- rbinom(n, size=1, prob=prop.xu)
  Y <- TT * Y1 + (1 - TT) * Y0
  
  return(list(Y1 = Y1, Y0 = Y0, Y = Y, X = X, T = TT))
}

## Generate confounded training data
n <- 5000
n_test <- 5000
p <- 10
Gamma <- 1.5
data <- data_model(n, p, Gamma, 2021)

## Generate confounded test data
test <- data_model(n_test, p, Gamma, 2022)

## Run sensitivity analysis with marginal guarantee
## grf pacakge needs to be installed
alpha <- 0.2
res <- cfsens_mgn(data$X, data$Y, data$T, 
                  alpha = alpha, null_type = "negative", 
                  ps_fun = regression_forest, 
                  pred_fun = quantile_forest)

out_res <- predict(res, test$X, Y1_test = test$Y1, type = "att")
out_res

## Run sensitivity analysis with PAC-type guarantee
## grf package needs to be installed
alpha <- 0.2
delta <- 0.1
res <- cfsens_pac(data$X, data$Y, data$T, 
                  alpha = alpha, delta = delta, null_type = "negative",
                  ps_fun = regression_forest, 
                  pred_fun = quantile_forest)

out_res <- predict(res, test$X, Y1_test = test$Y1, type = "att")
out_res




## Construct predictive intervals with marginal guarantees
## grf package needs to be installed
alpha <- 0.2
res <- cfsens_cf_mgn(data$X, data$Y, data$T, 
                     Gamma = Gamma, alpha = alpha,
                     ps_fun = regression_forest, 
                     pred_fun = quantile_forest)

out_res <- predict(res, test$X, estimand = "treated", type = "ate")
cover <- mean((test$Y1 >= out_res$Y_lo) * (test$Y1 <= out_res$Y_hi))
cover

## Construct predictive intervals with marginal guarantees
## grf package needs to be installed
alpha <- 0.2
delta <- 0.1
res <- cfsens_cf_pac(data$X, data$Y, data$T, 
                     Gamma = Gamma, alpha = alpha,
                     delta = delta,  
                     ps_fun = regression_forest, 
                     pred_fun = quantile_forest)

out_res <- predict(res, test$X, estimand = "treated", type = "ate")
cover <- mean((test$Y1 >= out_res$Y_lo) * (test$Y1 <= out_res$Y_hi))
cover
```
