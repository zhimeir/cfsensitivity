#' Robust conformal counterfactual inference with PAC-type guarantee
#'
#' \code{cfsens_cf_pac} constructs robust predictive intervals for 
#' counterfactuals with the PAC guarantee.
#'
#'
#' @details When \code{side = "two"}, the predictive interval takes the form [a,b];
#' when \code{side = "above"}, the predictive interval takes the form [a,Inf);
#' when \code{side = "below"}, the predictive interval takes the form (-Inf,a].
#'
#' When \code{null_type = "sharp"}, the null hypothesis is H_0: Y(1) - Y(0) = 0.
#' When \code{null_type = "negative"}, the null hypothesis is H_0: Y(1) - Y(0) <= 0.
#' When \code{null_type = "positive"}, the null hypothesis is H_0: Y(1) - Y(0) >= 0.
#'
#'
#' @param X covariates.
#' @param Y the observed outcome vector.
#' @param T the vector of treatment assignments.
#' @param Gamma The confounding level.
#' @param alpha the target confidence level.
#' @param delta the target confidence level over the randomness of calibration set.
#' @param side the type of predictive intervals that takes value in \{"two", "above", "below"\}. See details.
#' @param score_type the type of nonconformity scores. The default is "cqr".
#' @param ps_fun a function that models the treatment assignment mechanism. The default is "regression_forest".
#' @param ps a vector of propensity score. The default is \code{NULL}. 
#' @param pred_fun a function that models the potential outcome conditional on the covariates. The default is "quantile_forest".
#' @param train_prop proportion of units used for training. The default is 75\%. 
#' @param train_id The index of the units used for training. The default is \code{NULL}.
#'
#' @return an \code{cfpac} object.
#'
#' @seealso 
#' \code{\link{cfsens_cf_mgn}}
#'
#' @examples
#' \donttest{
#'
#' ## Data generating model
#' data_model <- function(n, p, Gamma, seed){
#'  set.seed(seed)
#'  beta <- matrix(c(-0.531,0.126,-0.312,0.018,rep(0,p-4)), nrow=p)
#'  X <- matrix(runif(n * p), n, p)
#'  U <- rnorm(n) * abs(1 + 0.5 * sin(2.5 * X[,1]))
#'  Y1 <-  X %*% beta + U
#'  Y0 <- X %*% beta - U
#'  prop.x <- exp(X %*% beta) / (1 + exp(X %*% beta))
#'  p.x <- 1-(1/(prop.x + (1-prop.x)/Gamma ) -1)/( 1/(prop.x + (1-prop.x)/Gamma) - 1/(prop.x + Gamma*(1-prop.x)))
#'  t.x <- qnorm(1-p.x/2) * abs(1+0.5*sin(2.5*X[,1]))
#'  prop.xu <- (prop.x/(prop.x+Gamma*(1-prop.x)))*(abs(U)<=t.x) + (prop.x/(prop.x+ (1-prop.x)/Gamma))*(abs(U)>t.x)
#'  TT <- rbinom(n, size=1, prob=prop.xu)
#'  Y <- TT * Y1 + (1 - TT) * Y0
#'  
#'  return(list(Y1 = Y1, Y0 = Y0, Y = Y, X = X, T = TT))
#'}
#'
#' ## Generate confounded training data
#' n <- 5000
#' n_test <- 5000
#' p <- 10
#' Gamma <- 1.5
#' data <- data_model(n, p, Gamma, 2021)
#'
#' ## Generate confounded test data
#' test <- data_model(n_test, p, Gamma, 2022)
#'
#' ## Construct predictive intervals with marginal guarantees
#' ## grf package needs to be installed
#' alpha <- 0.2
#' delta <- 0.1
#' res <- cfsens_cf_pac(data$X, data$Y, data$T, 
#'                      Gamma = Gamma, alpha = alpha,
#'                      delta = delta,  
#'                      ps_fun = regression_forest, 
#'                      pred_fun = quantile_forest)
#'
#' out_res <- predict(res, test$X, estimand = "treated", type = "ate")
#' cover <- mean((test$Y1 >= out_res$Y_lo) * (test$Y1 <= out_res$Y_hi))
#' cover
#' }
#'
#'
#'
#' @export

cfsens_cf_pac <- function(X, Y, T, 
                          Gamma, alpha, delta, 
                          side = c("two", "above", "below"),
                          score_type = c("cqr"),
                          ps_fun = regression_forest, 
                          ps = NULL, 
                          pred_fun = quantile_forest, 
                          train_pop = 0.75, train_id = NULL){
  
  ## Process the input
  side <- side[1]
  stopifnot(side %in% c("two", "above", "below"))
  score_type <- score_type[1]
  stopifnot(score_type %in% c("cqr"))

  ## Split the data into a training fold and a calibration fold
  n <- length(Y)
  if(is.null(train_id)){
    ntrain <- floor(n * train_pop)
    train_id <- sample(1:n, ntrain, replace = FALSE)
  }
  calib_id <- (1:n)[-train_id]

  X_train <- X[train_id,]
  T_train <- T[train_id]
  Y_train <- Y[train_id]
  ids1_train <- which(T_train == 1)

  X_calib <- X[calib_id,]
  T_calib <- T[calib_id]
  Y_calib <- Y[calib_id]
  ids1_calib <- which(T_calib == 1)

  ## Train the model of propensity scores on the training fold
  ps_mdl <- NULL
  if(is.null(ps)){
    ps_mdl <- ps_fun(X_train, T_train, num.threads = 1)
    ps <- predict(ps_mdl, X_calib)
  }

  ps0 <- ps[-ids1_calib,1]
  ps1 <- ps[ids1_calib,1]

  ## Estimate P(T=1) and P(T=0) 
  p1 <- mean(T_train == 1)
  p0 <- mean(T_train == 0)

  ## Train the prediction function on the training fold
  pred_mdl0 <- pred_fun(X_train[-ids1_train,], Y_train[-ids1_train], num.threads = 1)
  pred_mdl1 <- pred_fun(X_train[ids1_train,], Y_train[ids1_train], num.threads = 1)

  ## Compute the nonconformity scores on the calibration fold
  if(score_type == "cqr"){

    q_lo0 <- predict(pred_mdl0, X_calib[-ids1_calib,], alpha)
    q_hi0 <- predict(pred_mdl0, X_calib[-ids1_calib,], 1 - alpha)

    q_lo1 <- predict(pred_mdl1, X_calib[ids1_calib,], alpha)
    q_hi1 <- predict(pred_mdl1, X_calib[ids1_calib,], 1 - alpha)

    if(side == "two"){
      score0 <- pmax(q_lo0 - Y_calib[-ids1_calib], Y_calib[-ids1_calib] - q_hi0)
      score1 <- pmax(q_lo1 - Y_calib[ids1_calib], Y_calib[ids1_calib] - q_hi1) 
    }

    if(side == "above"){
      score0 <- q_lo0 - Y_calib[-ids1_calib] 
      score1 <- q_lo1 - Y_calib[ids1_calib]
    }

    if(side == "below"){
      score0 <- Y_calib[-ids1_calib] - q_hi0
      score1 <- Y_calib[ids1_calib] - q_hi1
    }
  }
  
  ## Attach the results to the output object
  obj <- list(ps_mdl = ps_mdl, 
              pred_mdl0 = pred_mdl0, 
              pred_mdl1 = pred_mdl1, 
              ps0 = ps0, ps1 = ps1,
              score0 = score0,
              score1 = score1,
              Gamma = Gamma, 
              alpha = alpha,
              delta = delta,
              side = side,
              p0 = p0, p1 = p1,
              score_type = score_type)

  class(obj) <- "cfpac"
  return(obj)

}

#' Predictive interval for cfpac objects
#'
#' Obtains predictive intervals on a test dataset based on 
#' an \code{cfpac} object from \code{link{cfsens_cf_pac}}.
#'
#' @details When \code{type = "treated"}, predictive intervals for Y(1)
#' are constructed; when \code{type = "control"}, predictive intervals
#' for Y(0) are constructed. 
#'
#' When \code{type = "ate"}, the inference is valid unconditionally;
#' when \code{type = "att"}, the inference is valid conditional on T=1, and 
#' \code{Y1_test} should be provided;
#' when \code{type = "atc"}, the inference is valid conditional on T=0, and
#' \code{Y0_test} should be provided.
#'
#' 
#' @param obj an object of class \code{cfpac}.
#' @param X_test testing covariates.
#' @param estimand the inferential target that 
#'                 takes value in \{"treated", "control"\}. See details.
#' @param type the type of inference target. 
#'             Takes value in \{"ate", "att", "atc"\}. See details.
#'
#' @return predictive intervals. A data.frame that contains \code{nrow(X_test)} rows and 
#'                               two columns: "Y_hi" refers the upper bound and "Y_lo" 
#'                               the lower bound.
#'
#' @export
predict.cfpac <- function(obj, X_test, 
                          estimand = c("treated", "control"), 
                          type = c("ate", "att", "atc"),
                          bnd_type = c("wsr", "hoeffding"), 
                          num_grids = 1000, rand_ind = NULL,
                          wsr_seed = 24601){
  ## Process the input
  estimand <- estimand[1]
  stopifnot(estimand %in% c("treated", "control"))
  type <- type[1]
  stopifnot(type %in% c("ate", "att", "atc"))
  bnd_type <- bnd_type[1]
  stopifnot(bnd_type %in% c("wsr", "hoeffding"))

  ps_mdl <- obj$ps_mdl
  pred_mdl0 <- obj$pred_mdl0
  pred_mdl1 <- obj$pred_mdl1
  ps0 <- obj$ps0
  ps1 <- obj$ps1
  score0 <- obj$score0
  score1 <- obj$score1
  Gamma <- obj$Gamma
  alpha <- obj$alpha
  delta <- obj$delta
  side <- obj$side
  p0 <- obj$p0
  p1 <- obj$p1
  score_type <- obj$score_type

  if(estimand == "treated"){
    pred_mdl <- pred_mdl1
    score <- score1
    ps <- ps1
  }else{
    pred_mdl <- pred_mdl0
    score <- score0
    ps <- ps0
  }

  ## Determine the dimension 
  n_calib <- length(score)
  if(is.null(nrow(X_test))){
    n_test <- length(X_test)
    X_test <- data.frame(X = X_test)
  }else{
    n_test <- nrow(X_test)
  }

  ## Compute the likelihood ratio bounds for the calibration fold
  bnds <- sapply(ps, lr_bnds, estimand = estimand, type = type, 
               Gamma = Gamma, p0 = p0, p1 = p1)
  lx_calib <- unlist(bnds[1,])
  ux_calib <- unlist(bnds[2,])
  M <- max(ux_calib)

  ## Assemble the calibration data
  data_calib <- data.frame(score = score, lx = lx_calib, ux = ux_calib)
  data_calib <- data_calib[order(data_calib$score),]

  ## Compute the cutoff for nonconformity scores
  if(bnd_type == "wsr"){

    if(is.null(rand_ind)){
      set.seed(wsr_seed)
      rand_ind <- sample(1:n_calib)
    }
    score_cutoff <- wsr_qtl(data_calib, delta, alpha, M, 
                      num_grids, rand_ind)
  }

  ## Generating the predictive interval
  if(score_type == "cqr"){
    
    Y_pred_lo <- predict(pred_mdl, X_test, alpha)
    Y_pred_hi <- predict(pred_mdl, X_test, 1-alpha)

    if(side == "two"){
      Y_lo <- Y_pred_lo - score_cutoff
      Y_hi <- Y_pred_hi + score_cutoff
    }

    if(side == "above"){
      Y_lo <- Y_pred_lo - score_cutoff 
      Y_hi <- NULL
    }

    if(side == "below"){
      Y_lo <- NULL
      Y_hi <- Y_pred_hi + score_cutoff
    }
 
  }

  return(list(Y_lo = Y_lo, Y_hi = Y_hi, score_cutoff = score_cutoff))
}
