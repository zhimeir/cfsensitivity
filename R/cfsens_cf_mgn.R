#' Robust conformal counterfactual inference with marginal guarantee
#'
#' \code{cfsens_cf_mgn} constructs robust predictive intervals for 
#' counterfactuals with the margianl guarantee.
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
#' @param side the type of predictive intervals that takes value in \{"two", "above", "below"\}. See details.
#' @param score_type the type of nonconformity scores. The default is "cqr".
#' @param ps_fun a function that models the treatment assignment mechanism. The default is "regression_forest".
#' @param ps a vector of propensity score. The default is \code{NULL}. 
#' @param pred_fun a function that models the potential outcome conditional on the covariates. The default is "quantile_forest".
#' @param train_prop proportion of units used for training. The default is 75\%. 
#' @param train_id The index of the units used for training. The default is \code{NULL}.
#'
#' @return an \code{cfmgn} object.
#'
#' @seealso 
#' \code{\link{cfsens_cf_pac}}
#'
#' @examples
#'
#'
#'
#' @export

cfsens_cf_mgn <- function(X, Y, T, 
                          Gamma, alpha,
                          side = c("two", "above", "below"),
                          score_type = c("cqr"),
                          ps_fun = regression_forest, 
                          ps = NULL,
                          pred_fun = quantile_forest,
                          train_prop = 0.75, train_id = NULL){
  
  ## Process the input
  score_type <- score_type[1]
  stopifnot(score_type %in% c("cqr"))
  side <- side[1]
  stopifnot(side %in% c("two", "above", "below"))

  ## Split the data into a training fold and a calibration fold
  n <- dim(X)[1]
  if(is.null(train_id)){
    ntrain <- floor(n * train_prop)
    train_id <- sample(1:n, ntrain, replace = FALSE)
  }
  calib_id <- (1:n)[-train_id]

  X_train <- X[train_id,]
  T_train <- T[train_id]
  Y_train <- Y[train_id]
  inds1_train <- which(T_train == 1)

  X_calib <- X[calib_id,]
  T_calib <- T[calib_id]
  Y_calib <- Y[calib_id]
  inds1_calib <- which(T_calib == 1)

  ## Train the model of propensity scores on the training fold
  ps_mdl <- NULL
  if(is.null(ps)){
    ps_mdl <- ps_fun(X_train, T_train, num.threads = 1)
    ps <- predict(ps_mdl, X_calib)
  }

  ps0 <- ps[-inds1_calib,1]
  ps1 <- ps[inds1_calib,1]

  ## Estimate P(T=1) and P(T=0) 
  p1 <- mean(T_train == 1)
  p0 <- mean(T_train == 0)

  ## Train the prediction function on the training fold
  pred_mdl0 <- pred_fun(X_train[-inds1_train,], Y_train[-inds1_train], num.threads = 1)
  pred_mdl1 <- pred_fun(X_train[inds1_train,], Y_train[inds1_train], num.threads = 1)

  ## Compute the nonconformity scores on the calibration fold
  if(score_type == "cqr"){
    q_lo0 <- predict(pred_mdl0, X_calib[-inds1_calib,], alpha)
    q_hi0 <- predict(pred_mdl0, X_calib[-inds1_calib,], 1 - alpha)
    q_lo1 <- predict(pred_mdl1, X_calib[inds1_calib,], alpha)
    q_hi1 <- predict(pred_mdl1, X_calib[inds1_calib,], 1 - alpha)

    ## If the predictive interval is two-sided: [a,b]
    if(side == "two"){      
      score0 <- pmax(q_lo0 - Y_calib[-inds1_calib], Y_calib[-inds1_calib] - q_hi0)
      score1 <- pmax(q_lo1 - Y_calib[inds1_calib], Y_calib[inds1_calib] - q_hi1)
    }

    ## If the predictive interval is one-sided: [a,Inf)
    if(side == "above"){
      score0 <- q_lo0 - Y_calib[-inds1_calib]
      score1 <- q_lo1 - Y_calib[inds1_calib]
    }

    ## If the predictive interval is one-sided: (-Inf,a]
    if(side == "below"){
      score0 <- Y_calib[-inds1_calib] - q_hi0
      score1 <- Y_calib[inds1_calib] - q_hi1
    }

  }
  
  ## Attach the results to the return object
  obj <- list(ps_mdl = ps_mdl, 
              pred_mdl0 = pred_mdl0, 
              pred_mdl1 = pred_mdl1, 
              score0 = score0, score1 = score1, 
              ps0 = ps0, ps1 = ps1, 
              Gamma = Gamma,
              p0 = p0, p1 = p1,
              score_type = score_type, 
              alpha = alpha, side = side)

  class(obj) <- "cfmgn"
  return(obj)
}


#' Prediction method for cfmgn objects
#'
#' Obtains predictive intervals on a test dataset based on 
#' an \code{cfmgn} object from \code{link{cfsens_cf_mgn}}.
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
#' @param obj an object of class \code{cfmgn}.
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

predict.cfmgn <- function(obj, X_test,  
                          estimand = c("treated", "control"),
                          type = c("ate", "att", "atc")){

  ## Process the input
  type <- type[1]
  estimand <- estimand[1]

  ps_mdl <- obj$ps_mdl
  pred_mdl0 <- obj$pred_mdl0
  pred_mdl1 <- obj$pred_mdl1
  ps0 <- obj$ps0
  ps1 <- obj$ps1
  score0 <- obj$score0
  score1 <- obj$score1
  Gamma <- obj$Gamma
  p0 <- obj$p0
  p1 <- obj$p1
  score_type <- obj$score_type
  alpha <- obj$alpha
  side <- obj$side

  ## Determine the score to use
  if(estimand == "treated"){
    pred_mdl <- pred_mdl1
    score <- score1
    ps <- ps1
  }else{
    pred_mdl <- pred_mdl0
    score <- score0
    ps <- ps0
  }

  ## Determine the dimension of X_test
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

  ## Compute the propensity scores for the test fold
  ps_test <- predict(ps_mdl, X_test)

  ## Compute the likelihood ratio bounds for the test fold
  bnds <- sapply(ps_test, lr_bnds, estimand = estimand, type = type, 
                 Gamma = Gamma, p0 = p0, p1 = p1)
  lx_test <- unlist(bnds[1,])
  ux_test <- unlist(bnds[2,])

  ## Find the cutoff kstar
  score_cutoff <- rep(NA, n_test)
  for(i in 1:n_test){
    score_cutoff[i] <- find_kstar(score, lx_calib, ux_calib, 
                                  ux_test[i], alpha)$score_cutoff
  } 

  ## Apply the predictive model to the test fold and
  ## invert the non-conformity score to predictive intervals 
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

  return(list(Y_lo = Y_lo, Y_hi = Y_hi))
}















