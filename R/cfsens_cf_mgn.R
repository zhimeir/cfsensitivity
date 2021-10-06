#' Robust conformal counterfactual inference with marginal guarantee
#'
#'
#' @export

cfsens_cf_mgn <- function(X, Y, T, 
                          Gamma,alpha,
                          estimand = c("ate", "att", "atc"),
                          score_type = c("cqr", "cmr", "cdr"),
                          ps_fun = regression_forest, 
                          ps_fun_args = NULL,
                          ps = NULL,
                          pred_fun = quantile_forest,
                          pred_fun_args = NULL,
                          train_pop = 0.75, train_id = NULL){
  
  ## Process the input
  estimand <- estimand[1]
  stopifnot(estimand %in% c("ate", "att", "atc"))
  score_type <- score_type[1]
  stopifnot(score_type %in% c("cqr", "cmr", "cdr"))

  ## Split the data into a training fold and a calibration fold
  n <- dim(X)[1]
  if(is.null(train_id)){
    ntrain <- floor(n * train_pop)
    train_id <- sample(1:n, ntrain, replace = FALSE)
  }
  calib_id <- (1:n)[-train_id]

  X_train <- X[train_id,]
  T_train <- T[train_id]
  Y_train <- Y[train_id]

  X_calib <- X[calib_id,]
  T_calib <- T[calib_id]
  Y_calib <- Y[calib_id]

  ## Train the model of propensity scores on the training fold
  ps_mdl <- NULL
  if(is.null(ps)){
    ps_mdl <- ps_fun(X_train, T_train, num.threads = 1)
    ps <- predict(ps_mdl, X_calib)
  }

  ## Compute P(T=1) and P(T=0) 
  p1 <- mean(T_train == 1)
  p0 <- mean(T_train == 0)

  ## Train the prediction function on the training fold
  pred_mdl <- pred_fun(X_train, Y_train, num.threads = 1)

  ## Compute the nonconformity scores on the calibration fold
  if(score_type == "cqr"){
    q_lo <- predict(pred_mdl, X_calib, alpha)
    q_hi <- predict(pred_mdl, X_calib, 1 - alpha)
    score <- pmax(q_lo - Y_calib, Y_calib - q_hi)
  }
  
  ## Attach the results to the return object
  obj <- list(ps_mdl = ps_mdl, 
              pred_mdl = pred_mdl, 
              ps = ps, score = score, 
              Gamma = Gamma,
              p0 = p0, p1 = p1,
              score_type = score_type)

  class(obj) <- "cfmgn"
  return(obj)
}


#' Predictive interval for cfmgn objects
#'
#' @export

predict.cfmgn <- function(obj, X_test, alpha = 0.1, 
                          estimand = c("treated", "control"),
                          type = c("ate", "att", "atc")){

  ## Process the input
  type <- type[1]
  estimand <- estimand[1]

  ps_mdl <- obj$ps_mdl
  pred_mdl <- obj$pred_mdl
  ps <- obj$ps
  score <- obj$score
  Gamma <- obj$Gamma
  p0 <- obj$p0
  p1 <- obj$p1
  score_type <- obj$score_type


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

    Y_lo <- Y_pred_lo - score_cutoff
    Y_hi <- Y_pred_hi + score_cutoff

  }

  return(list(Y_lo = Y_lo, Y_hi = Y_hi))
}















