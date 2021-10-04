#' Robust conformal counterfactual inference with marginal guarantee
#'
#'
#' @export

cfsens_cf_mgn <- function(X, Y, T, 
                          Gamma,alpha,
                          estimand = c("ate", "att", "atc"),
                          score = c("cqr", "cmr", "cdr"),
                          ps_fun = c("regression_forest"), 
                          ps_fun_args = NULL,
                          ps = NULL,
                          pred_fun = c("quantile_forest"),
                          pred_fun_args = NULL,
                          train_pop = 0.75, train_id = NULL){
  
  ## Process the input
  estimand <- estimand[1]
  stopifnot(estimand %in% c("ate", "att", "atc"))
  score <- score[1]
  stopifnot(score %in% c("cqr", "cmr", "cdr"))

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
  
  ## Train the prediction function on the training fold
  pred_mdl <- pred_fun(X_train, Y_train, num.threads = 1)

  ## Compute the nonconformity scores on the calibration fold
  if(score == "cqr"){
    q_lo <- predict(pred_mdl, X_calib, alpha)
    q_hi <- predict(pred_mdl, X_calib, 1 - alpha)
    score <- pmax(q_lo - Y_calib, Y_calib - q_hi)
  }
  
  ## Attach the results to the return object
  obj$ps_mdl <- ps_mdl
  obj$pred_mdl <- pred_mdl
  obj$ps <- ps
  obj$score <- score

  class(obj) <- cfmgn

}


#' Predictive interval for cfmgn objects
#'
#' @export

predict.cfmgn <- function(obj, Xtest, Ttest){

}
