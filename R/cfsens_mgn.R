#' Sensitivity analysis of individual treatment effect with marginal guarantee
#'
#'
#' @export

cfsens_mgn <- function(X, Y, T,  
                       alpha,
                       score_type = c("cqr", "cmr", "cdr"),
                       ps_fun = regression_forest, 
                       ps_fun_args = NULL,
                       ps = NULL,
                       pred_fun = quantile_forest,
                       pred_fun_args = NULL,
                       train_pop = 0.75, train_id = NULL){
  
  ## Process the input
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
              p0 = p0, p1 = p1,
              score_type = score_type)

  class(obj) <- "itemgn"
  return(obj)
}


#' Sensitivity analysis for itemgn objects
#'
#' @export

predict.itemgn <- function(obj, X_test, Y1_test, Y0_test,
                           alpha,
                           type = c("ate", "att", "atc"),
                           null_type = c("sharp", "directional"), 
                           Gamma_max = 10){

  ## Process the input
  type <- type[1]
  stopifnot(type %in% c("ate", "att", "atc"))
  null_type <- null_type[1]
  stopifnot(null_type %in% c("sharp", "directional"))

  ps_mdl <- obj$ps_mdl
  pred_mdl <- obj$pred_mdl
  ps <- obj$ps
  score <- obj$score
  p0 <- obj$p0
  p1 <- obj$p1
  score_type <- obj$score_type
  Gamma_grid <- seq(0.5, Gamma_max, by = 0.5)
  
  for (Gamma in Gamma_grid){
    gamma_obj <- list(ps_mdl = ps_mdl, pred_mdl = pred_mdl, 
                      ps = ps, score = score, 
                      p0 = p0, p1 = p1, 
                      score_type =  score_type)
    class(gamma_obj) <- "cfmgn"

    res_treated <- preduct(gamma_obj, X_test, alpha, "treated", type)
    res_control <- preduct(gamma_obj, X_test, alpha, "control", type)

  } 






  return(list(Y_lo = Y_lo, Y_hi = Y_hi))
}















