#' Robust conformal counterfactual inference with PAC-type guarantee
#'
#'
#' @export

cfsens_cf_pac <- function(X, Y, T, 
                          Gamma, alpha, 
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
  n <- length(Y)
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
  
  ## Attach the results to the output object
  obj <- list(ps_mdl = ps_mdl, 
              pred_mdl = pred_mdl, 
              ps = ps, score = score, 
              Gamma = Gamma, alpha = alpha,
              p0 = p0, p1 = p1,
              score_type = score_type)

  class(obj) <- "cfpac"
  return(obj)

}

#' Predictive interval for cfpac objects
#'
#' @export
predict.cfpac <- function(obj, X_test, delta = 0.1, 
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
  pred_mdl <- obj$pred_mdl
  ps <- obj$ps
  score <- obj$score
  Gamma <- obj$Gamma
  alpha <- obj$alpha
  p0 <- obj$p0
  p1 <- obj$p1
  score_type <- obj$score_type


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

    Y_lo <- Y_pred_lo - score_cutoff
    Y_hi <- Y_pred_hi + score_cutoff
 
  }

  return(list(Y_lo = Y_lo, Y_hi = Y_hi))
}
