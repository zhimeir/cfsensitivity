#' Sensitivity analysis of individual treatment effect with marginal guarantee
#'
#' \code{cfsens_mgn} conducts sensitivity analysis of individual treatment with marginal guarantee.
#' It currently supports testing the sharp null and the directional nulls.
#'
#' @param X covariates.
#' @param Y the observed outcome vector.
#' @param T the vector of treatment assignments.
#' @param alpha the target confidence level.
#' @param null_type the null to be tested. Takes 
#'
#'
#' @export

cfsens_mgn <- function(X, Y, T,  
                       alpha,
                       null_type = c("sharp", "negative", "positive"), 
                       score_type = c("cqr"),
                       ps_fun = regression_forest, 
                       ps = NULL,
                       pred_fun = quantile_forest,
                       train_pop = 0.75, train_id = NULL){
  
  ## Process the input
  score_type <- score_type[1]
  stopifnot(score_type %in% c("cqr"))
  null_type <- null_type[1]
  stopifnot(null_type %in% c("sharp", "negative","positive"))

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

  ## Compute P(T=1) and P(T=0) 
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

    ## H_0: Y(1) - Y(0) = 0
    if(null_type == "sharp"){
      score0 <- pmax(q_lo0 - Y_calib[-ids1_calib], Y_calib[-ids1_calib] - q_hi0)
      score1 <- pmax(q_lo1 - Y_calib[ids1_calib], Y_calib[ids1_calib] - q_hi1)
    }

    ## H0: Y(1) - Y(0) <= 0
    if(null_type == "negative"){
      score0 <-  Y_calib[-ids1_calib] - q_hi0
      score1 <- q_lo1 - Y_calib[ids1_calib]
    }

    ## H0: Y(1) - Y(0) >=0
    if(null_type == "positive"){
      score0 <- q_lo0 - Y_calib[-ids1_calib]
      score1 <- Y_calib[ids1_calib] - q_hi1
    }
  }
  
  ## Attach the results to the return object
  obj <- list(ps_mdl = ps_mdl, 
              pred_mdl0 = pred_mdl0,
              pred_mdl1 = pred_mdl1, 
              ps0 = ps0, ps1 = ps1, 
              score0 = score0,
              score1 = score1, 
              p0 = p0, p1 = p1,
              score_type = score_type, 
              alpha = alpha, null_type = null_type)

  class(obj) <- "itemgn"
  return(obj)
}


#' Sensitivity analysis for itemgn objects
#'
#' @export

predict.itemgn <- function(obj, X_test, Y1_test = NULL, Y0_test = NULL,
                           type = c("ate", "att", "atc"),
                           Gamma_max = 5, gamma_length = 51){

  ## Process the input
  type <- type[1]
  stopifnot(type %in% c("ate", "att", "atc"))

  alpha <- obj$alpha
  ps_mdl <- obj$ps_mdl
  pred_mdl0 <- obj$pred_mdl0
  pred_mdl1 <- obj$pred_mdl1
  ps0 <- obj$ps0
  ps1 <- obj$ps1
  score0 <- obj$score0
  score1 <- obj$score1
  p0 <- obj$p0
  p1 <- obj$p1
  score_type <- obj$score_type
  null_type <- obj$null_type
  Gamma_grid <- seq(1, Gamma_max, length.out = gamma_length)
  
  ## Assemble the input and pass to the cf_mgn function 
  gamma_obj <- list(ps_mdl = ps_mdl, 
                    pred_mdl0 = pred_mdl0,
                    pred_mdl1 = pred_mdl1, 
                    ps0 = ps0, ps1 = ps1, 
                    score0 = score0,
                    score1 = score1,
                    p0 = p0, p1 = p1, 
                    Gamma = NULL,
                    score_type =  score_type, 
                    alpha = alpha)
  class(gamma_obj) <- "cfmgn"

  ind_cover <- c()
  for (Gamma in Gamma_grid){
    
    gamma_obj$Gamma <- Gamma
    
    if(type == "ate"){
      gamma_obj$alpha <- alpha / 2
      
      if(null_type == "sharp"){
        gamma_obj$side <- "two"
        res_treated <- predict(gamma_obj, X_test, "treated", type)
        res_control <- predict(gamma_obj, X_test, "control", type)
        ite_lo <- res_treated$Y_lo - res_control$Y_hi
        ite_hi <- res_treated$Y_hi - res_control$Y_hi
        ind_cover <- cbind(ind_cover, (ite_lo <=0) * (ite_hi >=0))
      }

      if(null_type == "negative"){
        gamma_obj$side <- "above"
        res_treated <- predict(gamma_obj, X_test, "treated", type)
        gamma_obj$side <- "below"
        res_control <- predict(gamma_obj, X_test, "control", type)
        ite_lo <- res_treated$Y_lo - res_control$Y_hi
        ind_cover <- cbind(ind_cover, ite_lo <=0)
      }

      if(null_type == "positive"){
        gamma_obj$side <- "below"
        res_treated <- predict(gamma_obj, X_test, "treated", type)
        gamma_obj$side <- "above"
        res_control <- predict(gamma_obj, X_test, "control", type)
        ite_hi <- res_treated$Y_hi - res_control$Y_lo
        ind_cover <- cbind(ind_cover, ite_hi >=0)
      }

    }

    if(type == "att"){
      if(null_type == "sharp"){
        gamma_obj$side <- "two"
        res_control <- predict(gamma_obj, X_test, "control", type)
        ite_lo <- Y1_test - res_control$Y_lo
        ite_hi <- Y1_test - res_control$Y_hi
        ind_cover <- cbind(ind_cover, (ite_lo <=0) * (ite_hi >=0))
      }

      if(null_type == "negative"){
        gamma_obj$side <- "below"
        res_control <- predict(gamma_obj, X_test, "control", type)
        ite_lo <- Y1_test - res_control$Y_hi
        ind_cover <- cbind(ind_cover, ite_lo <=0) 
      }

      if(null_type == "positive"){
        gamma_obj$side <- "above"
        res_control <- predict(gamma_obj, X_test, "control", type)
        ite_hi <- Y1_test - res_control$Y_lo
        ind_cover <- cbind(ind_cover, ite_hi >=0) 
      }

    }
    
    if(type == "atc"){
      if(null_type == "sharp"){
        gamma_obj$side <- "two"
        res_treated <- predict(gamma_obj, X_test, "treated", type)
        ite_lo <- res_treated$Y_lo - Y0_test
        ite_hi <- res_treated$Y_hi - Y0_test 
        ind_cover <- cbind(ind_cover, (ite_lo <=0) * (ite_hi >=0))
      }

      if(null_type == "negative"){
        gamma_obj$side <- "above"
        res_treated <- predict(gamma_obj, X_test, "treated", type)
        ite_lo <- res_treated$Y_lo - Y0_test
        ind_cover <- cbind(ind_cover, ite_lo <=0)
      }

      if(null_type == "positive"){
        gamma_obj$side <- "below"
        res_treated <- predict(gamma_obj, X_test, "treated", type)
        ite_hi <- res_treated$Y_hi - Y0_test
        ind_cover <- cbind(ind_cover, ite_hi >=0) 
      }
    }

  } 

  hat_gamma <- apply(ind_cover, 1, get_hat_gamma, Gamma_grid)
  return(hat_gamma)
}

get_hat_gamma <- function(ind_cover, Gamma_grid){
  if(ind_cover[length(Gamma_grid)] == 0){
    hat_gamma <- max(Gamma_grid)
  }else{
    hat_gamma <- Gamma_grid[min(which(ind_cover == 1))]
  }

  return(hat_gamma)
}

