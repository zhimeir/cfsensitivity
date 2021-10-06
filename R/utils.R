#' Compute the bounds for the likelihood ratio 
#' under the marginal sensitivity model
#'
#' @export
lr_bnds <- function(estimand = c("treated", "control"), 
                    type = c("ate", "att", "atc"),
                    Gamma, ps, p1, p0){

  ## Preprocess the input
  estimand <- estimand[1]
  type <- type[1]
  
  ## Compute the bounds
  ## Y(1)
  if(estimand == "treated" & type == "ate"){
    lx <- p1 * (1 + (1 - ps) / ps / Gamma)
    ux <- p1 * (1 + (1 - ps) / ps * Gamma)
  }  

  ##  Y(1) | T = 1
  if(estimand == "treated" & type == "att"){
    lx <- 1
    ux <- 1
  }

  ## Y(1) | T = 0
  if(estimand == "treated" & type == "atc"){
    lx <- (p1 / p0) * ((1 - ps) / ps / Gamma)
    ux <- (p1 / p0) * ((1 - ps) / ps * Gamma)
  }

  ## Y(0)
  if(estimand == "control" & type == "ate"){
    lx <- p0 * (1 + ps / (1 - ps) / Gamma)
    ux <- p0 * (1 + ps / (1 - ps) * Gamma)
  }

  ## Y(0) | T = 1
  if(estimand == "control" & type == "att"){
    lx <- (p0 / p1) * (ps / (1 - ps) / Gamma)
    ux <- (p0 / p1) * (ps / (1 - ps) * Gamma)
  }

  ## Y(0) | T = 0
  if(estimand == "control" & type == "att"){
    lx <- 1
    ux <- 1
  }

  return(list(lx = lx, ux = ux))
}

#' Find the cutoff k* in the marginal conformal precedure
#'

find_kstar <- function(score, lx, ux, 
                       ux_test, alpha){

  ns <- length(score)

  ## Order the bounds according to the order of the non-conformity scores
  score_order <- order(score)
  lx <- lx[score_order]
  ux <- ux[score_order]

  ## Compute the partial sums in the numerator and the denominator
  sum_num <- rep(NA, ns)
  sum_den <- rep(NA, ns)
  sum_num[1] <- lx[1]
  sum_den[1] <- lx[1] + sum(ux[2:ns])

  for(k in 2:ns){
    sum_num[k] <- sum_num[k-1] + lx[k]
    sum_den[k] <- sum_den[k-1] - ux[k] + lx[k]
  }

  ## Find k*
  ratio <- sum_num / (sum_den + ux_test)
  kstar <- min(which(ratio >= 1 - alpha))
  score_cutoff <- score[score_order[kstar]]
  return(list(kstar = kstar, score_cutoff = score_cutoff))

}
