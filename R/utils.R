#' Compute the bounds for the likelihood ratio under the marginal sensitivity model
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
  if(estimand == "control" & type == "atc"){
    lx <- 1
    ux <- 1
  }

  return(list(lx = lx, ux = ux))
}

#' Find the cutoff k* in the marginal conformal precedure
#'

find_kstar <- function(score, lx, ux, ux_test, alpha){

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

  if(kstar != Inf){
    score_cutoff <- score[score_order[kstar]]
  }else{
    score_cutoff <- Inf
  }

  return(list(kstar = kstar, score_cutoff = score_cutoff))

}

#' Computing the lower bounds for G(t) via the Waudby-Smith and Ramdas inequality
#'

## Compute lower bound for E[l(x)1(V<=t)]
wsr_lower <- function(delta, x, num_grids){

  n <- length(x)
  mu_hat <- (1 / 2 + cumsum(x)) / (1:n + 1)
  sig_hat <- (1 / 4 + cumsum((x - mu_hat)^2)) / (1 : n + 1)
  nu <- pmin(1, sqrt(2 * log(1 / delta) / (n * sig_hat^2)))
  nu[2 : length(nu)] <- nu[1 : (length(nu) - 1)]
  nu[1] <- pmin(1, sqrt(2 * log(1 / delta) / (n / 4)))


  u_list <- (1:num_grids) / num_grids
  k_all <- sapply(u_list, FUN = compute_k_lower, x = x, nu = nu)
  u_ind <- min(which(k_all <= 1 / delta))
  if (u_ind == Inf){
    u_ind <- num_grids
  }
  bnd <- u_list[u_ind]

  return(bnd)

}

## Compute max K(g)
compute_k_lower <- function(x, g, nu){
  kterms <- 1 + nu * (x - g)
  ks <- rep(kterms[1], length(kterms))
  for (ii in 2:length(ks)){
    ks[ii] <- ks[ii-1] * kterms[ii]
  }
  return(max(ks))
}


## Lower bound for G(t) for a single value of t
wsr_cdf_single <- function(delta, lx, ux, M, i, 
                           num_grids, rand_ind = NULL){
  n <- length(lx)
  lxx <- c(lx[1:i], rep(0, n-i)) / M
  uxx <- 1 + (c(rep(0,i), -ux[(i+1):n])) / M
  if (!is.null(rand_ind)){
    wsr1 <- wsr_lower(delta / 2, lxx[rand_ind], num_grids) * M
    wsr2 <- wsr_lower(delta / 2, uxx[rand_ind], num_grids) * M + 1 - M
  }else{
    wsr1 <- wsr_lower(delta / 2, sample(lxx), num_grids) * M
    wsr2 <- wsr_lower(delta / 2, sample(uxx), num_grids) * M + 1 - M
  }
  
  return(pmax(wsr1, wsr2))

}

wsr_cdf <- function(delta, lx, ux, M, 
                    num_grids, rand_ind = NULL){

  n <- length(lx)
  w_all <- sapply(1:n, FUN = wsr.cdf.single, 
                  delta = delta, lx = lx, ux = ux, 
                  M = M, num_grids = num_grids, rand_ind = rand_ind)
  return(w_all)

}

#' Bi-search the cutoff for nonconformity scores 
#'

wsr_qtl <- function(data_calib, delta, alpha, M,
                    num_grids, rand_ind){

  if(is.null(nrow(data_calib))){
    n_calib <- length(data_calib)
  }else{
    n_calib <- nrow(data_calib)
  }

  lx <- data_calib$lx
  ux <- data_calib$ux

  ## Initialization
  l_i <- 1
  r_i <- n_calib - 1
  m_i <- floor((l_i + r_i) / 2)

  left_cdf <- wsr_cdf_single(delta, lx, ux, M, 
                             l_i, num_grids, rand_ind)
  right_cdf <- wsr_cdf_single(delta, lx, ux, M, 
                              r_i, num_grids, rand_ind)
  mid_cdf <- wsr_cdf_single(delta, lx, ux, M, 
                            m_i, num_grids, rand_ind)

  gap <- min(m_i - l_i, r_i - m_i)

  ## Bi-search the cutoff
  while(gap > 0){

    if (mid_cdf < 1 - alpha){
      l_i <- m_i
      m_i <- floor((l_i + r_i) / 2)
      left_cdf <- mid_cdf
      mid_cdf <- wsr_cdf_single(delta, lx, ux, M, m_i, 
                                num_grids, rand_ind)
    }else{
      r_i <- m_i
      m_i <- floor((l_i + r_i) / 2)
      right_cdf <- mid_cdf
      mid_cdf <- wsr_cdf_single(delta, lx, ux, M, m_i, 
                                num_grids, rand_ind)
    }

    gap <- min(m_i - l_i, r_i - m_i)
    
  }
  
  if (mid_cdf >= 1 - alpha){
    return(data_calib$score[m_i])
  }else{
    return(data_calib$score[r_i])
  }

}
