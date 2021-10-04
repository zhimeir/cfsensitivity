#' Robust conformal inference for counterfactuals
#' 
#' 
#'
#' @export

cfsens_cf <- function(X, Y, T, Gamma,
                      estimand = c("ate", "att", "atc"),
                      guarantee = c("mgn", "pac"),
                      score = c("cqr", "cmr", "cdr"),
                      cf_type = c("split"),
                      trainprop = 0.75, trainid = NULL){

  ## Check the input format
  type <- type[1]
  guarantee <- guarantee[1]
  cf_type <- cf_type[1]

  if(guarantee == "mgn" & cf_type == "split"){
      
  }

  if(guarantee == "pac" & cf_type == "split"){
  }
  


}
