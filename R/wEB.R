#' Obtain the EB-rMAP Weight
#'
#' @description
#' Obtain the EB-rMAP weight corresponding to an observed response in current trial
#'
#' @param obj "EBrMAP" class object returned by `EB_rMAP` function
#' @param y_obs Observed response in the current trial
#'
#' @details
#' `y_obs` for binomial and Poisson outcomes must be an integer.
#'
#' In the case of normal outcome, `y_obs` is matched to the closest break in `y_range` argument in `EB_rMAP` function, whose corresponding EB-rMAP weight is returned. Rounding issues may cause error occasionally. In this case, jittering `y_obs` might help. For example, try `wEB(obj, 30.001)` if `wEB(obj, 30)` fails.
#'
#' @returns
#' A scalar of EB-rMAP weight
#'
#' @export

wEB <- function(obj, y_obs){
  if(obj$dist=="beta"){
    w <- as.numeric(obj$wdt %>% filter(y==y_obs) %>% select(w_eb))
    return(w)
  }
  if(obj$dist=="norm"){
    w <- as.numeric(obj$wdt %>% mutate(ad=abs(y-y_obs)) %>% filter(ad==min(ad)) %>% select(w_eb))
    return(w)
  }
  if(obj$dist=="gamma"){
    w <- as.numeric(obj$wdt %>% filter(y==y_obs) %>% select(w_eb))
    return(w)
  }
}
