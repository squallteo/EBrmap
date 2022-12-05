#' @param obj "EBrMAP" class object returned by EB_rMAP function
#' @param y_obs Observed response in the current trial

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
