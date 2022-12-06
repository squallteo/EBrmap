#' Plot an EBrMAP Object
#'
#' @description
#' Plot EB-rMAP weight against various observed response in the current trial
#'
#' @param obj "EBrMAP" class object returned by `EB_rMAP` function
#'
#' @returns
#' A ggplot2 object which can be further customized
#'
#' @export

plot.EBrMAP <- function(obj){
  if(obj$dist=="beta"){
    p <-
    obj$wdt %>%
      ggplot(aes(x=y, y = w_eb)) + geom_line(size=1) +
      xlab("Observed Number of Events") + ylab("EB-rMAP Weight") + theme_bw() +
      ggtitle(paste("Endpoint: Binomial; Sample size (N): ", obj$n, "; Gamma: ", obj$ppp_cut, sep=""))
    return(p)
  }

  if(obj$dist=="norm"){
    p <-
      obj$wdt %>%
      ggplot(aes(x=y, y = w_eb)) + geom_line(size=1) +
      scale_x_continuous(breaks = seq(min(obj$wdt$y), max(obj$wdt$y), length=11)) +
      xlab("Observed Mean Response") + ylab("EB-rMAP Weight") + theme_bw() +
      ggtitle(paste("Endpoint: Normal; Sample size (N): ", obj$n, "; Gamma: ", obj$ppp_cut, sep=""))
    return(p)
  }

  if(obj$dist=="gamma"){
    p <-
      obj$wdt %>%
      ggplot(aes(x=y, y = w_eb)) + geom_line(size=1) +
      xlab("Observed Number of Events") + ylab("EB-rMAP Weight") + theme_bw() +
      ggtitle(paste("Endpoint: Poisson; Sample size (Total Exposure): ", obj$n, "; Gamma: ", obj$ppp_cut, sep=""))
    return(p)
  }
}

