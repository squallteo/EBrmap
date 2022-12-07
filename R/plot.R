#' Plot an EBrMAP Object
#'
#' @description
#' Plot EB-rMAP weight against various observed response in the current trial
#'
#' @param x "EBrMAP" class object returned by `EB_rMAP` function
#'
#' @returns
#' A ggplot2 object which can be further customized
#'
#' @export

plot.EBrMAP <- function(x, ...){
  if(x$dist=="beta"){
    p <-
    x$wdt %>%
      ggplot(aes(x=y, y = w_eb)) + geom_line(size=1) +
      xlab("Observed Number of Events") + ylab("EB-rMAP Weight") + theme_bw() +
      ggtitle(paste("Endpoint: Binomial; Sample size (N): ", x$n, "; Gamma: ", x$ppp_cut, sep=""))
    return(p)
  }

  if(x$dist=="norm"){
    p <-
      x$wdt %>%
      ggplot(aes(x=y, y = w_eb)) + geom_line(size=1) +
      scale_x_continuous(breaks = seq(min(x$wdt$y), max(x$wdt$y), length=11)) +
      xlab("Observed Mean Response") + ylab("EB-rMAP Weight") + theme_bw() +
      ggtitle(paste("Endpoint: Normal; Sample size (N): ", x$n, "; Gamma: ", x$ppp_cut, sep=""))
    return(p)
  }

  if(x$dist=="gamma"){
    p <-
      x$wdt %>%
      ggplot(aes(x=y, y = w_eb)) + geom_line(size=1) +
      xlab("Observed Number of Events") + ylab("EB-rMAP Weight") + theme_bw() +
      ggtitle(paste("Endpoint: Poisson; Sample size (Total Exposure): ", x$n, "; Gamma: ", x$ppp_cut, sep=""))
    return(p)
  }
}

