#' Empirical Bayes Robust Meta-Analytical-Predictive Prior
#'
#' @description
#' Implement the empirical Bayes robust meta-analytical-predictive (EB-rMAP) prior for outcomes distributions including binomial, normal (with known standard deviation) and Poisson.
#'
#' @param map_obj The mixture of conjugate priors approximating the MAP prior, derived by the *RBesT* package.
#' @param vague_prior The vague prior to form the robust MAP prior, must be in the same class of `map_obj`.
#' @param ppp_cut A scalar cutoff for the prior predictive p value.
#' @param n Current trial sample size: number of subjects for binary and normal endpoints; total exposure for time-to-event endpoint.
#' @param y_range A vector containing the lower and upper limits of the observed response in the current trial. Required for normal and time-to-event endpoints.
#' @param nbreak The number of breaks within the `y_range`, inclusive of the lower and upper limits. Required for normal endpoint. See details. Default is 100.
#'
#' @details
#' For discrete outcomes (binomial and Poisson), all possible outcome values can be enumerated and thus the `nbreak` option is needed. Note that `y_range` is still needed for a Poisson outcome.
#'
#' For the normal outcome, the interval `y_range` is evenly divided into `nbreak` breaks. EB-rMAP weight is computed at each break value. Default `nbreak` is 100. A finer grid can be used with a larger `nbreak` at the cost of more computation time.
#'
#' @returns
#' An object of S3 class "EBrMAP"
#'
#' @examples
#' \dontrun{
#' set.seed(712)
#'
#' dt <- RBesT::crohn
#' sigma <- 88
#' dt$se_yh <- sigma/sqrt(dt$n)
#' map_mcmc <- RBesT::gMAP(cbind(y, se_yh) ~ 1 | study, weight = n, data=dt,
#'                     family=gaussian, beta.prior=cbind(0, sigma),
#'                     tau.dist="HalfNormal",tau.prior=cbind(0,sigma/2))
#' map_hat <- RBesT::automixfit(map_mcmc)
#' sigma(map_hat) <- sigma
#' vague_prior <- RBesT::mixnorm(c(1, -50, 1), sigma=sigma, param="mn")
#'
#' obj <- EB_rMAP(map_hat, vague_prior, ppp_cut=0.85, n=100, y_range=c(-80, -20))
#' wEB(obj, -61)
#' plot(obj)
#' }
#'
#' @export

EB_rMAP <- function(map_obj, vague_prior, ppp_cut, n, y_range=NA, nbreak=100){
  if(class(map_obj)[3]=="betaMix"){
    assertthat::assert_that(class(vague_prior)[1]=="betaMix",
                            msg = "The vague prior distribution must be Beta.")

    a <- c(vague_prior)[2]
    b <- c(vague_prior)[3]
    for(y in 0:n){
      w <- seq(0, 1, 0.01)
      ppp <- rep(NA, length(w))
      for(i in 1:length(w)){
        rmap <- RBesT::robustify(map_obj, weight=w[i], mean = a/(a+b), n = a+b-1)
        rmap_pred <- RBesT::preddist(rmap, n=n)
        p_lower <- RBesT::pmix(rmap_pred, y)
        ppp[i] <- ifelse(p_lower < 0.5, 2*p_lower, 2*(1-p_lower))
      }
      pppdt <- tibble(w, ppp) %>% mutate(pass=(ppp >= ppp_cut)) %>% filter(pass)
      w_eb <- ifelse(nrow(pppdt) == 0, 1, min(pppdt$w))
      ppp_eb <- as.numeric(pppdt %>% filter(w==w_eb) %>% select(ppp))

      if(y==0){wdt <- tibble(y, w_eb, ppp_eb)}
      else {wdt <- rbind(wdt, tibble(y, w_eb, ppp_eb))}
    }
    out_obj <- list(dist="beta", n=n, ppp_cut=ppp_cut, wdt=wdt)
    class(out_obj) <- "EBrMAP"
    return(out_obj)
  }

  if(class(map_obj)[3]=="normMix"){
    assertthat::assert_that(class(vague_prior)[1]=="normMix",
                            msg = "The vague prior distribution must be Normal.")
    assertthat::assert_that(!all(is.na(y_range)),
                            msg = "Must specify the lower and upper limits of the observed response in the current trial.")

    y_vec <- seq(y_range[1], y_range[2], length=nbreak)
    m <- c(vague_prior)[2]
    s <- c(vague_prior)[3]
    for(j in 1:length(y_vec)){
      y <- y_vec[j]
      w <- seq(0, 1, 0.01)
      ppp <- rep(NA, length(w))
      for(i in 1:length(w)){
        rmap <- RBesT::robustify(map_obj, weight=w[i], mean=m, n=1, sigma=s)
        rmap_pred <- RBesT::preddist(rmap, n=n, sigma=sigma(rmap))
        p_lower <- RBesT::pmix(rmap_pred, y)
        ppp[i] <- ifelse(p_lower < 0.5, 2*p_lower, 2*(1-p_lower))
      }
      pppdt <- tibble(w, ppp) %>% mutate(pass=(ppp >= ppp_cut)) %>% filter(pass)
      w_eb <- ifelse(nrow(pppdt) == 0, 1, min(pppdt$w))
      ppp_eb <- as.numeric(pppdt %>% filter(w==w_eb) %>% select(ppp))

      if(j==1){wdt <- tibble(y, w_eb, ppp_eb)}
      else {wdt <- rbind(wdt, tibble(y, w_eb, ppp_eb))}
    }
    out_obj <- list(dist="norm", n=n, ppp_cut=ppp_cut, wdt=wdt)
    class(out_obj) <- "EBrMAP"
    return(out_obj)
  }

  if(class(map_obj)[3]=="gammaMix"){
    assertthat::assert_that(class(vague_prior)[1]=="gammaMix",
                            msg = "The vague prior distribution must be Gamma.")
    assertthat::assert_that(!all(is.na(y_range)),
                            msg = "Must specify the lower and upper limits of the observed response in the current trial.")
    assertthat::assert_that(y_range[1]%%1==0 & y_range[2]%%1==0 & y_range[1] >=0,
                            msg = "The lower and upper limits must be non-negative integers.")

    m <- c(vague_prior)[2]
    s <- c(vague_prior)[3]
    for(y in y_range[1]:y_range[2]){
      w <- seq(0, 1, 0.01)
      ppp <- rep(NA, length(w))
      for(i in 1:length(w)){
        rmap <- RBesT::robustify(map_obj, weight=w[i], mean=m, n=s)
        rmap_pred <- RBesT::preddist(rmap, n=n)
        p_lower <- RBesT::pmix(rmap_pred, y)
        ppp[i] <- ifelse(p_lower < 0.5, 2*p_lower, 2*(1-p_lower))
      }
      pppdt <- tibble(w, ppp) %>% mutate(pass=(ppp >= ppp_cut)) %>% filter(pass)
      w_eb <- ifelse(nrow(pppdt) == 0, 1, min(pppdt$w))
      ppp_eb <- as.numeric(pppdt %>% filter(w==w_eb) %>% select(ppp))

      if(y==y_range[1]){wdt <- tibble(y, w_eb, ppp_eb)}
      else {wdt <- rbind(wdt, tibble(y, w_eb, ppp_eb))}
    }
    out_obj <- list(dist="gamma", n=n, ppp_cut=ppp_cut, wdt=wdt)
    class(out_obj) <- "EBrMAP"
    return(out_obj)
  }
}
