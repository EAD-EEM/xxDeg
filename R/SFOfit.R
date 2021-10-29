##' @include xxDeg-utils.R
NULL

##' Calculate curve fits to pesticide degradation study data
##'
##' Fit an exponential (SFO) curve to data
##'
##' @param day the independent variable, time of observarion
##' @param conc the dependent variable, the concentration at each time
##' @param start an initial guess of parameter values. Used only for method="nls"
##' @param method any value other than "glm" will make R use the nls function
##' @param x logical value passed to the glm function
##' @return an R glm object, with added properties chisq, giving the chisq
##'   error level, eq, a plotmath expression for the fitted equation and
##'   parameters, and dtx, a vector contining the dt50 and dt90 of the fitted
##'   curve.
##' @export
fit.sfo <- function(day, conc, start=NULL, method="glm", x=TRUE) {
  ## x is set to TRUE so it can be used to calc. untransformed residuals
  day <- day[!is.na(conc)]
  conc <- conc[!is.na(conc)]
  if(is.null(start)) {
    p <- which.min(abs(conc-0.25*max(conc)))
    start <- c(c0=max(conc), r=-log(conc[p]/100)/day[p])
  }
  if(method == "glm")  sfo <- try(glm(conc~day, family=gaussian(link=log),x=x))
  else sfo <- try(nls(conc~c0*exp(-r*day),start=start, x=x))
  if(inherits(sfo, "try-error")) return(BadFit(c("(Intercept)","day"),length(day)))
  else {
    co <- coef(sfo)
    if(method=="glm") co[1] <- exp(co[1]) #glm(gaussian(link=log)) gives log(c0)
    names(co) <- c("c0","r")
    pred <- co[1]*exp(co[2]*day)
    sfo[["chisq"]] <- focus.chisq(day, conc, pred, 2)
    sfo[["dtx"]] <- -c(log(2)/co[2], log(10)/co[2])
    names(sfo[["dtx"]]) <- c("dt50", "dt90")
    sfo[['eq']] <- substitute(C(t)==c0*e^{r*t}, as.list(signif(co,3)))
  }
  sfo
}

##' A function which fits an exponential curve to data by taking the log of the
##' dependent variable and using linear regression
##'
##' @param day The independent variable, time of observarion
##' @param conc The dependent variable, the concentration at each time
##' @return an R lm object, with added properties chisq, giving the chisq
##'   error level, eq, a plotmath expression for the fitted equation and
##'   parameters, and dtx, a vector contining the dt50 and dt90 of the fitted
##'   curve.
##' @export
fit.loglmsfo <- function(day, conc) {
  fit <- lm(log(conc)~day)
  co <- coef(fit)
  names(co) <- c("c0","r")
  co[1] <- exp(co[1])
  pred <- co[1]*exp(co[2]*day)
  fit[['chisq']] <- focus.chisq(day, conc, pred, 2)
  fit[['eq']] <- substitute(C(t)==c0*e^{r*t}, as.list(signif(co,3)))
  fit[['dtx']] <- -log(c(dt50=2,dt90=10))/co[2]
  fit
}
