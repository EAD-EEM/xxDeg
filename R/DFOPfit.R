##' @include xxDeg-utils.R
NULL

##' ##' A function for fitting a dual exponential to data using \code{nls}
##'
##' @param day the independent variable, time of observarion
##' @param conc the dependent variable, the concentration at each time
##' @param start an initial guess of parameter values. Used only for method="nls"
##' @param ... Additional parameters passed to \code{nls}
##' @return an R nls object, with added properties chisq, giving the chisq
##'   error level, eq, a plotmath expression for the fitted equation and
##'   parameters, and dtx, a vector contining the dt50 and dt90 of the fitted
##'   curve.
##' @export
fit.dfop <- function(day, conc, start=NULL, ...) {
  if(is.null(start)) {  ## estimate starting parameters
     ## first fit a single exponential to the bounding box around the data:
     x <- range(day)
     y <- rev(range(conc))  # this is expected to be decreasing, so reverse it
     rexp <- (log(y[2]) - log(y[1]))/(x[1]-x[2])
     Cexp <- y[1]*exp(rexp*x[1])
     ## Next find the ratio of a sum of two exponentials with rates five times
     ## apart, that still fit the same box
     r <- c(sqrt(5), 1/sqrt(5))* rexp
     f <- (y[2]/Cexp - exp(-r[2]*x[2])) / (exp(-r[1]*x[2]) - exp(-r[2]*x[2]))
     start=c(c0=Cexp*f, r0=-r[1], c1=Cexp*(1-f), r1=-r[2])
  }
  else {
    ## the values may be meaningless, but we'll make sure the names are right
    names(start) <- c("c0","r0","c1","r1")
  }
  ## Fit the curve using the optim function.
  ## If the fit fails, adjust the rates by a factor and try again
  start1 <- start
  for(V in c(1:3))  {
     start1[c(2,4)] <- start[c(2,4)]*c(V, 1/V)
     OP <- try(optim(start1, function(P, time, conc)
                 sum((P[1]*exp(P[2]*time)+P[3]*exp(P[4]*time)-conc)^2),
                 time=day,conc= conc))
     if(!inherits(OP,"try-error")) break
  }
  start <- if(inherits(OP, "try-error")) start else OP$par
  dfop <- try(nls(conc~c0*exp(r0*day)+c1*exp(r1*day), start=start, ...))
  ## If it fails, try making the rates more different, one last time
  if(inherits(dfop, "try-error"))  {
    start=start*c(1,3,1,.2)
    dfop <- try(nls(conc~c0*exp(r0*day)+c1*exp(r1*day), start=start, ...))
  }
  ## Set the return value, and run nls again to make sure r0 and r1 are fast
  ## and slow rates, respectively.
  if(inherits(dfop,"try-error") & inherits(OP,"try-error"))
      dfop <- BadFit(c("c0","r0","c1","r1"),length(day))
  else { # One of the two methods succeeded
     if(!inherits(dfop, "try-error")) { # nls worked
        co <- coef(dfop)
        if(co[2]>co[4]) {#The first half should be the fast part or redo flipped
           names(co) <- names(co)[c(3,4,1,2)]
           dfop <- nls(conc~c0*exp(r0*day)+c1*exp(r1*day),
                       start=co[c(3,4,1,2)],...)
           co <- coef(dfop)
        }
     }
     else  { # nls failed; return the optim result, but convert it back
             # to nls by fitting with zero iterations
        if(OP$par[2] > OP$par[4])  {  # rerunning this with the swapped starting
                              ## parameters seems to give a better fit
           OP <- try(optim(OP$par[c(3,4,1,2)], function(P, time, conc)
                 sum((P[1]*exp(P[2]*time)+P[3]*exp(P[4]*time)-conc)^2),
                 time=day,conc= conc))
        }
        start <- OP$par
        names(start) <- c('c0','r0','c1','r1')
        dfop <- nls(conc~c0*exp(r0*day)+c1*exp(r1*day), start=start,
                    control=list(warnOnly=TRUE, maxiter=0),
                    algorithm="plinear",  ...) # plinear my be less likely to give singular gradient error
        co <- coef(dfop)
     }
     pred <- predict(dfop, newdata=data.frame(day=day))

     ## The following applies to either nls or optim fits
     ## handle the case where the slow rate is an increase by finding minimum
     dfop[["chisq"]] <- focus.chisq(day, conc, pred, 4)
     mx <- if(sum(sign(co[c(2,4)]))==0)log(co[1]/co[3])/(co[4]-co[2])
     else  -3/co[4]
     dt90 <- try(uniroot(function(x) co[1]*exp(co[2]*x)+co[3]*exp(co[4]*x)
                         - 0.1*(co[1]+co[3]), interval= c(-2/co[2], mx)))
     if(inherits(dt90, "try-error")) {
        dt90 <- list(root=NA)
        if(sum(sign(co[c(2,4)]))==0) ## find upper end of interval in which to
         ## find dt50. If one rate is positive, use the minimum of fitted curve:
            xmax <- (log(co[3]*co[4])+log(co[2]*co[2]))/(co[2]-co[4])
        else
            xmax <- 100 * (-log(2)/co[2]) ## otherwise use 100 * half-life
     }
     else xmax <- dt90$root          ## If dt90 found, use it for upper end.
     ## here start interval a bit below one fast half-life.
     dt50 <- try(uniroot(function(x) co[1]*exp(co[2]*x)+co[3]*exp(co[4]*x)
                         -0.5*(co[1]+co[3]), interval= c(-0.4/co[2], xmax)))
     if(inherits(dt50, "try-error")) dt50 <- list(root=NA)
     dfop[["dtx"]] <- c(dt50=dt50$root, dt90=dt90$root)
     co <- as.list(signif(coef(dfop),3))
     dfop[['eq']] <- substitute(C(t)==c0*e^{r0*t}+c1*e^{r1*t}, co)
  }
  dfop
}

##' A coef function for for the opdfopfit class, which is returned by the
##' fit.dfop function in cases where the nls function fails to fit the dfop
##' model, but the optim function succeeds
##' @param fit An object of class opdfopfit, which is a list returned by the
##' optim function, with added elements dtx, chisq and eq
##' @param ... Additional paramters. Ignored
##' @return The coefficients of the fitted DFOP model, i.e. the  $par
##' element from the list returned by optim.
##' @export
coef.opdfopfit <- function(fit, ...) return(fit[['par']])

##' A residuals function for for the opdfopfit class, which is returned
##' by the fit.dfop function in cases where the nls function fails to
##' fit the dfop model, but the optim function succeeds.
##' @param fit An object of class opdfopfit, which is a list returned by the
##' optim function, with added elements dtx, chisq, residuals and eq
##' @param ... Additional parameters. Ignored.
##' @return The coefficients of the fitted DFOP model, i.e.
##' the $residuals element added to the list returned by optim
##' @export
residuals.opdfopfit <- function(fit, ...) return(fit[['residuals']])

##' A predict function for for the opdfopfit class, which is returned
##' by the fit.dfop function in cases where the nls function fails to
##' fit the dfop model, but the optim function succeeds.
##' @param fit An object of class opdfopfit, which is a list returned by the
##' optim function, with added elements dtx, chisq, residuals and eq
##' @param newdata An optional data frame in which to look for
##' variables with which to predict. If omitted, the fitted values are used.
##' @param ... Additional parameters. Ignored.
##' @return A vector of predictions
##' @export
predict.opdfopfit <- function(fit, newdata=NULL, ...) {
   if(is.null(newdata)) return(fit[['pred']])
   co <- fit$par
   co[1]*exp(co[2]*newdata$day) + co[3]*exp(co[4]*newdata$day)
}
