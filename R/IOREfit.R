##' @include xxDeg-utils.R
NULL

##' A function for fitting an nth order rate equation to data using \code{nls}
##' @param day the independent variable, time of observarion
##' @param conc the dependent variable, the concentration at each time
##' @param start an initial guess of parameter values. Parameters are the initial
##'   concentration, n, the order of the equation and k the rate
##' @param gradient If TRUE asks the iore function to calculate a gradient.
##' Default FALSE
##' @param ... Other parameters passed to \code{nls}
##' @return an R nls object, with added properties chisq, giving the chisq
##'   error level, eq, a plotmath expression for the fitted equation and
##'   parameters, and dtx, a vector contining the dt50 and dt90 of the fitted
##'   curve.
##' @export
fit.iore <- function(day, conc, start=NULL, gradient=FALSE, ...) {
   ## start, if entered, should be a 3-vector with c0, n and k in that order
   ## if start is not entered, then make a guess either by fitting
   ## log(dC/dt)=log(k) + n*log(C) (only when mean(C) at each time is
   ## monotonically decreasing) or by assuming n=2.5 and estimating k from that
   ## and the endpoints of the curve
   ISS <- function(P, day, conc) { # IORE sum of squares; P=c(c0, n, k)
      if(P[2] ==1) R <- sum((P[1]*exp(-P[3]*day)-conc)^2)
      else {
         pred <- P[1]*(P[1]^(P[2]-1)*(P[2]-1)*P[3]*day +1)^(1/(1-P[2]))
         R <- sum((pred-conc)^2)
      }
      return(R)
   }
   ## Sort the data for the starting values estimate:
   ord <- order(day)
   dsort <- day[ord]
   csort <- conc[ord]
   d1 <- unique(dsort)
   c1 <- as.vector(tapply(csort, dsort, mean, na.rm=TRUE))
   if(is.null(start)) {
      if(all(sign(diff(c1))==-1)) { # monotonically decreasing average conc.
         ## fit the log ODE: log(C') = log(k) + N log(C),
         ## where y=log(C') and x is the midpoint between times
         y <- log(-diff(c1)/diff(d1))
         x <- log((c1[-1]+c1[-length(c1)])/2)
         defit <- lm(y~x)
         start <- c(max(c1), coef(defit)[2], exp(coef(defit)[1]))
      }
      else {
         n <- 2.15
         k <- (min(c1)^(1-n) - max(c1)^(1-n))/((n-1) * d1[length(d1)])
         start <- c(max(c1), n, k)
      }
   }
   ## the values may be meaningless, but we'll make sure the names are right
   names(start) <- c("c0", "n", "k")
   fit.op <- try(optim(start, ISS, day=day, conc=conc))
   if(!inherits(fit.op,"try-error")) {     # Use params from the optim fit, if
     start=c(fit.op$par[1], fit.op$par[2], fit.op$par[3])  # possible
     names(start) <- c("c0", "n", "k")
  }
   ## Now use the nls function to do the same fit
   ## The "port" algorithm seems more robust than the default of plinear
   print(start)
   # start <- log(start)  # Not sure why this was here (Apr 2011)
   ##fit.nls <- try(nls(conc~c0*(c0^(n-1)*(n-1)*k*day +1)^(1/(1-n)),start=start,
   fit.nls <- try(nls(conc~iore(day, c0, n, k, gradient), start=start,
                   algorithm="port"))
   ## If the fit fails, try other values of n:
   if(inherits(fit.nls, "try-error") | inherits(fit.op,"try-error")) {
      cn <- min(c1)
      c0 <- c1[1]
      ns <- exp(seq(log(.5), log(30), length.out=20))
      ks <- ((c0/cn)^(ns-1)-1)/(c0^(ns-1)*(ns-1) * d1[length(d1)])
      ## Calculate the sums of squares of all these
      ss <- mapply(function(a, b) ISS(c(100,a,b),day, conc), ns,ks)
      i <- which.min(ss)
      start1 <- c(c0, ns[i], ks[i])
      names(start1) <- c("c0","n","k")
      print("trying starting parameters")
      print(start1)
      if(inherits(fit.nls,"try-error"))
      ##    fit.nls <- try(nls(conc~c0*(c0^(n-1)*(n-1)*k*day +1)^(1/(1-n)),
          fit.nls <- try(nls(conc~iore(day, c0, n, k, gradient),
                           start=start1, algorithm="port"))
      if(inherits(fit.op, "try-error")) # make another optim try if needed
        fit.op <- optim(start1, ISS, day=day, conc=conc)
   }

   if(inherits(fit.nls, "try-error") & inherits(fit.op,"try-error"))
       return(BadFit(c("c0","n","k"),length(day)))
   else if(inherits(fit.nls,"nls")) { # use the nls fit first if it exists
      co <- coef(fit.nls)
      pred <- predict(fit.nls, newdata=data.frame(day=day))
      iore <- fit.nls
   }
   else { # the optim fit is the only one that worked
          ## Run nls with zero iterations, to get the fit into an nls object
      co <- fit.op$par
      print(names(co))
      print(co)
      iore <- nls(conc~iore_g(day, c0, n, k), start=co,
                  control=list(maxiter=0, warnOnly=TRUE),algorithm="plinear")
      pred <- predict(iore, newdata=data.frame(day=day))
   }
   n1 <- co[2] - 1 # n-1
   dt90 <- (10^n1-1)/(n1*co[3]*co[1]^n1)
   dt50 <- (2^n1-1)/(n1*co[3]*co[1]^n1)
   iore[["dtx"]] <- c(dt50, dt90)
   names(iore[["dtx"]]) <- c("dt50","dt90") # using = above doesn't work
   iore[["chisq"]] <- focus.chisq(day, conc, pred, 3)
   iore[['eq']] <- substitute(list(frac(1,C^{n-1}) == frac(1,C[0]^{n-1}) + (n-1)*k*t,~~where~~C[0]==c0, n==N, k==K), list(c0 = signif(co[1],2), N=signif(co[2],2), K=signif(co[3],2)))
    ## substitute(C(t) == frac(c0, bgroup("[",frac(t,bt)+1,"]")^al), list(c0=signif(co[[1]],3), al=-signif(co[[2]],3), bt=signif(1/co[[3]],3)))
    return(iore)
}

##' A coef function for for the opiorefit class, which is returned by the
##' fit.iore function in cases where the nls function fails to fit the iore
##' model, but the optim function succeeds
##' @param fit An object of class opiorefit, which is a list returned by the
##' optim function, with added elements dtx, chisq and eq
##' @param ... Other arguments that might be passed to a coef function. Ignored.
##' @return The coefficients of the fitted IORE model, i.e. the  $par
##' element from the list returned by optim.
coef.opiorefit <- function(fit, ...) return(fit[['par']])

##' A residuals function for for the opiorefit class, which is returned
##' by the fit.iore function in cases where the nls function fails to
##' fit the iore model, but the optim function succeeds.
##' @param fit An object of class opiorefit, which is a list returned by the
##' optim function, with added elements dtx, chisq, residuals and eq
##' @param ... Other parameters that might be passed to a residuals function.
##' Ignored.
##' @return The coefficients of the fitted IORE model, i.e.
##' the $residuals element added to the list returned by optim
residuals.opiorefit <- function(fit, ...) return(fit[['residuals']])

##' A predict function for for the opiorefit class, which is returned
##' by the fit.iore function in cases where the nls function fails to
##' fit the iore model, but the optim function succeeds.
##' @param fit An object of class opiorefit, which is a list returned by the
##' optim function, with added elements dtx, chisq, residuals and eq
##' @param newdata An optional data frame in which to look for
##' variables with which to predict. If omitted, the fitted values are used.
##' @param ... Other parameters that might be passed to a predict function.
##' Ignored.
##' @return A vector of predictions
predict.opiorefit <- function(fit, newdata=NULL, ...) {
   if(is.null(newdata)) return(fit[['pred']])
   co <- fit$par
   co[1]*(co[1]^(co[2]-1)*(co[2]-1)*co[3]*newdata$day +1)^(1/(1-co[2]))
}

##' This is a method logLik method for class opiorefit. It is basically
##' copied from logLik.nls, except that it ignores weights (since there
##' will be no weights).
##' @param object The opiorefit object
##' @param REML Must be FALSE (the default)
##' @param ... Other paramters used by other LogLik functions. Ignored.
##' @return The log Likelihood for the model
logLik.opiorefit <- function(object, REML=FALSE, ...) {
   if (REML)
       stop("cannot calculate REML log-likelihood for \"nls\" objects")
   res <- residuals(object)
   N <- length(res)
   val <- -N * (log(2*pi)+1 - log(N) + log(sum(res^2))) /2.0
   attr(val, "df") <- 1L + length(coef(object))
   attr(val, "nobs") <- attr(val, "nall") <- N
   class(val) <- "logLik"
   val
}

##' Calculates the IORE function and gradient for use by nls
##'
##' @param T The time value for IORE
##' @param c0 The starting concentration parameter
##' @param n The order of the reaction
##' @param k The rate constant for the IORE equation
##' @param gradient Whether to return the gradient in an attribute ("gradient")
##' @return The value of the IORE function at each time in T, together with
##' a gradient attribute with the gradient along each parameter
##' @export
iore <- function(T, c0, n, k, gradient=FALSE){
   model.func <- ( c0^(1-n)- (1-n)*k*T)^(1/(1-n))
   if(gradient) {
      a <- c0^(1-n) - k*(1-n)*T
      Z <- cbind(a^(n/(1-n))/c0^n,   # d/dc0
              -((a^(1/(1-n))) * (k*(1-n)*T*c0*n + (k*(n-1)*T*c0^n+c0)*
                               log(a)+c0*(n-1)*log(c0))) /
              ((n-1)^2 * (k*(n-1)*T*c0^n + c0)),
              T*(a)^(n/(1-n))) #d/dk
      dimnames(Z) <- list(NULL, c("c0","n","k"))
      attr(model.func, "gradient") <- Z
   }
   model.func
}

##' Calculates the IORE function and gradient for use by nls. Always includes
##' The gradient
##'
##' @param T The time value for IORE
##' @param c0 The starting concentration parameter
##' @param n The order of the reaction
##' @param k The rate constant for the IORE equation
##' @return The value of the IORE function at each time in T, together with
##' a gradient attribute with the gradient along each parameter
##' @export
iore_g<- function(T, c0, n, k){
   model.func <- ( c0^(1-n)- (1-n)*k*T)^(1/(1-n))
   a <- c0^(1-n) - k*(1-n)*T
   Z <- cbind(a^(n/(1-n))/c0^n,   # d/dc0
           -((a^(1/(1-n))) * (k*(1-n)*T*c0*n + (k*(n-1)*T*c0^n+c0)*
                            log(a)+c0*(n-1)*log(c0))) /
           ((n-1)^2 * (k*(n-1)*T*c0^n + c0)),
           T*(a)^(n/(1-n))) #d/dk
   dimnames(Z) <- list(NULL, c("c0","n","k"))
   attr(model.func, "gradient") <- Z

   model.func
}
