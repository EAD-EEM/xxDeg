##' Read the clipboard (or other file) into a list consisting of a data frame
##' and a header.
##'
##' Reads the contents of the clipboard (or another connection) into
##' a header and a data frame. It works only for two column, tab-separated
##' text, and assumes that the last non-numeric line in the text is the header.
##' @param x Either "clipboard" (the default) or the name of a file
##' @return A list containing two elements \code{head}, the first line of text
##' before the table, and \code{data}, a data frame containing the tabular data.
##' @export
clipread <- function(x="clipboard") {
  if(x=="clipboard")
    clip <- readClipboard()
  else
    clip <- readLines(x)

  s <- strsplit(clip, "\\t")
  ## Break out the data, which is assumed to be the set of lines which have
  ## the maximum number of columns
  columns <- sapply(s,length)
  maxcolumns <- max(columns)
  D <- do.call(rbind, s[columns == maxcolumns])
  header <- D[1,c(1,maxcolumns)] # save this for the redefinition of D below
  if(grepl("French",Sys.getlocale())) {
     D <- apply(D, 2,gsub, pattern=',',replacement='.')
  }
  C1 <- suppressWarnings(as.numeric(D[,1]))
  C2 <- suppressWarnings(as.numeric(D[,maxcolumns]))
  Numeric <- !is.na(C1) & !is.na(C2)
  ## Remove any lines, except the header, which are not numeric
  D[-1,] <- D[-1,][Numeric[-1]]
  ## Redefine D as a Data frame with header saved above in "header"
  D <- data.frame(C1[-1], C2[-1])
  names(D) <- header
  shortlines <- columns < maxcolumns
  if(any(shortlines)) title <- sub("\\s*$","",clip[shortlines][1], perl=TRUE)
  else title <- NA
  return(list(head=title, data=D))
}


##   ## this determines which lines have pairs of numbers. It should skip anything
##   ## else (such as NA values, blanks, three values, etc)
##   nums <- grep("^[0-9.eE\\-]+\\t[0-9.eE\\-]+$", clip, perl=TRUE)
##   if(length(nums) <= 4) {
##     warning("Insufficient data (fewer than five points) found")
##     return(NULL)
##   }
##   cols <- grep("^[^\\t]+\\t[^\\t]+$", clip, perl=TRUE)
##   data <- strsplit(clip[nums], "\\t", perl=TRUE)
##   data <- matrix(as.numeric(unlist(data)), nrow=length(nums), byrow=TRUE)
##   data <- as.data.frame(data)
##   if( (nums[1]-1) == cols[1])
##     names(data) <- strsplit(clip[cols[1]], "\\t", perl=TRUE)[[1]]
##   else names(data) <- c("Day","Conc")
##   if(cols[1] == 1) return(list(head=NA, data=data))
##   else {
##     head <- sub("\\s$","",clip[1:(cols[1]-1)], perl=TRUE)
##     return(list(head=head, data=data))
##   }
## }

##' Convert a 2-D character array, such as obtained from RExcel,  to a vector
##' of strings and a data frame
##'
##' @param a The character array
##' @return A 2-element list with a vector of strings and a data frame. Only
##' the first and last columns of a are returned
##' @export
arrayRead <- function(a) {
  ## similar to clipread. Reads a 2-D character array. Takes the first line
  ## before numerical lines to be the header and takes the stuff before the
  ## header as additional information to be printed above the plots

  ## remove lines with nothing in first column, and remove middle columns
  a <- a[length(a[,1])>0,c(1, dim(a)[2])]

  ## find last line with any letter or number other than "E" (which could be
  ## exponential notation) or a "-" sign after a digit (which could be a
  ## submission number)
  header.line <- which.max(grep("([A-DF-Za-df-z])|([0-9]-)",a[,1]))
  data <- as.data.frame(apply(a[(header.line+1):dim(a)[1],],2,as.numeric))
#  data <- lapply(data, as.numeric)
  if(any(length(a[header.line,])==0)) { # last text line is not header
                                        ## so make up a header line
    names(data) <- c("Time","Concentration") # default values
    head <- a[1:header.line, 1]
  }
  else {
    names(data) <- a[header.line,] # might make sure this has length 2
    head <- a[1:(header.line-1),1]
  }
  list(head=head, data=data)
}

##' This function converts parameters from IORE to FOMC
##' @param C0 A vector of parameters, or C0. If C0, then n and k must also be
##' provided. If a vector, must be of length 3 and have names "C0", "n" and "k".
##' or be in the order (C0, n, k).
##' @param n A value for the IORE n parameter, if only C0 is passed for C0.
##' @param k A value for the IORE k parameter, if only C0 is passed for C0.
##' @return The converted parameters, in a vector with names C0, alpha, beta.
##' @export
iore2fomc <- function(C0, n=NULL, k=NULL) {
   if(is.null(n) | is.null(k)) { # Assume all parameters are in C0
      if(inherits(C0, "nls") | inherits(C0, "opiorefit")) {
         p <- coef(C0) # if passed the whole model, get coefs
      } else if(length(C0) == 3) {
         if(is.null(names(C0))) names(C0) <- c("c0","n","k")
         p <- C0
      } else stop("Invalid values passed to iore2fomc")
   } else  p <- c(c0=C0, n=n, k=k)
   ret <- c(p['c0'], 1/(p['n']-1), p['c0']^(1-p['n'])/(p['k']*(p['n']-1)))
   names(ret) <- c("c0", "alpha","beta")
   ret
}

##' This function converts parameters from FOMC to IORE
##' @param C0 A vector of parameters, or C0. If C0, then alpha and beta must
##' also be provided. If a vector, must be of length 3 and have names
##' "C0", "alpha" and "beta", or be in the order (C0, alpha, beta).
##' @param alpha  A value for the FOMC alpha parameter, if only C0 is passed
##' as C0.
##' @param beta  A value for the FOMC beta parameter, if only C0 is passed
##' as C0.
##' @return The converted parameters.
##' @export
fomc2iore <- function (C0, alpha=NULL, beta=NULL) {
   if(is.null(alpha) | is.null(beta)) { # Assume all parameters are in C0
      if(inherits(C0, "nls")) { # Assume it is an FOMC fit (unlikely to be used)
         p <- coef(C0) # if passed the whole model, get coefs
      } else if(length(C0) == 3) {
         if(is.null(names(C0))) names(C0) <- c("c0","alpha","beta")
         p <- C0
      } else stop("Invalid values passed to iore2fomc")
   } else  p <- c(c0=C0, alpha=alpha,beta=beta)
   n <- 1 + 1/p['alpha']
   k <- p['c0']^(1-n)/(p['beta']*(n-1))
   ret <- c(p['c0'], n, k)
   names(ret) <- c("c0","n","k")
   ret
}

##' Apply the chi squared test as described by FOCUS
##'
##' @param time a vector of observation times
##' @param obs a vector of observed values
##' @param calc a vector of calculated values
##' @param nparm the number of parameters in the fit used for the
##' calulated values
##' @param level the level of the test, defaults to 0.05
##' @return a vector containing the error level and tabulated chi squared value
##' @export
focus.chisq <- function(time, obs, calc, nparm, level=0.05) {
  require(stats)
  av.obs <- mean(obs, na.rm=TRUE)
  obs <- tapply(obs, time, mean, na.rm=TRUE)
  if(length(obs) <= nparm) return(c(NA, NA))
  calc <- tapply(calc, time, mean , na.rm=TRUE)
  ssresid <- sum((calc - obs)^2)
  chi2table <- qchisq(1-level, length(obs)-nparm)
  error.level <- 100 * sqrt(1/chi2table * ssresid/av.obs^2)
  c(error.level,chi2table)
}
##' This function runs the summary function on all fits in fits and returns
##' data frame made by concatenating the parameters or coefficients elements
##' returned by R's summary() function
##'
##' @param fits A list of fits as returned by degfit()
##' @param file A file name. If blank the results will be returned and printed
##' only. If set to "clipboard" the results will be written to the clipboard.
##' @return A matrix with rows for each parameter in all the fits and
##' columns as returned by R's summary() functions
##' @export
parameterstats <- function(fits, file="") {
   R <- data.frame()
   for(f in fits) {
      if(!inherits(f, c("lm","nls","FitEmpty"))) next
      if(inherits(f, "lm")) {
         S <- summary(f)$coefficients
         row.names(S) <- c("log c0","k")
      }
      else if(inherits(f,c("nls", "FitEmpty"))) {
         S <- summary(f)$parameters
      }
      R <- rbind(R, S)
   }
   if(file != "") write.table(R, file, sep="\t")
   R
}

##' This function sets up a class for a failed curve fit. It is a list with
##' those elements needed by other functions, but with all values set to NA
##'
##' @param P A character vector with the names of the parameters in the model
##' @param N The number of parameters in the model that failed.
##' @return A list with elements chisq and dtx, which are NA vectors of
##' length 2,  spfits, which is set to NA, and par, which is a named vector
##' using the  parameter names passed in P, but with all values set to NA
BadFit <- function(P, N=length(P)) {
   x <- list(chisq=c(NA, NA), dtx=c(dt50=NA, dt90=NA), spfits=NA,
             par=c(rep(NA, length(P))), N=N)
   names(x[['par']]) <- P
   class(x) <- "FitEmpty"
   x
}

##' A summary method for the FitEmpty class.
##'
##' @param object A FitEmpty object, such as created by the BadFit function.
##' @param ... Other arguments that might be passed to a summary function.
##' Ignored.
##' @return A list with one element, a data frame named "parameters". The values
##' of the parameters vector are all NA, and the names depend on the number of
##' parameters for the FitEmpty object
summary.FitEmpty <- function(object, ...) {
   parnames <- names(object[['par']])
   print(parnames)
   S <- matrix(rep(NA, 4*length(parnames)), nc=4)
   colnames(S) <- c("Estimate",  "Std. Error", "t value", "Pr(>|t|)")
   print(S)
   rownames(S) <- parnames
       ## if(object$N == 2) return(c(c0=NA, k=NA))
       ## else if(object$N == 3) return(c(c0=NA, n=NA, k=NA))
       ## else if(object$N==4) return(c(c0=NA, r0=NA, c1=NA, r1=NA))
       ## else NULL
   S
}

##' A coef function for the FitEmpty class.
##' @param object A FitEmpty object such as created by the BadFit function
##' @param ... Other arguments that might be passed to a coef function. Ignored.
##' @return The par element of the FitEmpty object
coef.FitEmpty <- function(object, ...) return(object$par)

##' A predict function for the FitEmpty class.
##'
##' @param object A FitEmpty object, not used, present for compatibility
##' @param newdata A data frame. The number of rows us used to determine the
##' length of the return vector
##' @param ... Any other parameters that might be passed to "predict".
##' @return A vector of NA values with length the same as the number of rows in
##' \code{newdata}
predict.FitEmpty <- function(object, newdata, ...) return(rep(NA,nrow(newdata)))

##' A residuals function for the FitEmpty class. Returns NA values
##'
##' @param object A FitEmpty object, which would be created by a failed DFOP
##' or IORE fit.
##' @param ... Any other parameters which could be passed to a predict function.
##' Ignored
residuals.FitEmpty <- function(object, ...) return(rep(NA, object$N))



##' This is a function for generating a density plot of all possible GUS
##' scores for a set of DT50 and KOC values
##' @param dt50 A vector of DT50 values, in days
##' @param koc A vector of Koc values, in L/kg
##' @return The vector of GUS scores used to generate the plot, i.e. all
##' possible GUS scores.
##' @export
GUS <- function(dt50, koc) {
  dt50 <- dt50[!is.na(dt50)]
  koc <- koc[!is.na(koc)]
  gus <- as.vector(outer(dt50,koc, function(d,k) log(d,10)*(4-log(k,10))))
  ## The x-axis must span at least the range from 1.4 to 3.2 since the GUS
  ## cutoffs are at 1.8 and 2.8. Do this the easy way by making a dummy plot
  den <- density(gus)
  X <- range(c(1.4, 3.2, den[[1]]))

  par(mar=c(3.5,0.5,0.5,0.5), mgp=c(2,.25,0))
  plot(den, xaxt='n', yaxt='n', main="", xlab="GUS score", ylab="",
       xlim=c(X[1],X[2]))
  axis(1,c((par('usr')[1]+1.8)/2, 2.3, (par('usr')[2]+2.8)/2),
       c("Non-leacher","Borderline","Leacher"),tick=FALSE)
  axis(1,c(1.8,2.8),cex.axis=0.7)
  abline(v=c(1.8,2.8), lty=2)
  rug(gus)
  return(gus)
}
##' Runs the GUS function from data on the clipboard
##'
##' A convenience function for reading data off the clipboard and
##' using it to call the GUS function. This will assuming column
##' headings starting with numbers are numeric values which should
##' be included in the calculation.
##' @return An vector of all possible GUS scores
##' @export
GUSplot <- function() {
   decimal <- ifelse(grepl("French",Sys.getlocale()), ',','.')
   d <- read.delim("clipboard", dec=decimal)
   if(all(substr(colnames(d),1,1) == 'X'))
       d <- read.delim("clipboard", header=FALSE, dec=decimal)
   if(ncol(d) != 2)
       print("GUSplot must be called with exactly two columns, for halflife and Koc, on the clipboard")
   else
       return(GUS(na.omit(d[,1]), na.omit(d[,2])))
}


