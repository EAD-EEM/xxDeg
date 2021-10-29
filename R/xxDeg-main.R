##' @include DegPlot.R
##' @include DFOPfit.R
##' @include IOREfit.R
##' @include SFOfit.R
NULL

##' A function for fitting degradation functions to data.
##' Based on ead.plot, but does not include plotting, instead returning a list
##' consisting of the original data and objects for each fit.
##'
##' @param DF a data frame with time in the first column and concentrations in
##'   the second. DF can also be set to clipboard in which case the contents of
##'   the clipboard will be converted to a data frame using \code{clipread}
##' @param pagehead a vector of strings to print at the top of the output page.
##'   These will be taken from any non-tabular lines in the clipboard if
##'   \code{DF} is set to clipboard, overwriting this parameter. The first
##'   element of \code{pagehead} is added as the attribute "title" to the
##'   Data element (i.e. \code{DF}) of the returned list.
##' @param tofit the models to fit. Allowed values are \code{sfo},
##' \code{logomsfo},  \code{dfop}, \code{iore} and \code{fomc}.
##' @return A list including the data (in an element named "Data") and the
##' fitted models, which are standard R fit types (e.g. lm, nls) with some
##' added data.
##' @export
degfit <- function(DF="clipboard",pagehead=NULL, tofit=c("sfo","dfop","iore")) {
   ## 1. Get the data
   if(inherits(DF, "character")) {
      input <- clipread(DF)
      DF <- input[['data']]
      pagehead <- input[['header']]
      if(is.null(DF)) return("No data found")
   }
   DF <- DF[!is.na(DF[,1]) & !is.na(DF[,2]), ] # remove lines with NA values
   attr(DF, "title") <- pagehead[1]
   ## Make sure the models are the ones allowed:
   mymodels <- c("sfo","dfop", "iore", "fomc", "loglmsfo")
   models <- as.list(mymodels[mymodels %in% tofit])
   nfit <- length(models)
   if(length(grep("^[0-9.eE\\-]$",names(DF)))>0)
       warning("Be sure first line of data is the column headings")
   timeunits <- strsplit(names(DF)[1], "[ ,(\\[)\\]]+", perl=TRUE )[[1]][2]
   if(is.na(timeunits)) timeunits <- ""
   ## 2. Do the fitting
   fits <- lapply(models, function(x)
                  eval(call(paste("fit",x,sep="."),DF[,1], DF[,2])))
   names(fits) <- models
   c(list(Data=DF), fits) # return a list, with DF as a single element
}

##' A convenience function to read in data from the clipboard or a file, fit
##' it and plot the results.
##' @param file The name of a file that can be read into a data frame. May have
##' descriptives lines of text before the (tab-delimited) data. If file is a
##' data frame, the first and last columns in that data frame will be used.
##' @param title An optional title for the plot. Default of NULL allows the
##' title to be taken from file, or set to a default if there is none in file.
##' @return a lst of fits, as returned by degFit
##' @export
##'

fitandplot <- function(file="clipboard", title=NULL) {
   if(is.character(file)) {
      DF <- clipread(file) # this will also cover file == "clipboard"
      D <- DF[[2]]
      title <- if(is.null(title)) DF[[1]]
   }
   else if(is.data.frame(file)) {
      D <- file  # title needs no change in this case
   }
   else stop("No valid data found for fitandplot\n")
   f <- degfit(D, title)
   degplot(f)
   f
}


