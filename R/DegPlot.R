##' Sets up a Grid viewport with appropriate spacing for axis labels
##' @param maxX The maximum x-value
##' @param maxY The maximum y-value
##' @param vpname A name for the viewport
##' @param MainPlt If TRUE, the viewport is for the main plot. This is used for
##' setting spacing for axis titles on the plots
##' @param xax Whether or not to include x-axis labels, default TRUE
##' @return The created viewport
deg.vp <- function(maxX, maxY, vpname, MainPlt=TRUE, xax=TRUE) {
   require(grid)
   if(length(maxX)>1) maxX <- max(maxX)
   if(length(maxY)>1) maxY <- max(maxY)
   yax <- if(MainPlt) 1 else 0
   viewport(name=vpname, x=unit(yax*3,"lines"), y=unit(xax*3,"lines"),
            width=unit(1,"npc") - unit(3,"lines"),
            height = unit(1,"npc")-unit(xax*3,"lines"),just=c("left","bottom"),
            xscale=c(0,maxX) + c(-0.05, 0.05) * maxX,
            #Note the yscale varies between main and residual plots:
            yscale=maxY*c(MainPlt-1, 1) + c(-0.05, 0.05) * maxY*2*(1-MainPlt))
}

##' A function for making a grob for the main plot, which consists of one set
##' of data point and any number of fits.
##' @param x The x (time) values for the data
##' @param y The y (concentration) data values
##' @param x2, the x (time) values for the fits. Fit data must all be on the
##' same time values
##' @param y2 A matrix with fitted values for each fit in the columns.
##' @param clr A vector of colors for the plotted lines
##' @param name A name for the returned gTree
##' @param draw Boolean. Whether to draw the plot. Ignored.
##' @param gp An object of class gpar; typically the output from a call to the
##' function gpar.
##' @param vp A viewport to use for the returned gTree
##' @return A gTree for the main plot
mainPlot <- function(x,y, x2, y2, clr=c("cyan","magenta","green2"),
                     name=NULL, draw=FALSE, gp=gpar(), vp=NULL) {
   require(grid)
   maxX <- max(c(x,x2), na.rm=TRUE)
   maxY <- max(c(y, as.vector(y2)), na.rm=TRUE) # y2 is a matrix
   fitnames <- colnames(y2)
   ## determine location for legend based on Y values for the last 1/3 of data
   Xc <- max(x, na.rm=TRUE)*0.7
   endY <- max(c(y[x>Xc], as.vector(y2[x2>Xc,])), na.rm=TRUE)
   Loc <- endY < 0.7*maxY
   spg <- gTree(x=x, y=y, name=name,
                childrenvp=deg.vp(maxX, maxY, 'dataregion'),
                children=gList(rectGrob(name="border", vp="dataregion",
                                         gp=gpar(lwd=1.6)),
                xaxisGrob(name="xaxis", vp="dataregion",gp=gpar(fontsize=9, lwd=0.1)),
                yaxisGrob(name="yaxis", vp="dataregion",gp=gpar(fontsize=9)),
                pointsGrob(x, y, name="points", pch=19, size=unit(0.6,"char"),
                           vp="dataregion"),
                textGrob("Time", y=unit(-3,"lines"), name="xlab",
                         gp=gpar(fontsize=10), vp="dataregion"),
                textGrob("Concentration", x=unit(-3,"lines"), name="ylab",
                         gp=gpar(fontsize=10), rot=90, vp="dataregion"),
                degLegend(toupper(fitnames), clr, 1:3, upper=Loc)),
                gp=gp, vp=vp, cl="splot")
   for(i in 1:ncol(y2)) {
      spg <- addGrob(spg,  linesGrob(x2, y2[,i], name=fitnames[i],
                                     vp="dataregion",default.units="native",
                                     gp=gpar(col=clr[i], lty=i, lwd=2)))
   }
   spg
}

##' Generate a legend for the degradation plot
##'
##' @param names The names of the fits
##' @param colors The colors to use for the lines
##' @param linetypes The linetypes to use for the lines
##' @param upper If TRUE (the default) place the legend in the upper right
##' If FALSE, place it in the lower left
##' @param box If TRUE draw a box around the legend
##' @param gp An object of class gpar; typically the output from a call to the
##' function gpar.
##' @return A grob for the legend
degLegend <- function(names, colors, linetypes, upper=TRUE, box=TRUE,
                      gp=gpar()) {
   require(grid)
   if(upper) {
      J <- c("right","top")     # J is the "just" parameter for the viewport
      os <- -2                  # os is the offset from the parent edge
   } else {
      J <- c("left","bottom")
      os <- 2
   }
   upper <- as.integer(upper) # use as numeric  value for location
   LC <- unit(3.6*(1-upper),"lines") # compensation for 3-line axes in lower pos
   leg <- gTree(childrenvp=viewport(x=unit(upper,"npc") + LC + unit(os,"mm"),
                y=unit(upper,"npc") + LC + unit(os,"mm"), width=unit(25,"mm"),
                height=unit(3,"lines"),just=J, name="legendvp"),
                children=gList(rectGrob(vp="legendvp")),
                gp=gpar(fontsize=10))
   thelayout <- grid.layout(3,3, heights=unit(rep(1,3),"lines"),
                            widths=unit(c(11, 2, 1),"mm"))
   fg <- frameGrob(layout=thelayout, vp="legendvp")
   for(i in 1:length(names)) {
      fg <- placeGrob(fg, row=i, col=1, grob = linesGrob(y=unit(0.5,"npc"),
                             gp=gpar(col=colors[i],lty=linetypes[i],lwd=2)))
      fg <- packGrob(fg, row=i, col=3,grob=textGrob(names[i],x=0,just=c(0,0.4)))
   }
   leg <- addGrob(leg, fg)
   leg
}

##' A function for drawing the residual plot. The plot has the y-axis on the
##' right, and only optionally has an x-axis.
##' @param x1 The x values for the residuals
##' @param y1 The y values for the residuals
##' @param ymax The absolute value of the y scale. The plot is done with
##' zero in the middle.
##' @param vpname A name for the viewport in which to draw the plot
##' @param xaxis Whether to draw an x-axis
##' @param splKnots The number of knots to use in the natural spline curve
##' @param name A name for the gTree returned by the function
##' @param lgp A set of gpar parameters for the spline
##' @param gp A set of gpar parameters for all other aspects of the plot
##' @param vp An optional viewport. Currently unused
##' @return a gTree
residPlot <- function(x1,y1, ymax, vpname, xaxis=FALSE,
                      splKnots=3, name=NULL, lgp=gpar(), gp=gpar(), vp=NULL) {
   require(grid)
   ## make sure number of unique x values is more than knots
   doSpline <- splKnots >= 3 & sum(!is.na(unique(x1)))> splKnots+1
   if(doSpline) SplCrv <- myspline(x1, y1, k=splKnots[1])
   rpg <- gTree(x1=x1, y1=y1, name=name,
                childrenvp=deg.vp(x1,ymax,vpname=vpname,MainPlt=FALSE,
                                 xax=xaxis),
                children=gList(rectGrob(name="border", vp=vpname,
                                        gp=gpar(lwd=1.6)),
                yaxisGrob(name="yaxis", vp=vpname, main=FALSE,
                          gp=gpar(fontsize=7)),
                pointsGrob(x1, y1, pch=19, size=unit(0.5,"char"),
                           name="points", vp=vpname),
                linesGrob(y=unit(0.5,"npc"),vp=vpname,name="center")),
                gp=gp, vp=vp, cl="rplot")
   if(doSpline ) {
      rpg <- addGrob(rpg, linesGrob(SplCrv$xSpline, SplCrv$ySpline, gp=lgp,
                               name="lines", vp=vpname, default.units="native"))
   }
   if(xaxis) {
      rpg <- addGrob(rpg,
                     xaxisGrob(name="xaxis", vp=vpname, gp=gpar(fontsize=9)))
   }
   if(all(is.na(y1))) {
      rpg <- addGrob(rpg, textGrob("(Failed)", y=unit(.75, "npc"),hjust=1.2,
                                   gp=gpar(fontsize=8)))
   }
   rpg
}
##' A function which fits a natural spline with k knots to data
##' @param x X values of the data
##' @param y Y values of the data
##' @param k Number of knots to use for the natural spline
##' @return A data frame with 121 (x, y) pairs
myspline <- function(x, y, k=3) {
   require(splines)
   xSpline <- seq(0, max(x), length.out=121)
   Spl <- try(lm(y~ns(x, df=k)))
   if(inherits(Spl, "try-error")) return(NA)
   ySpline <- predict(Spl, data.frame(x=xSpline))
   data.frame(xSpline, ySpline)
}

##' A function for plotting the results of a list of curve fits
##'
##' Creates a new plotting window if none is available
##' @param fits A list containing a data frame with the original data and a
##' each of the fits to plot.
##' @param Title A title for the plot. A default title of "this is there the
##' title goes" will be used if no title is supplied
##' @return the plot
##' @export
degplot <- function(fits, Title=NULL, draw=TRUE) {
   require(grid)
   if(dev.cur() < 2) x11(width=6.5, height=4) # 1/2 8.5x11 page
   if(is.null(Title)) Title <- attr(fits[["Data"]], "title") # from Data attr
   if(is.null(Title)) Title <- "This is where the title goes" #default value
   tbl <- dataTableGrob(fits[["Data"]])
   DF <- fits[["Data"]]   # Separate the Data element from the rest of fits
   fits[["Data"]] <- NULL # Data used, no longer needed, so remove it
   stattbl <- summaryTable(summarizeFits(fits, level=0.5))
   title=textGrob(Title, x=0, just=c(0,0),gp=gpar(cex=1.5, lineheight=0.5))

   x2 <- seq(0, max(DF[,1]), length.out=121)
   # the sfo (and loglmsfo) models return a log of the data, so the predicted
   # values have to be exponentiated.
   fitCrv <- sapply(fits, function(fit, newdata) {
      x <- predict(fit, newdata)
      if(inherits(fit, "lm")) x <- exp(x)
      x
      }, newdata=data.frame(day=x2))

   rdl <- sapply(fits, residuals) #get data frame of residuals
   resRange <- max(abs(rdl),na.rm=TRUE)

   clr <- c("cyan","magenta","green2")
   gplot <- mainPlot(DF[,1], DF[,2], x2, fitCrv, clr, "main")
   rplot1 <- residPlot(DF[,1], rdl[,1], resRange, "resid1",
                       lgp=gpar(col=clr[1], lty=1), splKnots=4)
   rplot2 <- residPlot(DF[,1], rdl[,2], resRange, "resid2",
                       lgp=gpar(col=clr[2], lty=2))
   rplot3 <- residPlot(DF[,1], rdl[,3], resRange, "resid3", xaxis=TRUE,
                       lgp=gpar(col=clr[3], lty=3),splKnots=4)
   lay1 <- grid.layout(6,3,widths=unit(c(1,2,1), c("grobwidth","null","null"),
                           list(tbl, NULL, NULL)),
                      heights=unit(c(2,1,1,1,3,1),
                      c("grobheight","null","null","null","lines","grobheight"),
                      list(title, NULL, NULL, NULL,NULL, stattbl)))

   grid.newpage()
   fg <- frameGrob(layout=lay1)
   fg <- placeGrob(fg, title, row=1, col=1:3)
   fg <- placeGrob(fg, tbl, col=1, row=2:6)
   fg <- placeGrob(fg, gplot, col=2, row=2:5)
   fg <- placeGrob(fg, rplot1, col=3, row=2)
   fg <- placeGrob(fg, rplot2, col=3, row=3)
   fg <- placeGrob(fg, rplot3, col=3, row=4:5)
   fg <- placeGrob(fg, stattbl, col=2:3, row=6)
   if(draw) {
      grid.draw(fg)
   }
   return(fg)
}

##' Extracts summary values from a list of fits
##' @param fits A list of curve fits. If an element named "Data" is present
##' it will be ignored
##' @param level The level for calculating the Sc comparison value.
##' Default is 0.5
##' @return a list containing DT50, DT90, Chi squared and fitted parameters.
##' Each element is a vector, except for the fitted parameters which are
##' returned as a list with a named vector for each curve fit. Also returns
##' an elemen "Chosen" which designates which model should be used under
##' EAD-EFED guidance
##' @export
summarizeFits <- function(fits, level=0.5)  {
   fits["Data"] <- NULL # We don't use the data component of fits
   models <- names(fits)
   dtx <- sapply(fits, "[[", "dtx")
   chi <- sapply(fits, function(x) x[["chisq"]][1])
   bic <- sapply(fits, BIC)
   ioreN <- coef(fits[["iore"]])[2]
   ## Get the coefficients, converting any fits done with an lm type model
   co <- lapply(fits, function(x) {
      y <- coef(x)
      if(inherits(x, "lm")) {
         y[1] <- exp(y[1])
         y[2] <- -y[2]
         names(y) <- c("c0","r")
      }
      y
   })
   S <- sapply(fits, function(x) sum(residuals(x)^2))
   if("iore" %in% models) {
      npar <- 3
      N <- length(residuals(fits[["iore"]]))
      F <- qf(1-level, npar, N-npar)
      Sc <- S["iore"]*(1+npar/(N-npar)*F)
   }
   else Sc <- NA
   ## redefine the coefficients for the dfop curve to use fraction, f, and
   ## make rates positive.
   if("dfop" %in% models) {
      co[["dfop"]] <- c(c0 = co[["dfop"]]["c0"] + co[["dfop"]]["c1"],
                     f = co[["dfop"]]["c0"]/sum(co[["dfop"]][c("c0","c1")]),
                     -co[["dfop"]][c("r0","r1")])
      names(co[["dfop"]]) <- c("c0","f","r0","r1")
      dfopf <- co[["dfop"]]["f"]
      slowThalf <- log(2)/min(-coef(fits[["dfop"]])[c("r0","r1")])
   }
   else slowThalf <- NA
   ## Calculate the EFED-EAD method chosen value:
   if(ioreN < 1) chosen <- "sfo"
   else if(is.na(Sc) | is.na(S['sfo'])) chosen <- "none"
   else if(S['sfo'] < Sc) chosen <- "sfo"
   else if(is.na(slowThalf) & is.na(dtx['dt90', 'iore'])) chosen <- "none"
   else if(is.na(slowThalf) | slowThalf <= 0 | dfopf>1 | dfopf<0)
       chosen <- "iore"
   else if(is.na((dtx['dt90', 'iore']))) chosen <- "dfop"
   else if(slowThalf < dtx['dt90','iore'] * log(2)/log(10)) chosen <- "dfop"
   else chosen <- "iore"
   list(Models=models, DT50=dtx[1,], DT90=dtx[2,], ChiSq=chi, P=co,
        slowThalf=slowThalf, S=S, Sc=Sc, Chosen=chosen)
}

##' Generates a grob for a summary table, including DT50, DT90 and fitted
##' parameters
##' @param S The output from summarizeFits
##' @return a grob representing the table
summaryTable <- function(S) {
   require(grid)
   ## Set up some standard gpar values for rectangles, headings and body
   hgp <- gpar(fill="gray90", col="white")
   bgp <- gpar(fill="gray96", col="white")
   htgp <- gpar(fontsize=10, fontface="bold") # header text
   btgp <- gpar(fontsize=9) # body text
   dt50 <- lapply(formatC(S$DT50, digits=3, format="fg"), textGrob, gp=btgp)
   dt90 <- lapply(formatC(S$DT90, digits=3, format="fg"), textGrob, gp=btgp)
   chi <- lapply(formatC(S$ChiSq, digits=2, format="fg"), textGrob, gp=btgp)

   c0 <- sapply(S$P, "[", "c0")
###   names(c0) <- sub("\\.c0","",names(c0))
   c0 <- lapply(formatC(c0, digits=3, format="fg"), textGrob, gp=btgp)
   ## Write the remaining coefficients as a list of equations
   co <- lapply(S$P, formatC, digits=3, format="g")
   co <- list(substitute(k==X, list(X=co[["sfo"]][2])),
              substitute(list(f==X, k[0]==Y, k[1]==Z),
                  list(X=co[["dfop"]][2],Y=co[["dfop"]][3],Z=co[["dfop"]][4])),
              substitute(list(N==X, k==Y),
                  list(X=co[["iore"]][2], Y=co[["iore"]][3])))
   co <- lapply(co, textGrob, gp=btgp)
   headings <- list(expression(bold(DT[50])), expression(bold(DT[90])),
                   expression(bold(chi^2)), expression(bold(C[0])),"Parameters")
   headings <- lapply(headings, textGrob, gp=htgp)
   mNames <- lapply(toupper(S$Models), textGrob, gp=htgp)

   ## Values for the right hand column of the right hand table
   R <- formatC(c(S$Sc, S$S["sfo"], S$slowThalf, S$DT90['iore']*log(2)/log(10)),
                digits=3, format="g")
   RG <- lapply(R, textGrob, gp=btgp)
   ## highlight the EFED-EAD chosen value:
   highlight <- gpar(col="red2")
   if(S$Chosen == "sfo")
       dt50$sfo <- editGrob(dt50$sfo, gp=highlight)
   else if(S$Chosen == "dfop")
       RG[[3]] <- editGrob(RG[[3]], gp=highlight)
   else if(S$Chosen == "iore")
       RG[[4]] <- editGrob(RG[[4]], gp=highlight)
   ## Otherwise don't highlight anything because NA values preclude the method

   ## Start assembling the frameGrob, into which we will put the above values.
   ## Note that it is divided into two tables, the first with six columns and
   ## the second with two columns, with a blank column between them. There
   ## are also small buffer columns around each data column
   Widths <- unit(c(0.2,10,rep(c(0.4, 1, 0.1),5),0.2,10,rep(c(0.6, 1, 0.6),2)),"mm")
   thelayout <- grid.layout(4,25, heights=unit(c(1.1, 1, 1, 1), "lines"),
                         widths=Widths)
   fg <- frameGrob(name="summaryTable", gp=gpar(fill="white"), layout=thelayout)
   for(i in 1:5) {   # Header of left (main) table
      fg <- placeGrob(fg, grid.rect(gp=hgp), row=1, col=(i*3+1):(i*3+3))
      fg <- packGrob(fg, headings[[i]], row=1, col=3*i+2)
   }
   for(i in 2:4) {   # Body of left (main) table
      fg <- placeGrob(fg, grid.rect(gp=hgp), row=i, col=1:3)
      fg <- packGrob(fg, mNames[[i-1]], row=i, col=2)
      fg <- placeGrob(fg, grid.rect(gp=bgp), row=i, col=4:6)
      fg <- packGrob(fg, dt50[[i-1]], row=i, col=5)
      fg <- placeGrob(fg, grid.rect(gp=bgp), row=i, col=7:9)
      fg <- packGrob(fg, dt90[[i-1]], row=i, col=8)
      fg <- placeGrob(fg, grid.rect(gp=bgp), row=i, col=10:12)
      fg <- packGrob(fg, chi[[i-1]], row=i, col=11)
      fg <- placeGrob(fg, grid.rect(gp=bgp), row=i, col=13:15)
      fg <- packGrob(fg, c0[[i-1]],  row=i, col=14)
      fg <- placeGrob(fg, grid.rect(gp=bgp), row=i, col=16:18)
      fg <- packGrob(fg, co[[i-1]], row=i, col=17)
   }
   ## Right hand table. Has no header and two columns, label (left) column
   fg <- placeGrob(fg, grid.rect(gp=hgp), row=1, col=20:22)
   fg <- packGrob(fg, textGrob(expression(bold(S[C])), gp=htgp),row=1,col=21)
   fg <- placeGrob(fg, grid.rect(gp=hgp), row=2, col=20:22)
   fg <- packGrob(fg, textGrob(expression(bold(S[SFO])), gp=htgp),row=2,col=21)
   fg <- placeGrob(fg, grid.rect(gp=hgp), row=3, col=20:22)
   fg <- packGrob(fg, textGrob(expression(bold(Slow~~t[1/2])), gp=htgp),row=3,col=21)
   fg <- placeGrob(fg, grid.rect(gp=hgp), row=4, col=20:22)
   fg <- packGrob(fg, textGrob(expression(bold(t[R~~IORE])), gp=htgp),row=4,col=21)
   ## Right (value) column. The RG variable is set just before assembly of fg
   for(i in 1:4){
      fg <- placeGrob(fg, grid.rect(gp=bgp), row=i, col=23:25)
      fg <- packGrob(fg, RG[[i]], row=i,col=24)
   }

   fg
}

##' Generates a grob of a table listing a two column data frame
##'
##' @param DF A data frame. Only the first two columns used
##' @param gp An object of class gpar; typically the output from a call to the
##' function gpar.
##' @return A frameGrob with four cells. All data are squeezed into single
##' cells using newlines.
dataTableGrob <- function(DF, gp=gpar()) {
   require(grid)
   hgp <- gpar(fill="gray90", col="white")
   bgp <- gpar(fill="gray96", col="white")
   H1 <- textGrob(names(DF)[1], gp=gpar(fontface="bold", fontsize=10))
   H2 <- textGrob(names(DF)[2], gp=gpar(fontface="bold", fontsize=10))
   N1 <- paste(formatC(DF[,1], digits=3,format="fg"), collapse="\n")
   N2 <- paste(formatC(DF[,2], digits=3,format="fg",flag="#"), collapse="\n")
   B1 <- textGrob(N1, gp=gpar(lineheight=0.8, fontsize=10))
   B2 <- textGrob(N2, gp=gpar(lineheight=0.8, fontsize=10))
   thelayout <- grid.layout(3, 6,     heights=unit(c(1,1,0.1),"lines"),
                            widths=unit(c(.2,1,.6,.2,1,.2), "mm"))
   fg <- frameGrob(name="datatable", gp=gpar(fill="white"),layout=thelayout)
   fg <- placeGrob(fg, grid.rect(gp=hgp), row=1, col=1:3)
   fg <- packGrob(fg, H1, row=1, col=2)
   fg <- placeGrob(fg, grid.rect(gp=hgp), row=1, col=4:6)
   fg <- packGrob(fg, H2, row=1, col=5)
   fg <- placeGrob(fg, grid.rect(gp=bgp), row=2:3, col=1:3)
   fg <- packGrob(fg, B1, row=2, col=2)
   fg <- placeGrob(fg, grid.rect(gp=bgp), row=2:3, col=4:6)
   fg <- packGrob(fg, B2, row=2, col=5)
   fg
}

##' Convert the results of \code{summarizeFits} to a data frame
##'
##' @param S the results of a call to \code{summarizeFits}. If S is a list
##' with a first element named "Data", then it will be taken as a set of fits
##' (as generated by \code{degfit} and \code{summarizeFits} will be run on it.
##' @param writeto The name of a file or connection to write the results to
##' @return A data frame in the same giving the DT50, DT90 ChiSquared statistic
##' and parameters for each fit. The parameters of each model are different
##' but are lumped into the last four columns of the data frame with column
##' labels such as "r/f/n" meaning \code{r} for the first (sfo) fit, \code{f}
##' for the second and \code{n} for the third.
##' @export
fitsummary2DF <- function(S, writeto=NULL) {
   if(names(S)[1] == "Data")
       S <- summarizeFits(S)
   DF <- data.frame(S[c("DT50","DT90","ChiSq")])
   ## Add in a matrix of parameters. This lists c0 in the first column, but
   ## after that the parameters vary for each model. All parameters are stored
   ## as lists S$P and so must be converted to a matrix before adding to DF
   params <- t(sapply(S$P, function(x) c(x, rep(NA,4-length(x)))))
   DF <- cbind(DF, params)
   names(DF)[5:ncol(DF)] <- c("r/f/n","-/r0/k","-/r1/-")
   DF
}

##' Generate a 1 row data frame with the values from the rightmost column of
##' table in the summary plot
##'
##' @param S the results of a call to \code{summarizeFits}. If S is a list
##' with a first element named "Data", then it will be taken as a set of fits
##' (as generated by \code{degfit} and \code{summarizeFits} will be run on it.
##' @param writeto The name of a file or connection to write the results to
##' @return A data frame with a single line, containing values for Sc, Ssfo,
##' the slow half-life from the DFOP fit, the half-life of an exponential
##' curve passing through the DT90 of the IORE fit, and the model chosen by the
##' NAFTA project agreed method.
##' @export
fitsummary.other <- function(S, writeto="NULL") {
   if(names(S)[1] == "Data")
       S <- summarizeFits(S)
   DF <- data.frame(Sc=S$Sc, Ssfo=S[['S']]['sfo'], dfopslow=S$slowThalf,
              TrIORE=S[["DT50"]]["iore"], chosenmodel=toupper(S$Chosen),
              row.names='value')
   if(!is.null(writeto))
       write.table(DF, writeto, sep="\t", row.names=FALSE,quote=FALSE)
   DF
}

## DF <- data.frame(Day=rep(c(0,1,3,7,14,28,60,90,120,150,180),each=2),
##                  Conc = rep(seq(100,8, length.out=11),each=2)+rnorm(22))
## DF[,2] <- round(DF[,2], 4-ceiling(log(max(DF[,2]),10))) # round to cnst figs
## fits <- degfit(DF, "This is my title")

##' Make a ggplot2 plot of a set of fits, similar to the main plot
##' created by degPlot
##'
##' @param fits A named list of curve fits, as returned by degfit
##' @return a ggplot2 plot of the data and (usually) three curve fits
##' @export
gdegplot <- function(fits, nrow=NULL) {
   ## fits should be a list of results from the degfit function,
   ## and names should be a vector of names the elements of fits
   x <- mapply(function(x,y) cbind(Soil=x, y$Data, Model="Data"),
               names(fits), fits, SIMPLIFY=FALSE)
   d <- do.call(rbind, x)
   names(d) <- c("Soil","Day","Conc","Model")
#   pred <- mapply(evalall, fits, names, tmax=max(d[,1]))

   # generate the curve fits, first as a list and then rbind together
   predL <- mapply(degeval, fits, names(fits),
               MoreArgs = list(tmax=max(d[,2]), n=10), SIMPLIFY=FALSE)
   pred <- do.call(rbind, predL)
   P <- qplot(Day, Conc, data=d, colour=Model, linetype=Model)
   P <- P + geom_line(data=pred, size=I(1.4))
   P <- P + scale_colour_manual(values=c("black","cyan", "magenta","green3"))
   P <- P+scale_linetype_manual(values=c(4,1,2,3))
   if(is.null(nrow))
      thefacets <- facet_wrap(~Soil)
   else
      thefacets <- facet_wrap(~Soil, nrow=nrow)
   P  + thefacets +  ylim(0,103) + theme(legend.position="top")
}

##' Make a dataframe for plotting a set of fits, similar to the main
##' plot created by degPlot. Only works with three (sfo, dfop iore)
##' fits
##'
##' @param fit A curve fit, as returned by degfit
##' @param name The name of the soil (or other identifier) for the fit
##' @param tmax The maximum time in the output data frame
##' @param n the number of points to generate
##' @return a data frame with Soil, Model, Day and Conc columns
degeval <- function(fit, name, tmax, n=101) {
   trange <- data.frame(day = seq(0, tmax, length.out=n))
   sfo <- predict(fit$sfo, newdata=trange)
   dfop <- predict(fit$dfop, newdata=trange)
   iore <- predict(fit$iore, newdata=trange)
   data.frame(Soil=name, Model=rep(c("SFO","DFOP","IORE"), each=n),
              Day=rep(trange$day,3),  Conc = c(exp(sfo), dfop, iore))

}
