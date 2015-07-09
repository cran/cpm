\name{getStatistics}
     \alias{getStatistics}
     \title{Returns the Test Statistics Associated With A CPM S4 Object}
     \author{Gordon J. Ross \email{gordon@gordonjross.co.uk}}
     \description{
       Returns the \eqn{D_{k,t}} statistics associated with an existing Change Point Model (CPM) S4 object. These statistics depend on the state of the object, which depends on the observations which have been processed to date. Calling this function returns the most recent set of statistics, which were generated after the previous observation was processed. 
    
	Note that this function is part of the S4 object section of the \code{cpm} package, which allows for more precise control over the change detection process. For many simple change detection applications this extra complexity will not be required, and the \code{detectChangePoint} and \code{processStream} functions should be used instead.
    
    For a fuller overview of this function including a description of the CPM framework and examples of how to use the various functions, please consult the package manual "Parametric and Nonparametric Sequential Change Detection in R: The cpm Package" available from www.gordonjross.co.uk
    }
    
    \usage{
     getStatistics(cpm)
     }
     \arguments{
       \item{cpm}{The CPM S4 object for which the test statistics are to be returned.}
     }
     \value{
       A vector containing the \eqn{D_{k,t}} statistics generated after the previous observation was processed.
     }
     \seealso{
       \code{\link{makeChangePointModel}, \link{processObservation}, \link{changeDetected}}.
     }
     \examples{
#generate a sequence containing two change points
x <- c(rnorm(200,0,1),rnorm(200,1,1),rnorm(200,0,1))

#vectors to hold the result
detectiontimes <- numeric()
changepoints <- numeric()

#use a Lepage CPM
cpm <- makeChangePointModel(cpmType="Lepage", ARL0=500)


i <- 0
while (i < length(x)) {
  i <- i + 1
  #process each observation in turn
  cpm <- processObservation(cpm,x[i])
  
  #if a change has been found, log it, and reset the CPM
  if (changeDetected(cpm) == TRUE) {
    print(sprintf("Change detected at observation \%d", i))
    detectiontimes <- c(detectiontimes,i)
    
    #the change point estimate is the maximum D_kt statistic 
    Ds <- getStatistics(cpm)
    tau <- which.max(Ds)
    
    if (length(changepoints) > 0) {
      tau <- tau + changepoints[length(changepoints)]
    }
    changepoints <- c(changepoints,tau)
    
    #reset the CPM
    cpm <- cpmReset(cpm)
    
    #resume monitoring from the observation following the 
    #change point
    i <- tau
  }
}
}