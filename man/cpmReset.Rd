\name{cpmReset}
     \alias{cpmReset}
     \title{Resets a CPM S4 Object}
     \author{Gordon J. Ross \email{gordon@gordonjross.co.uk}}
     \description{
       Resets an existing Change Point Model (CPM) S4 object. Typically it will be used in situations where the CPM is processing a stream, and has encountered a change point. If the stream may contain multiple change points, the typical reaction procedure after detecting a change is to reset the CPM. This involves clearing its memory of all the observations prior to the change point, so that monitoring can be resumed starting with the next observation after the change point.
    
        Note that this function is part of the S4 object section of the \code{cpm} package, which allows for more precise control over the change detection process. For many simple change detection applications this extra complexity will not be required, and the \code{detectChangePoint} and \code{processStream} functions should be used instead.
        
	For a fuller overview of this function including a description of the CPM framework and examples of how to use the various functions, please consult the package manual "Parametric and Nonparametric Sequential Change Detection in R: The cpm Package" available from www.gordonjross.co.uk
    }
    
    \usage{
     cpmReset(cpm)
     }
     \arguments{
       \item{cpm}{The CPM S4 object which is to be reset.}
     }
     \value{
       The reinitialised CPM.
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