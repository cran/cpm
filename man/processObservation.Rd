\name{processObservation}
     \alias{processObservation}
     \title{Process an Observation Using a CPM S4 Object.}
     \author{Gordon J. Ross \email{gordon@gordonjross.co.uk}}
     \description{
       Updates the state of an existing Change Point Model (CPM) S4 object, by processing a single observation. This effectively computes the \eqn{D_{k,t+1}} statistics, for a CPM that had previously seen \eqn{t} observations. When the function is called, several events happen. First, the function returns a CPM object which is identical the CPM object passed to the function, except that the observation passed as an argument has been processed and added to the state. Second, the CPM computes the \eqn{D_{t+1}} statistic and compares it to its stored sequence of thresholds. If a change is detected, then this is stored in the state of the CPM, and a call to \code{changeDetected} will now return \code{TRUE}.
       
           Note that this function is part of the S4 object section of the \code{cpm} package, which allows for more precise control over the change detection process. For many simple change detection applications this extra complexity will not be required, and the \code{detectChangePoint} and \code{processStream} functions should be used instead.
           
	For a fuller overview of this function including a description of the CPM framework and examples of how to use the various functions, please consult the package manual "Parametric and Nonparametric Sequential Change Detection in R: The cpm Package" available from www.gordonjross.co.uk
    }
    
    \usage{
     processObservation(cpm,x)
     }
     \arguments{
       \item{cpm}{The CPM S4 object which is to be updated.}
       \item{x}{The observation which is to be processed.}
     }
     \value{
     The updated CPM. If a stream is being processed, then this should be stored, and used to process the next observation in the sequence.
     }
     \seealso{
       \code{\link{makeChangePointModel}, \link{changeDetected}}.
     }
     \examples{
#generate a sequence containing a single change point
x <- c(rnorm(100,0,1),rnorm(100,1,1))

#use a Student CPM
cpm <- makeChangePointModel(cpmType="Student", ARL0=500)

for (i in 1:length(x)) {

  #process each observation in turn
  cpm <- processObservation(cpm,x[i])
  
  if (changeDetected(cpm)) {
    print(sprintf("change detected at observation \%s",i))
    break
  }
}
}

     