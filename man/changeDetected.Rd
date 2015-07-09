\name{changeDetected}
     \alias{changeDetected}
     \title{Tests Whether a CPM S4 Object Has Encountered a Change Point}
     \author{Gordon J. Ross \email{gordon@gordonjross.co.uk}}
     \description{
       Tests whether an existing Change Point Model (CPM) S4 object has encountered a change point. It returns \code{TRUE} if a change has been encountered, otherwise \code{FALSE}.
       
           Note that this function is part of the S4 object section of the \code{cpm} package, which allows for more precise control over the change detection process. For many simple change detection applications this extra complexity will not be required, and the \code{detectChangePoint} and \code{processStream} functions should be used instead.
           
	For a fuller overview of this function including a description of the CPM framework and examples of how to use the various functions, please consult the package manual "Parametric and Nonparametric Sequential Change Detection in R: The cpm Package" available from www.gordonjross.co.uk
    }
    
    \usage{
     changeDetected(cpm)
     }
     \arguments{
       \item{cpm}{The CPM S4 object which is to be tested for whether a change has occurred.}
     }
     
     \value{
       \code{TRUE} if a change has been detected, otherwise \code{FALSE}.
     }
     \seealso{
       \code{\link{makeChangePointModel}, \link{processObservation}}.
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

