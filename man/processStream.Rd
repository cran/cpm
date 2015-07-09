\name{processStream}
     \alias{processStream}
     \title{Detects Multiple Change Points in a Sequence}
     \author{Gordon J. Ross \email{gordon@gordonjross.co.uk}}
     \description{
       This function is used to detect a multiple change points in a sequence of observations using the Change Point Model (CPM) framework for sequential (Phase II) change detection. The observations are processed in order, starting with the first, and a decision is made after each observation whether a change point has occurred. A full description of the CPM framework can be found in the papers cited in the reference section.
       
       Unlike the \code{detectChange} function, \code{processStream} does not terminate and return when a change point is encountered. Instead, a new CPM is initialised immediately following the change point, with all previous observations being discarded. The monitoring then continues, starting from the first observation after the change point. If more change points are discovered later in the sequence, the CPM is again reinitialised after each one. In this way, the whole sequence of observations will be processed and multiple change points may be detected.        
       
       For a fuller overview of this function including a description of the CPM framework and examples of how to use the various functions, please consult the package manual "Parametric and Nonparametric Sequential Change Detection in R: The cpm Package" available from www.gordonjross.co.uk
    }
     \usage{
     processStream(x, cpmType, ARL0=500, startup=20, lambda=NA)
     }
     \arguments{
       \item{x}{A vector containing the univariate data stream to be processed.}
       
       \item{cpmType}{The type of CPM which is used. Possible arguments are: 
         \itemize{
           \item \code{Student}: Student-t test statistic, as in [Hawkins et al, 2003]. Use to detect mean changes in a Gaussian sequence.
           \item \code{Bartlett}: Bartlett test statistic, as in [Hawkins and Zamba, 2005]. Use to detect variance changes in a Gaussian sequence.
           \item \code{GLR}: Generalized Likelihood Ratio test statistic, as in [Hawkins and Zamba, 2005b]. Use to detect both mean and variance changes in a Gaussian sequence.
        	\item \code{Exponential}: Generalized Likelihood Ratio test statistic for the Exponential distribution, as in [Ross, 2013]. Used to detect changes in the parameter of an Exponentially distributed sequence.
           \item \code{GLRAdjusted} and \code{ExponentialAdjusted}: Identical to the GLR and Exponential statistics, except with the finite-sample correction discussed in [Ross, 2013] which can lead to more powerful change detection.
           \item \code{FET}: Fishers Exact Test statistic, as in [Ross and Adams, 2012b]. Use to detect parameter changes in a Bernoulli sequence.
           \item \code{Mann-Whitney}: Mann-Whitney test statistic, as in [Ross et al, 2011]. Use to detect location shifts in a stream with a (possibly unknown) non-Gaussian distribution.
           \item \code{Mood}: Mood test statistic, as in [Ross et al, 2011]. Use to detect scale shifts in a stream with a (possibly unknown) non-Gaussian distribution.
           \item \code{Lepage}: Lepage test statistics in [Ross et al, 2011]. Use to detect location and/or shifts in a stream with a (possibly unknown) non-Gaussian distribution.
           \item \code{Kolmogorov-Smirnov}: Kolmogorov-Smirnov test statistic, as in [Ross et al 2012]. Use to detect arbitrary changes in a stream with a (possibly unknown) non-Gaussian distribution.
           \item \code{Cramer-von-Mises}: Cramer-von-Mises test statistic, as in [Ross et al 2012]. Use to detect arbitrary changes in a stream with a (possibly unknown) non-Gaussian distribution.
         }
    	}     
       \item{ARL0}{Determines the \eqn{ARL_0} which the CPM should have, which corresponds to the average number of observations before a false positive occurs, assuming that the sequence does not undergo a chang. Because the thresholds of the CPM are computationally expensive to estimate, the package contains pre-computed values of the thresholds corresponding to several common values of the \eqn{ARL_0}. This means that only certain  values for the \eqn{ARL_0} are allowed. Specifically, the \eqn{ARL_0} must have one of the following values: 370, 500, 600, 700, ..., 1000, 2000, 3000, ..., 10000, 20000, ..., 50000.}
       \item{startup}{The number of observations after which monitoring begins. No change points will be flagged during this startup period. This should be set to at least 20.}
         
       \item{lambda}{A smoothing parameter which is used to reduce the discreteness of the test statistic when using the FET CPM. See [Ross and Adams, 2012b] in the References section for more details on how this parameter is used.  Currently the package only contains sequences of ARL0 thresholds corresponding to lambda=0.1 and lambda=0.3, so using other values will result in an error. If no value is specified, the default value will be 0.1.}
     }
     \value{
       \item{x}{The sequence of observations which was processed.}       \item{detectionTimes}{A vector containing the points in the sequence at which changes were detected, defined as the first observation after which \eqn{D_t} exceeded the test threshold.}
       \item{changePoints}{A vector containing the best estimates of the change point locations, for each detecting change point. If a change is detected after the \eqn{t^{th}} observation, then the change estimate is the value of \eqn{k} which maximises \eqn{D_{k,t}}.}
     }
     
     \seealso{
       \code{\link{detectChangePoint}}.
     }
     \examples{
     ## Use a Student-t CPM to detect several mean shift in a stream of 
     ## Gaussian random variables
     x <- c(rnorm(100,0,1),rnorm(100,1,1), rnorm(100,0,1), rnorm(100,-1,1))
     result <- processStream(x,"Student",ARL0=500,startup=20)
     plot(x)
     for (i in 1:length(result$changePoints)) {
         abline(v=result$changePoints[i], lty=2)
     }

     ## Use a Mood CPM to detect several scale shifts in a stream of 
     ##Student-t random variables
     x <- c(rt(100,3),rt(100,3)*2, rt(100,3), rt(100,3)*2)
     result <- processStream(x,"Mood",ARL0=500,startup=20)
     plot(x)
     for (i in 1:length(result$changePoints)) {
         abline(v=result$changePoints[i], lty=2)
     }
     
     ## Use a FET CPM to detect several parameter shifts in a stream of 
     ## Bernoulli observations. In this case, the lambda parameter acts to 
     ## reduce the discreteness of the test statistic.
     x <- c(rbinom(300,1,0.1),rbinom(300,1,0.4), rbinom(300,1,0.7))
     result <- processStream(x,"FET",ARL0=500,startup=20,lambda=0.3)
     plot(x)
     for (i in 1:length(result$changePoints)) {
         abline(v=result$changePoints[i], lty=2)
     }
     }
     
     
\references{
  Hawkins, D. , Zamba, K. (2005) -- A Change-Point Model for a Shift in Variance, \emph{Journal of Quality Technology}, 37, 21-31

  Hawkins, D. , Zamba, K. (2005b) -- Statistical Process Control for Shifts in Mean or Variance Using a Changepoint Formulation, \emph{Technometrics}, 47(2), 164-173
  
  Hawkins, D., Qiu, P., Kang, C. (2003) -- The Changepoint Model for Statistical Process Control, \emph{Journal of Quality Technology}, 35, 355-366. 
    
  Ross, G. J., Tasoulis, D. K., Adams, N. M. (2011) -- A Nonparametric Change-Point Model for Streaming Data, \emph{Technometrics}, 53(4)
  
  Ross, G. J., Adams, N. M. (2012) -- Two Nonparametric Control Charts for Detecting Arbitary Distribution Changes, \emph{Journal of Quality Technology}, 44:102-116
  
  Ross, G. J., Adams, N. M. (2013) -- Sequential Monitoring of a Proportion, \emph{Computational Statistics}, 28(2)

  Ross, G. J., (2014) -- Sequential Change Detection in the Presence of Unknown Parameters, \emph{Statistics and Computing} 24:1017-1030
  
  Ross, G. J., (2015) -- Parametric and Nonparametric Sequential Change Detection in R: The cpm Package, \emph{Journal of Statistical Software}, forthcoming
}     