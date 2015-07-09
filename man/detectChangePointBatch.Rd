\name{detectChangePointBatch}
     \alias{detectChangePointBatch}
     \title{Detect a Single Change Point in a Sequence}
     \author{Gordon J. Ross \email{gordon@gordonjross.co.uk}}
     \description{
       This function is used to detect a single change point in a sequence of observations using the Change Point Model (CPM) framework for batch (Phase I) change detection. The observations are processed in one batch
       and information is returned regarding whether the sequence contains a change point. A full description of the CPM framework can be found in the papers cited in the reference section.
    
	For a fuller overview of this function including a description of the CPM framework and examples of how to use the various functions, please consult the package manual "Parametric and Nonparametric Sequential Change Detection in R: The cpm Package" available from www.gordonjross.co.uk
    }
    
     \usage{
     detectChangePointBatch(x, cpmType, alpha=0.05, lambda=NA) 
     }
     \arguments{
       \item{x}{A vector containing the univariate data stream to be processed.}
       
       \item{cpmType}{The type of CPM which is used. With the exception of the FET, these CPMs are all implemented in their two sided forms,
and are able to detect both increases and decreases in the parameters monitored. Possible arguments are: 
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
    	\item{alpha}{the null hypothesis of no change is rejected if \eqn{D_n > h_n} where n is the length of the sequence and \eqn{h_n} is the upper alpha percentile of the test statistic distribution. Because computing the values of \eqn{h_n} associated with each value of alpha is a laborious task, the package includes a function \code{getBatchThreshold} which returns the threshold associated with a particular choice of alpha and n. This function is called automatically whenever \code{detectChangePointBatch} is called, so the user need only provide the desired value of alpha. The allowable values for this argument are 0.05, 0.01, 0.005, 0.001. If a different value is required then the user will need to compute it manually.}
       \item{lambda}{A smoothing parameter which is used to reduce the discreteness of the test statistic when using the FET CPM. See [Ross and Adams, 2012b] in the References section for more details on how this parameter is used.  Currently the package only contains sequences of ARL0 thresholds corresponding to lambda=0.1 and lambda=0.3, so using other values will result in an error. If no value is specified, the default value will be 0.1.}
     }
     \value{
       \item{x}{The sequence of observations which was processed.}
       \item{changeDetected}{TRUE if \eqn{D_n} exceeds the value of \eqn{h_n} associated with alpha, otherwise FALSE.}
       \item{changePoint}{assuming a change was detected, this stores the most likely location of the change point, defined as the value of k which maximized \eqn{D_kt}. If no change is detected, this will be equal to 0.}
       \item{threshold}{The value of \eqn{h_n} which corresponds to the specified alpha.}
       \item{Ds}{The sequence of \eqn{D_kt} test statistics.}
     }
    
     \seealso{
       \code{\link{detectChangePoint}}.
     }
     \examples{
     ## Use a Student-t CPM to detect a mean shift in a stream of Gaussian 
     ## random variables which occurs after the 100th observation
     x <- c(rnorm(100,0,1),rnorm(1000,1,1))
     detectChangePointBatch(x,"Student",alpha=0.05)


     ## Use a Mood CPM to detect a scale shift in a stream of Student-t 
     ## random variables which occurs after the 100th observation
     x <- c(rt(100,5),rt(1000,5*2))
     detectChangePointBatch(x,"Mood",alpha=0.05)

     ## Use a FET CPM to detect a parameter shift in a stream of Bernoulli 
     ##observations. In this case, the lambda parameter acts to reduce the 
     ##discreteness of the test statistic.
     x <- c(rbinom(100,1,0.2), rbinom(1000,1,0.5))
     detectChangePointBatch(x,"FET",alpha=0.05,lambda=0.3)
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