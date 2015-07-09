\name{getBatchThreshold}
     \alias{getBatchThreshold}
     \title{Returns the Threshold Associated with a Type I Error Probability .}
     \author{Gordon J. Ross \email{gordon@gordonjross.co.uk}}
     \description{
	
     When performing Phase I analysis within the CPM framework for a sequence of length n, the null hypothesis of no change is rejected if \eqn{D_n > h_n} for some threshold \eqn{h_n}. Typically this threshold is chosen to be the upper alpha quantile of the distribution of \eqn{D_n} under the null hypothesis of no change. Given a particular choice of alpha and n, this function returns the associated \eqn{h_n} threshold. Because these thresholds are laborious to compute, the package contains pre-computed values of \eqn{h_n} for alpha = 0.05, 0.01, 0.005 and 0.001, and for \eqn{n < 10000}.

	For a fuller overview of this function including a description of the CPM framework and examples of how to use the various functions, please consult the package manual "Parametric and Nonparametric Sequential Change Detection in R: The cpm Package" available from www.gordonjross.co.uk
    }
    
    \usage{
	getBatchThreshold(cpmType, alpha, n, lambda=0.3)
    }

\arguments{
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
\item{alpha}{the null hypothesis of no change is rejected if \eqn{D_n > h_n} where n is the length of the sequence and \eqn{h_n} is the upper alpha percentile of the test statistic distribution.} 
\item{n}{the sequence length the value should be calculated for, i.e. the value of n in \eqn{D_n}.} 
 \item{lambda}{A smoothing parameter which is used to reduce the discreteness of the test statistic when using the FET CPM. See [Ross and Adams, 2012b] in the References section for more details on how this parameter is used.  Currently the package only contains sequences of ARL0 thresholds corresponding to lambda=0.1 and lambda=0.3, so using other values will result in an error. If no value is specified, the default value will be 0.1.}
}


     \seealso{
       \code{\link{detectChangePointBatch}}.
     }


     \examples{
     ## Returns the threshold for n=1000, alpha=0.05 and the Mann-Whitney CPM
     h <- getBatchThreshold("Mann-Whitney", 0.05, 1000)
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



     