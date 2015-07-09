\name{cpm-package}
\title{The Change Point Model Package}
\alias{cpm-package}
\author{Gordon J. Ross \email{gordon@gordonjross.co.uk}}
\docType{package}
\description{An implementation of several different change point models (CPMs) for performing both parametric and nonparametric change detection on univariate data streams.}
\details{The CPM framework is an approach to sequential change detection (also known as Phase II process monitoring) which allows standard statistical hypothesis tests to be deployed sequentially. The main two general purpose functions in the package are \code{detectChangePoint} and \code{processStream} for detecting single and multiple change points respectively. The remainder of the functions allow for more precise control over the change detection procedure. To cite this R package in a research paper, please use \code{citation('cpm')} to obtain the reference, and BibTeX entry.

Note: this package has a manual titled "Parametric and Nonparametric Sequential Change Detection in R: The cpm Package" available from www.gordonjross.co.uk, which contains a full description of all the functions and algorithms in the package, as well as detailed instructions on how to use it. 

If you would like to cite this package, the citation information is "G. J. Ross - Parametric and Nonparametric Sequential Change Detection in R: The cpm Package, Journal of Statistical Software, 2015, 66(3), 1-20"


A Brief CPM Overview

Given a sequence \eqn{X_1,...,X_n} of random variables, the CPM works by evaluating a two-sample test statistic at every possible split point. Let \eqn{D_{k,n}} be the value of the test statistic when the sequence is split into the two samples \eqn{\{X_1, X_2,..., X_k\}} and \eqn{\{X_{k+1}, X_{k+2} ,..., X_{n}\}},  and define \eqn{D_n} to be the maximum of these values. \eqn{D_n} is then compared to some threshold, with a change being detected if the threshold is exceeded. 

In the sequential context, the observations are processed one-by-one, with \eqn{D_t} being computed based on the first \emph{t} observations, \eqn{D_{t+1}} being computed based on the first \emph{t+1} observations, and so on. The change detection time is defined as the first value of \eqn{t} where the threshold is exceeded. Supposing this occurs at time \eqn{t=T}, then the best estimate of the location of the change point is the value of \eqn{k} which maximised \eqn{D_{k,T}}. Writing \eqn{\hat{\tau}} for this, we have that \eqn{\hat{\tau} \leq T}.

The thresholds are chosen so that there is a constant probability of a false positive occurring after each observation. This leads to control of the Average Run Length (\eqn{ARL_0}), defined as the expected number of observations received before a change is falsely detecting, assuming that no change has occurred.

The choice of test statistic in the CPM defines the class of changes which it is optimised towards detecting. This package implements CPMs using the following statistics. More details can be found in the references section:

\itemize{
           \item \emph{Student}: Student-t test statistic, as in [Hawkins et al, 2003]. Use to detect mean changes in a Gaussian sequence.
           \item \emph{Bartlett}: Bartlett test statistic, as in [Hawkins and Zamba, 2005]. Use to detect variance changes in a Gaussian sequence.
           \item \code{GLR}: Generalized Likelihood Ratio test statistic, as in [Hawkins and Zamba, 2005b]. Use to detect both mean and variance changes in a Gaussian sequence.
           \item \code{Exponential}: Generalized Likelihood Ratio test statistic for the Exponential distribution, as in [Ross, 2013]. Used to detect changes in the parameter of an Exponentially distributed sequence.
           \item \code{GLRAdjusted} and \code{ExponentialAdjusted}: Identical to the GLR and Exponential statistics, except with the finite-sample correction discussed in [Ross, 2013] which can lead to more powerful change detection.
           \item \emph{FET}: Fishers Exact Test statistic, as in [Ross and Adams, 2012b]. Use to detect parameter changes in a Bernoulli sequence.
           \item \emph{Mann-Whitney}: Mann-Whitney test statistic, as in [Ross et al, 2011]. Use to detect location shifts in a stream with a (possibly unknown) non-Gaussian distribution.
           \item \emph{Mood}: Mood test statistic, as in [Ross et al, 2011]. Use to detect scale shifts in a stream with a (possibly unknown) non-Gaussian distribution.
           \item \emph{Lepage}: Lepage test statistics in [Ross et al, 2011]. Use to detect location and/ort shifts in a stream with a (possibly unknown) non-Gaussian distribution.
           \item \emph{Kolmogorov-Smirnov}: Kolmogorov-Smirnov test statistic, as in [Ross et al 2012]. Use to detect arbitrary changes in a stream with a (possibly unknown) non-Gaussian distribution.
           \item \emph{Cramer-von-Mises}: Cramer-von-Mises test statistic, as in [Ross et al 2012]. Use to detect arbitrary changes in a stream with a (possibly unknown) non-Gaussian distribution.
}

For a fuller overview of the package which includes a description of the CPM framework and examples of how to use the various functions, please consult the full package manual titled  "Parametric and Nonparametric Sequential Change Detection in R: The cpm Package" 

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