\name{ForexData}
\docType{data}
\alias{ForexData}
\title{Foreign Exchange Data}
\description{
This dataset consists of a historical sequence of the exchange rates between the Swiss Franc (CHF) and the British Pound (GBP). The value of the exchange rate was recorded at three hour intervals running from October 21st 2002, to May 15th 2007. In total, 9273 observations were made.

For an example of how to use the \code{cpm} package to detect changes in this data set, please see the "Streams with multiple change points' section of the package vignette, available by typing: \code{vignette("cpm")}
    
}
\usage{data(ForexData)}
\format{A three column matrix. Each observation is one row. The first column is the data, the second is the time of day at which the observation was recorded (1400 = 14:00pm, and so on), and the third is observed exchange rate.}
\keyword{datasets}