% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare_forecast.R
\name{Cforecast}
\alias{Cforecast}
\title{Forecasting time series}
\usage{
Cforecast(
  object,
  h = ifelse(frequency(object) > 1, 2 * frequency(object), 10),
  level = c(80, 95),
  fan = FALSE,
  robust = FALSE,
  lambda = NULL,
  biasadj = FALSE,
  find.frequency = FALSE,
  allow.multiplicative.trend = FALSE,
  model = NULL,
  plot = TRUE,
  dates,
  ...
)
}
\arguments{
\item{object}{a time series object or a matrix containing the time series values.}

\item{h}{number of periods for forecasting.}

\item{model}{an object describing a time series model; e.g. of class \code{\link{CArima}}.}

\item{plot}{if \code{TRUE}, it will show all the plots in execution.}

\item{dates}{a vector or an array with the dates as \code{character} of each observation of the time series. This vector can be generated using the function \code{\link{dateSeq}}.}

\item{...}{additional arguments to be passed to \code{\link[forecast]{forecast}}.}
}
\value{
See the \code{\link[forecast]{forecast}} function in the forecast package. The additional objects returned are: 
   \item{forecasts}{The forecasts with their confidence intervals as a \code{data.frame}.}
   \item{realForecPlot}{A real + forecast plot.}
   \item{realFitPlot}{A real vs fitted values plot.}
}
\description{
Largely a wrapper for the \code{\link[forecast]{forecast}} function in the forecast package. It adds some usefull datas and graphs to validate and compare models.
}
\details{
See the \code{\link[forecast]{forecast}} function in the forecast package.
}
\examples{
#An easy example with AirPassengers
dates <- dateSeq(from = time(AirPassengers)[1], length.out = length(AirPassengers),
                 by = "month")
model <- CArima(as.vector(log(AirPassengers)), order = c(0, 1, 1), 
                 seasonal = list(order=c(0,1,1), period = 12))
model   
                
predModel <- Cforecast(as.vector(log(AirPassengers)), model = model, h = 24, 
                       dates = dates, plot = F) 
predModel
summary(predModel)
}
\seealso{
\code{\link[ggplot2]{ggplot}}
}
