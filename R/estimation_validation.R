#' @title Fit ARIMA model to univariate time series
#'
#' @description Largely a wrapper for the \code{\link[forecast]{Arima}} function in the forecast package. It adds some usefull statistics and tests to validate and compare models.
#' @usage CArima(
#' y,
#' order = c(0,0,0),
#' seasonal = list(order = c(0,0,0), period = 1),
#' plot = T
#' xreg = NULL,
#' include.constant = TRUE,
#' lambda = model$lambda,
#' biasadj = FALSE,
#' method = c("CSS-ML", "ML", "CSS"),
#' model = NULL,
#' x = y,
#' ...
#' )
#' @param y a univariate time series of class \code{ts} or a matrix with the values of the time series.
#' @param order a specification of the non-seasonal part of the ARIMA model: the three integer components \code{(p, d, q)} are the AR order, the order of differencing, and the MA order.
#' @param seasonal a specification of the seasonal part of the ARIMA model, plus the period (which defaults to frequency(x) if the \code{y} parameter is a \code{ts} object). This may be a list with components order and period, or just a numeric vector of length 3 which specifies the seasonal order. In the latter case the default period is used.
#' @param include.constant for undifferenced series it fits the mean of the time series, for differenced series it fits the mean of the differenced time series. Note that if there is more than one difference taken, no constant is included regardless of the value of this argument.
#' @param plot if \code{TRUE}, it will show all the plots in execution.
#' @param ... additional arguments to be passed to \code{\link[forecast]{Arima}}.
#' @details See the \code{\link[forecast]{Arima}} function in the forecast package.
#' @return See the \code{\link[forecast]{Arima}} function in the forecast package. The additional objects returned are: 
#'    \item{period}{The time series period of the seasonal part.}
#'    \item{SBC}{The SBC value corresponding to the log-likelihood.}
#'    \item{cor.coef}{Correlation matrix of the parameters of the model.}
#'    \item{resid.acf}{ACF of the residuals obtained with the \code{\link[forecast]{Acf}} function of the forecast package.}
#'    \item{resid.pacf}{PACF of the residuals obtained with the \code{\link[forecast]{Pacf}} function of the forecast package.}
#'    \item{lbtests.df}{Table of the Ljung-Box tests for some lags of the residuals, obtained with the \code{\link[TSA]{LB.test}} function of the TSA package.}
#'    \item{lbtests.plot}{Barplot of the Ljung-Box pvalues of the residuals.}
#'    \item{residRM}{Range-mean plot of the residuals obtained with the \code{\link{rmplot}} function.}
#'    \item{residQQp}{Qqplot of the residuals.}
#'    \item{residDensity}{Plot of the residual density compared with the normal density.}
#'    \item{residMuTest}{T-test for the residual's mean, obtained with \code{\link[stats]{t.test}} function of the stats package.}
#'    \item{residShapiro}{Shapiro-Wilk normality test for the residuals, obtained with \code{\link[stats]{shapiro.test}} function of the stats package.}
#' 
#' @seealso \code{\link[lmtest]{coeftest}}, \code{\link[gridExtra]{grid.arrange}}
#' @examples 
#' #Some examples with AirPassengers
#' model1 <- CArima(as.vector(log(AirPassengers)), order = c(0, 1, 1), 
#'                  seasonal = list(order=c(0,0,1), period = 12))
#' model1 
#' summary(model1, plot=F)
#' 
#' model2 <- CArima(as.vector(log(AirPassengers)), order = c(0, 1, 1), 
#'                  seasonal = list(order=c(0,1,1), period = 12))
#' model2 
#' summary(model2, plot=F) 
#' @import lmtest
#' @import TSA
#' @import gridExtra
#' @import ggplot2
#' @export
CArima <- function(y, order=c(0, 0, 0), seasonal=c(0, 0, 0),plot = T, xreg=NULL,
                   include.constant = TRUE, lambda=model$lambda, biasadj=FALSE,
                   method=c("CSS-ML", "ML", "CSS"), model=NULL, x=y, ...){
  
  series <- deparse(substitute(y))
  digits <- max(3L, getOption("digits") - 3L)
  
  
  if(missing(seasonal)){
    period <- 1
    D <- 0
    method <- "ML"
  }
  else if (!is.list(seasonal)){
    D <- seasonal[2]
    if (frequency(x) > 1){
      period <- frequency(x)
    }
  } 
  else {
    D <- seasonal$order[2]
    period <- seasonal$period
  }
  suppressWarnings(tmp <- forecast::Arima(y = y, order = order, seasonal = seasonal, 
                                            xreg = xreg, include.constant = include.constant,
                                            lambda = lambda, biasadj = biasadj,
                                            method = method, model = model))
  
  tmp$period <- period
  if (D > 0 & !is.na(tmp$coef['drift'])){
    tmp$coef['drift'] <- tmp$coef['drift']*period
    tmp$var.coef[,'drift'] <- tmp$var.coef[,'drift']*period
    tmp$var.coef['drift','drift'] <- tmp$var.coef['drift','drift']*period
  }
  
  if(is.element("drift", names(tmp$coef))){
    names(tmp$coef)[length(names(tmp$coef))] <- "mu"
    colnames(tmp$var.coef)[length(colnames(tmp$var.coef))]<-"mu"
    rownames(tmp$var.coef)[length(rownames(tmp$var.coef))]<-"mu"
  }
  
  if (length(tmp$coef)) {
    #t tests for the coefs. 
    degrees<-length(tmp$residuals)-length(tmp$coef)
    tmp$coef<-coeftest(tmp, df=degrees)
    
    #SBC
    tmp$SBC<-(-2*tmp$loglik + log(length(tmp$residuals))*length(tmp$coef))
    
    #correlations 
    tmp$cor.coef<-cov2cor(tmp$var.coef)
  }
  
  lag.max <- period*round((length(tmp$residuals)/4)/period)
  if(plot){
    par(mfrow = c(2,1))
  }
  #residual PACF and ACF
  tmp$resid.acf <- Acf(tmp$residuals, lag.max = lag.max, main = "Residual ACF", plot = plot)
  tmp$resid.pacf <- Pacf(tmp$residuals, lag.max = lag.max, main = "Residual PACF", plot = plot)
  
  #Ljung-Box tests for the residuals
  lb.df <- data.frame(Lag = c(), ChiSq = c(), DF = c(), pval = c())
  
  for(i in seq(from = 6, to = round(length(tmp$x)/4)+period, by = 6)){
    lbtest <- LB.test(tmp, type = "Ljung-Box", lag=i)
    lb.df <- rbind(lb.df, c(lbtest$lag, lbtest$statistic, lbtest$parameter, lbtest$p.value))
  }
  colnames(lb.df) <- c("Lag", "ChiSq", "DF", "pval")
  tmp$lbtests.df <- lb.df
  
  lb.plot <- ggplot(lb.df, aes(x = factor(Lag), y = pval))+
    geom_bar(stat="identity", color = "black", fill = "lightblue")+
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 0.7)+
    xlab("Lag")+
    ggtitle("Ljung-Box p-values for residuals")
  
  tmp$lbtests.plot <- lb.plot
  
  #residual range-mean plot
  if(period != 1)
    residRM <- rmplot(tmp$residuals, n = period)
  else
    residRM <- rmplot(tmp$residuals)
    
  tmp$residRM <- residRM
  if(plot){
    suppressMessages(grid.arrange(lb.plot, residRM, nrow = 2))
  }
  #residual qqplot
  residQQplot <- ggplot(tmp$residuals, aes(sample = tmp$residuals))+
    stat_qq()+
    stat_qq_line()+
    xlab("Quantiles")+
    ylab("Residuals")+
    ggtitle("Residuals QQ-Plot")
  
  tmp$residQQp <- residQQplot
  
  #residual density plot
  X<-seq(min(tmp$residuals),max(tmp$residuals),length.out=length(tmp$residuals))
  df<-data.frame(x = X, y = dnorm(X, 0, sd(tmp$residuals)))
  
  residDensity <- ggplot(tmp$residuals, aes(x=tmp$residuals))+
    geom_histogram(aes(y = ..density..), color = "black", fill = "lightblue")+
    geom_density(alpha = 0.2, aes(color = "Resid."), size = 0.75)+
    geom_line(data = df, aes(x = x, y = y, color = "Normal"), size = 0.75)+
    xlab("Residual")+
    ggtitle("Residuals Density")+
    labs(color='Density Curve') 
  
  tmp$residDensity <- residDensity
  if(plot){
    suppressMessages(grid.arrange(residQQplot, residDensity, nrow = 2))
  }
  #test residual mean = 0
  tmp$residMuTest <- t.test(tmp$residuals, mu = 0)
  
  #Shapiro Wilks normality test for residuals
  tmp$residShapiro <- shapiro.test(tmp$residuals)
  
  tmp$series <- series
  class(tmp) <- "CArima"
  
  return(tmp)
}

#Modified version of function arima.string from StMoMo package
arima.string <- function (object, padding = FALSE) 
{
  order <- object$arma[c(1, 6, 2, 3, 7, 4, 5)]
  result <- paste("ARIMA(", order[1], ",", order[2], ",", order[3], ")", sep = "")
  
  if (order[7] > 1 & sum(order[4:6]) > 0) 
    result <- paste(result, "(", order[4], ",", order[5], ",", order[6], ")[", order[7], "]", sep = "")
  
  if (is.element("mu", rownames(object$coef)) | is.element("intercept", rownames(object$coef)))
    result <- paste(result, "with constant")
  else result <- paste(result, "                  ")
  
  if (!padding) 
    result <- gsub("[ ]*$", "", result)
  
  return(result)
}

# Modified version of function print.forecast_ARIMA from forecast package
#'@export
print.CArima <- function(x){
  
  digits <- max(3L, getOption("digits") - 3L)
  cat("Series:", x$series, "\n")
  cat(arima.string(x, padding = FALSE), "\n")
  if (!is.null(x$lambda)) {
    cat("Box Cox transformation: lambda=", x$lambda, "\n")
  }
  print(x$coef)
  
  cat("\nsigma^2 = ", format(x$sigma2, digits = digits), sep="")
  if(!is.na(x$loglik))
    cat(":  log likelihood = ", format(round(x$loglik, 2L)), sep = "")
  cat("\n")
  cat("AIC=", format(round(x$aic, 2L)), sep = "")
  cat("   BIC=", format(round(x$bic, 2L)), sep = "")
  cat("   SBC=", format(round(x$SBC, 2L)), "\n", sep = "")
  
  cat("\nCorrelations of the estimated parameters:\n")
  print.default(x$cor.coef)
}

# Modified version of function summary for the CArima function output
#'@export
summary.CArima <- function(x, plot = T){
  
  print(x)
  cat("\n\nLjung-Box tests for residuals:\n")
  print(x$lbtests.df, row.names = F)
  
  if(plot){
    
    par(mfrow = c(2,1))
  }
  lag.max <- x$period*round((length(x$residuals)/4)/x$period)
  acftmp <- Acf(x$residuals, lag.max = lag.max, main = "Residual ACF", plot = plot)
  pacftmp <- Pacf(x$residuals, lag.max = lag.max, main = "Residual PACF", plot = plot)
  
  if(plot){
    suppressMessages(grid.arrange(x$lbtests.plot, x$residRM, nrow = 2))
  }
  
  if(plot){
    suppressMessages(grid.arrange(x$residQQp, x$residDensity, nrow = 2))
  }
  
  cat("\n\nT-test for residuals mean:\n\n")
  cat("H0: mean = 0\n")
  cat("H1: mean != 0\n")
  cat("t = ", x$residMuTest$statistic, ", df = ", x$residMuTest$parameter, ", p-value = ", x$residMuTest$p.value, "\n")
  cat("Confidence interval 95%:\n")
  cat(x$residMuTest$conf.int[1], "  ", x$residMuTest$conf.int[2])
  
  cat("\n\nShapiro-Wilk normality test for residuals:\n\n")
  cat("W = ", x$residShapiro$statistic, ", p-value = ", x$residShapiro$p.value, "\n")
}

