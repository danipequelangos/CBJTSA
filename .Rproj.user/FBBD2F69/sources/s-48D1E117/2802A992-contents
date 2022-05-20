#' @title Forecasting Performance in ARIMA Models
#' 
#' @description Performs more complete analysis of the forecasting performance of an ARIMA model by fitting the model with a train set with \code{h} less points and then use those \code{h} points as a test set.
#' @param object a time series object or a matrix containing the time series values.
#' @param h number of points that will not be used to fit the model, but to use as testing. The points will be taken from the bottom.
#' @param model an object describing a time series model; e.g. of class \code{\link{CArima}}.
#' @param dates a vector or an array with the dates as \code{character} of each observation of the time series. This vector can be generated using the function \code{\link{dateSeq}}.
#' @param plot if \code{TRUE}, it will show all the plots in execution.
#' @return A list containing at least the following components:
#'    \item{forecasts}{The forecasts of the \code{h} points extracted from the original data object, with their confidence intervals as a \code{data.frame}.}
#'    \item{SSEtable}{A \code{data.frame} containing the fitted values, real values, residual values and SSE values of the \code{h} points extracted from the original data object.}
#'    \item{forecPlot}{A plot of the forecasts values and the confidence intervals calculated.}
#'    \item{bandWidth}{The width of the forecasts confidence intervals.}
#'    \item{bandWidthPlot}{A plot of the forecasts confidence bands width.}
#' @import forecast
#' @import ggplot2
#' @seealso \code{\link[ggplot2]{ggplot}}.
#' @examples 
#' #Easy example with AirPassengers
#' dates <- dateSeq(from = time(AirPassengers)[1], length.out = length(AirPassengers),
#'                  by = "month")
#' model <- CArima(as.vector(log(AirPassengers)), order = c(0, 1, 1), 
#'                 seasonal = list(order=c(0,1,1), period = 12))
#' model
#' 
#' predCapModel <- predCap(as.vector(log(AirPassengers)), h = 12, model, dates, 
#'                         plot = F) 
#' predCapModel
#' @export
predCap <- function(object, h = 10, model, dates, plot = T){
  
  result <- list()
  
  newObject <- object[1:(length(object)-h)]
  order <- model$arma[c(1, 6, 2, 3, 7, 4, 5)]
  if (is.element("mu", rownames(model$coef)))
    include.constant <- TRUE
  else
    include.constant <- FALSE
  
  model <- suppressWarnings(forecast::Arima(newObject, order = c(order[1], order[2], order[3]), 
                                            seasonal = list(order = c(order[4], order[5], order[6]), period = order[7]),
                                            include.constant = include.constant))
  
  tmp <- suppressWarnings(forecast::forecast(newObject, model = model, h = h))
  
  predictions <- as.data.frame(tmp)
  rownames(predictions) <- dates[(length(dates)-h+1):length(dates)]
  result$forecasts <- predictions
  
  #Fitted
  SSEtable <- cbind(predictions$`Point Forecast`)
  #Real
  SSEtable <- cbind(SSEtable, object[(length(object)-h+1):length(object)])
  #Residuals
  SSEtable <- cbind(SSEtable, abs(SSEtable[,1] - SSEtable[,2]))
  SSE <- c(SSEtable[1,3]**2)
  for(i in 2:h){
    SSE <- append(SSE, ((SSEtable[i,3]**2) + SSE[i-1]))
  }
  SSEtable <- cbind(SSEtable, SSE)
  colnames(SSEtable) <- c("Fit", "Real", "Res", "SSE")
  
  result$SSEtable <- as.data.frame(SSEtable, row.names = 
                                     dates[(length(dates)-h+1):length(dates)])
  
  #Forecast plot
  if(h > 24){
    n.unlabeled <- 1
  }else{
    n.unlabeled <- 0
  }
  
  predDates <- dates[(length(dates)-h+1):length(dates)]
  datesOrd <- ordered(predDates, levels = predDates)
  forecPlot <- ggplot(predictions, aes(x=factor(datesOrd), y=`Point Forecast`, group = 1)) +
    geom_line()+
    theme(axis.text.x=element_text(angle=60,hjust=1))+
    scale_x_discrete(breaks = levels(datesOrd)[c(T, rep(F, n.unlabeled))])+
    ggtitle("Forecasts plot")+
    xlab("dates")+
    ylab("forecasts")+
    geom_ribbon(aes(ymin=`Lo 95`, ymax = `Hi 95`), fill = "dodgerblue", alpha = 0.6)+
    geom_ribbon(aes(ymin=`Lo 80`, ymax = `Hi 80`), fill = "dodgerblue4", alpha = 0.5)
  
  if(plot){
    print(forecPlot)
  }
  result$forecPlot <- forecPlot
  
  bandWidth <- predictions$`Hi 95`-predictions$`Lo 95`
  result$bandWidth <- bandWidth
  
  bandWidthPlot <- ggplot(predictions, aes(x=factor(datesOrd), y = bandWidth, group=1))+
    geom_line()+
    geom_point()+
    theme(axis.text.x=element_text(angle=60,hjust=1))+
    scale_x_discrete(breaks = levels(datesOrd)[c(T, rep(F, n.unlabeled))])+
    ggtitle("Forecasts band width plot")+
    ylab("band width")+
    xlab("dates")
  
  result$bandWidthPlot <- bandWidthPlot
  
  return(result)
}

#' @title Compare ARIMA Models
#' 
#' @description Compare up to three ARIMA models of the class \code{\link{CArima}}, with theirs metrics (Sigma square, AIC and SBC) and their forecasting performance.
#' @param object a time series object or a matrix containing the time series values.
#' @param dates a vector or an array with the dates as \code{character} of each observation of the time series. This vector can be generated using the function \code{\link{dateSeq}}.
#' @param m1 an object describing a time series model; e.g. of class \code{\link{CArima}}.
#' @param m2 an object describing a time series model; e.g. of class \code{\link{CArima}}.
#' @param m3 an object describing a time series model; e.g. of class \code{\link{CArima}}.
#' @param h number of points that will not be used to fit the model, but to use as testing. The points will be taken from the bottom.
#' @param plot if \code{TRUE}, it will show all the plots in execution.
#' @return An object of class \code{compareARIMA} containing:
#'    \item{metrics}{Some metrics of the models, such as AIC, Sigma square, SBC, SSE of prediction and two Ljung-Box test values.}
#'    \item{SSEtable}{A merge of the \code{SSEtable} of m1 and m2 obtained with \code{\link{predCap}} function.}
#'    \item{forecasts}{A merge of the \code{forecasts} of m1 and m2 obtained with \code{\link{predCap}} function.}
#'    \item{forecPlot1}{A plot of the forecasts values and the confidence intervals calculated of m1 obtained with \code{\link{predCap}} function.}
#'    \item{forecPlot2}{A plot of the forecasts values and the confidence intervals calculated of m2 obtained with \code{\link{predCap}} function.}
#'    \item{forecPlot3}{A plot of the forecasts values and the confidence intervals calculated of m3 obtained with \code{\link{predCap}} function.}
#'    \item{bandWidthComp}{A data frame with the forecasts band widths of each model in a column.}
#'    \item{bandWidthPlot}{A plot of the forecasts band widths of each model.}
#' @import gridExtra
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link{predCap}}.
#' @examples 
#' #Some examples with AirPassengers
#' dates <- dateSeq(from = time(AirPassengers)[1], length.out = length(AirPassengers),
#'                  by = "month")
#' model1 <- CArima(as.vector(log(AirPassengers)), order = c(0, 1, 1), 
#'                  seasonal = list(order=c(0,0,1), period = 12))
#' model1
#' model2 <- CArima(as.vector(log(AirPassengers)), order = c(0, 1, 1), 
#'                  seasonal = list(order=c(0,1,1), period = 12))
#' model2 
#' model3 <- CArima(as.vector(log(AirPassengers)), order = c(1, 1, 0), 
#'                  seasonal = list(order=c(1, 1, 0), period = 12))
#' model3
#' 
#' comp <- compareARIMA(as.vector(log(AirPassengers)), dates, model1, model2, 
#'                      model3, h = 12, plot = F)
#' comp
#' @export
compareARIMA <- function(object, dates, m1, m2, m3, h = 10, plot = T){
  
  result <- list()
  metrics <- data.frame()
  
  predCap1 <- predCap(object, h, m1, dates, plot = F)
  result$forecPlot1 <- predCap1$forecPlot+
    ggtitle("Forecasts plot for m1")
  
  predCap2 <- predCap(object, h, m2, dates, plot = F)
  result$forecPlot2 <- predCap2$forecPlot+
    ggtitle("Forecasts plot for m2")
  
  bandWidthComp <- cbind(predCap1$bandWidth, predCap2$bandWidth)
  
  predDates <- dates[(length(dates)-h+1):length(dates)]
  datesOrd <- ordered(predDates, levels = predDates)
  
  if(h > 24){
    n.unlabeled <- 1
  }else{
    n.unlabeled <- 0
  }
  
  if(!missing(m3)){
    predCap3 <- predCap(object, h, m3, dates, plot = F)
    result$forecPlot3 <- predCap3$forecPlot+
      ggtitle("Forecasts plot for m3")
    
    bandWidthComp <- cbind(bandWidthComp, predCap3$bandWidth)
    colnames(bandWidthComp) <- c("Band Width m1", "Band Width m2", "Band Width m3")
    rownames(bandWidthComp) <- predDates
    result$bandWidthComp <- bandWidthComp
    
    bandWidthPlot <- ggplot(predCap1$forecasts, aes(x=factor(datesOrd), y = predCap1$bandWidth, group=1, color = "m1"))+
      geom_line()+
      geom_point()+
      geom_line(aes(x=factor(datesOrd), y=predCap2$bandWidth, color="m2"))+
      geom_point(aes(x=factor(datesOrd), y=predCap2$bandWidth, color="m2"))+
      geom_line(aes(x=factor(datesOrd), y=predCap3$bandWidth, color="m3"))+
      geom_point(aes(x=factor(datesOrd), y=predCap3$bandWidth, color="m3"))+
      theme(axis.text.x=element_text(angle=60,hjust=1))+
      scale_x_discrete(breaks = levels(datesOrd)[c(T, rep(F, n.unlabeled))])+
      ggtitle("Forecasts band width plot")+
      ylab("band width")+
      xlab("dates")
    result$bandWidthPlot <- bandWidthPlot
  }
  else{
    colnames(bandWidthComp) <- c("Band Width m1", "Band Width m2")
    rownames(bandWidthComp) <- predDates
    result$bandWidthComp <- bandWidthComp
    
    bandWidthPlot <- ggplot(predCap1$forecasts, aes(x=factor(datesOrd), y = predCap1$bandWidth, group=1, color = "m1"))+
      geom_line()+
      geom_point()+
      geom_line(aes(x=factor(datesOrd), y=predCap2$bandWidth, color="m2"))+
      geom_point(aes(x=factor(datesOrd), y=predCap2$bandWidth, color="m2"))+
      theme(axis.text.x=element_text(angle=60,hjust=1))+
      scale_x_discrete(breaks = levels(datesOrd)[c(T, rep(F, n.unlabeled))])+
      ggtitle("Forecasts band width plot")+
      ylab("band width")+
      xlab("dates")
    result$bandWidthPlot <- bandWidthPlot
  }
  
  if(plot){
    if(!missing(m3)){
      suppressMessages(grid.arrange(result$forecPlot1, result$forecPlot2, 
                                    result$forecPlot3, nrow = 3))
    }
    else{
      suppressMessages(grid.arrange(result$forecPlot1, result$forecPlot2, nrow = 2))
    }
    print(result$bandWidthPlot)
  }
  
  metrics <- cbind(c(m1$sigma2, m1$aic, m1$SBC, 
                     predCap1$SSEtable$SSE[NROW(predCap1$SSEtable$SSE)], 
                     m1$lbtests.df$pval[1], m1$lbtests.df$pval[length(m1$lbtests.df$pval)-1]))
  
  metrics <- cbind(metrics, c(m2$sigma2, m2$aic, m2$SBC, 
                     predCap2$SSEtable$SSE[NROW(predCap2$SSEtable$SSE)],
                     m2$lbtests.df$pval[1], m2$lbtests.df$pval[length(m1$lbtests.df$pval)-1]))
  
  if(!missing(m3)){
    metrics <- cbind(metrics, c(m3$sigma2, m3$aic, m3$SBC, 
                                predCap3$SSEtable$SSE[NROW(predCap3$SSEtable$SSE)],
                                m3$lbtests.df$pval[1], m3$lbtests.df$pval[length(m1$lbtests.df$pval)-1]))
  }
  
  rownames(metrics) <- c("Sigma2", "AIC", "SBC", "SSEpred", 
                         paste("Ljung-Box Lag",m1$lbtests.df$Lag[1]), 
                         paste("Ljung-Box Lag",m1$lbtests.df$Lag[length(m1$lbtests.df$DF)-1]))
  if(!missing(m3))
    colnames(metrics) <- c("m1", "m2", "m3")
  else
    colnames(metrics) <- c("m1", "m2")
  
  result$metrics <- metrics
  
  if(!missing(m3)){
    SSEtable <- cbind(predCap1$SSEtable$Fit, predCap2$SSEtable$Fit, predCap3$SSEtable$Fit, 
                      predCap1$SSEtable$Real, predCap1$SSEtable$Res, predCap2$SSEtable$Res, 
                      predCap3$SSEtable$Res, predCap1$SSEtable$SSE, predCap2$SSEtable$SSE,
                      predCap3$SSEtable$SSE)
    colnames(SSEtable) <- c("Fit_m1", "Fit_m2", "Fit_m3", "Real", "Res_m1", "Res_m2",
                            "Res_m3", "SSE_m1", "SSE_m2", "SSE_m3")
  }
  else{
    SSEtable <- cbind(predCap1$SSEtable$Fit, predCap2$SSEtable$Fit, predCap1$SSEtable$Real, 
                      predCap1$SSEtable$Res, predCap2$SSEtable$Res, predCap1$SSEtable$SSE, 
                      predCap2$SSEtable$SSE)
    colnames(SSEtable) <- c("Fit_m1", "Fit_m2", "Real", "Res_m1", "Res_m2", "SSE_m1", 
                            "SSE_m2")
  }
  
  rownames(SSEtable) <- predDates
  
  result$SSEtable <- SSEtable
  
  colnames(predCap1$forecasts)[1]<-"Forecast"
  colnames(predCap2$forecasts)[1]<-"Forecast"
  if(!missing(m3)){
    colnames(predCap3$forecasts)[1]<-"Forecast"
    forecasts <- data.frame(m1 = predCap1$forecasts, space = rep("", NROW(predCap1$forecasts)),
                            m2 = predCap2$forecasts, space = rep("", NROW(predCap1$forecasts)),
                            m3 = predCap3$forecasts)
    colnames(forecasts)[12] <- ""
  }
  else{
    forecasts <- data.frame(m1 = predCap1$forecasts, space = rep("", NROW(predCap1$forecasts)),
                            m2 = predCap2$forecasts)
  }
  colnames(forecasts)[6] <- ""
  rownames(forecasts) <- predDates
  
  result$forecasts <- forecasts
  
  class(result) <- "compareARIMA"
  return(result)
}

# Modified version of function print for the class compareARIMA
#' @export
print.compareARIMA <- function(x){
  cat("\n\nComparison of metrics:\n\n")
  print(x$metrics)
  
  cat("\n\nComparison of SSE tables:\n\n")
  print(x$SSEtable)
  
  cat("\n\nComparison of forecasts:\n\n")
  print(x$forecasts)
  
  cat("\n\nComparison of forecasts band widths:\n\n")
  print(x$bandWidthComp)
  
  if(NCOL(x$metrics) == 3){
    suppressMessages(grid.arrange(x$forecPlot1, x$forecPlot2, x$forecPlot3, nrow = 3))
    print(x$bandWidthPlot)
  }
  else{
    suppressMessages(grid.arrange(x$forecPlot1, x$forecPlot2, nrow = 2))
    print(x$bandWidthPlot)
  }
}

#' @title Forecasting time series
#' 
#' @description Largely a wrapper for the \code{\link[forecast]{forecast}} function in the forecast package. It adds some usefull datas and graphs to validate and compare models.
#' @param object a time series object or a matrix containing the time series values.
#' @param h number of periods for forecasting.
#' @param model an object describing a time series model; e.g. of class \code{\link{CArima}}.
#' @param plot if \code{TRUE}, it will show all the plots in execution.
#' @param dates a vector or an array with the dates as \code{character} of each observation of the time series. This vector can be generated using the function \code{\link{dateSeq}}.
#' @param ... additional arguments to be passed to \code{\link[forecast]{forecast}}.
#' @details See the \code{\link[forecast]{forecast}} function in the forecast package.
#' @return See the \code{\link[forecast]{forecast}} function in the forecast package. The additional objects returned are: 
#'    \item{forecasts}{The forecasts with their confidence intervals as a \code{data.frame}.}
#'    \item{realForecPlot}{A real + forecast plot.}
#'    \item{realFitPlot}{A real vs fitted values plot.}
#' @seealso \code{\link[ggplot2]{ggplot}}
#' @examples 
#' #An easy example with AirPassengers
#' dates <- dateSeq(from = time(AirPassengers)[1], length.out = length(AirPassengers),
#'                  by = "month")
#' model <- CArima(as.vector(log(AirPassengers)), order = c(0, 1, 1), 
#'                  seasonal = list(order=c(0,1,1), period = 12))
#' model   
#'                 
#' predModel <- Cforecast(as.vector(log(AirPassengers)), model = model, h = 24, 
#'                        dates = dates, plot = F) 
#' predModel
#' summary(predModel)
#' @export
Cforecast <- function(object, h=ifelse(frequency(object) > 1, 2 * frequency(object), 10),
                      level=c(80, 95), fan=FALSE, robust=FALSE, lambda = NULL, 
                      biasadj = FALSE, find.frequency = FALSE, allow.multiplicative.trend=FALSE, 
                      model=NULL, plot = TRUE, dates, ...){
  
  result <- model
  order <- model$arma[c(1, 6, 2, 3, 7, 4, 5)]
  if (is.element("mu", rownames(model$coef)))
    include.constant <- TRUE
  else
    include.constant <- FALSE
  model <- suppressWarnings(forecast::Arima(object, order = c(order[1], order[2], order[3]), 
                                            seasonal = list(order = c(order[4], order[5], order[6]), period = order[7]),
                                            include.constant = include.constant))
  tmp <- suppressWarnings(forecast::forecast(object, model = model, h = h, level = level,
                                             fan = fan, robust = robust, lambda = lambda,
                                             biasadj = biasadj, find.frequency = find.frequency,
                                             allow.multiplicative.trend = allow.multiplicative.trend))
  
  result$forecasts <- as.data.frame(tmp)
  
  datesExtra <- dateSeq(from = dates[length(dates)], length.out = h+1, by = "month")
  datesExtra <- datesExtra[2:length(datesExtra)]
  rownames(result$forecasts) <- datesExtra
  
  datesOrd <- ordered(append(dates,datesExtra), levels = append(dates,datesExtra))
  xTot <- append(result$x, result$forecasts$`Point Forecast`)
  realForecPlot <- ggplot(as.data.frame(xTot), aes(x=factor(datesOrd), y=xTot, group = 1)) +
    geom_line()+
    theme(axis.text.x=element_text(angle=60,hjust=1))+
    scale_x_discrete(breaks = levels(datesOrd)[c(T, rep(F, 9))])+
    ggtitle("Real + Forecasts plot")+
    xlab("dates")+
    ylab("real and forecasts")+
    geom_ribbon(aes(ymin=append(result$x, result$forecasts$`Lo 95`), ymax = append(result$x, result$forecasts$`Hi 95`)), fill = "dodgerblue", alpha = 0.6)+
    geom_ribbon(aes(ymin=append(result$x, result$forecasts$`Lo 80`), ymax = append(result$x, result$forecasts$`Hi 80`)), fill = "dodgerblue4", alpha = 0.5)
  
  if(plot){
    print(realForecPlot)
  }
  result$realForecPlot <- realForecPlot

  datesOrdReal <- ordered(dates, levels = dates)
  xPred <- append(result$fitted, result$forecasts$`Point Forecast`)
  realFitPlot <- ggplot(as.data.frame(xPred), aes(x=factor(datesOrd), y=xPred, group = 1)) +
    geom_line(aes(color = "Fitted"), size= 0.6)+
    geom_line(data = as.data.frame(result$x), aes(x=factor(datesOrdReal), y=result$x, color="Real"), size= 0.6)+
    theme(axis.text.x=element_text(angle=60,hjust=1))+
    scale_x_discrete(breaks = levels(datesOrd)[c(T, rep(F, 9))])+
    xlab("dates")+
    ylab("values")+
    ggtitle("Real vs Fitted plot")+
    labs(color='Series')   
  
  if(plot){
    print(realFitPlot)
  }
  result$realFitPlot <- realFitPlot
  
  class(result) <- "Cforecast"
  return(result)
  
}

# Modified version of function print.forecast_ARIMA from forecast package
#'@export
print.Cforecast <- function(x){
  
  print.CArima(x)
}

# Modified version of function summary for the Cforecast function output
#'@export
summary.Cforecast <- function(x, plot = T){
  
  summary.CArima(x)
  cat("\n\nForecasts:\n\n")
  print(x$forecasts)
  
  if(plot){
    print(x$realForecPlot)
    print(x$realFitPlot)
  }
}
