#' @title Time Series Plot
#' 
#' @description \code{tsplot} creates a \code{\link[ggplot2]{ggplot}} object showing the main plot of a time series.
#' @param y a vector or an array with the values of the time series.
#' @param x a vector or an array with the dates as \code{character} of each observation of the time series. This vector can be generated using the function \code{\link{dateSeq}}.
#' @param n.unlabel number of unlabeled dates between labeled dates in the plot.
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @examples 
#' # Easy example with AirPassengers time series
#' dates <- dateSeq(from = time(AirPassengers)[1], length.out = length(AirPassengers), 
#'                  by = "month")
#' tsplot(dates, as.vector(AirPassengers))
#' @import ggplot2
#' @export
tsplot <- function(x, y, n.unlabel = 9){
  datesOrd <- ordered(x, levels=x)
  df <- data.frame(x = factor(datesOrd), y = y)
  p <- ggplot(df, aes(x = x, y = y, group = 1)) +
    geom_line()+
    geom_point(color='red') +
    theme(axis.text.x=element_text(angle=60,hjust=1))+
    scale_x_discrete(breaks = levels(datesOrd)[c(T, rep(F, n.unlabel))])+
    xlab("Dates")+
    ylab("Values")+
    ggtitle("Time Series Plot")
  return(p)
}

#' @title Time Series Frequency Spectrum
#' 
#' @description \code{tsspec} generates a ggplot showing the frequency spectrum of a time series, and a table with the value of each frequency.
#' @param y a vector or an array with the values of the time series.
#' @return A list containing the following components: 
#'    \item{table}{A data frame with the values of the spectrum for each frequency represented.}
#'    \item{plot}{A \code{\link[ggplot2]{ggplot}} object of the frequency spectrum.}
#' @seealso \code{\link[TSA]{periodogram}}
#' @examples 
#' # Easy example with AirPassengers time series
#' tsspec(as.vector(AirPassengers))
#' @import TSA
#' @export
tsspec <- function(y){
  spec <- periodogram(y, plot = F)
  table <- data.frame(freq = spec$freq, value = spec$spec)
  
  plot <- ggplot(table, aes(x = freq, y = value))+
    geom_line()+
    geom_point(color = 'red')+
    ggtitle("Time Series Frequency Spectrum")
  return(list(table = table, plot = plot))
}

#' @title Time Series Range-Mean Plot
#' 
#' @description \code{rmplot} generates a ggplot object showing a range-mean plot of a time series.
#' @param y a vector or an array with the values of the time series.
#' @param n size of the subsamples to calculate the mean and range, this parameter should be equals to the period of the time series.
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @examples 
#' # Easy example with AirPassengers time series
#' rmplot(as.vector(AirPassengers), n = 12)
#' rmplot(log(as.vector(AirPassengers)), n = 12)
#' @export
rmplot <- function(y, n=5){
  groups <- split(y, rep(1:round(length(y)/n), each=n))
  means <- lapply(groups,mean)
  ranges <- lapply(groups,var)
  mrdata <- data.frame(mean=unlist(means), range=unlist(ranges))
  
  p <- ggplot(mrdata, aes(x=mean,y=range))+
    geom_point()+
    ggtitle("Range-Mean Plot")
  return(p)
}

#' @title Simulate ARIMA Model and Compare ACF and PACF
#' 
#' @description Simulate an ARIMA model of a specified size, showing its ACF, PACF and time series plot, and compare it with a time series sample if wanted.
#' @param dates optional: a vector or an array with the dates as \code{character} of each observation of the time series. This vector can be generated using the function \code{\link{dateSeq}}.
#' @param sample optional: a vector or an array with the values of the time series without differencing.
#' @param order a specification of the non-seasonal part of the ARIMA model: the three components (p, d, q) are the AR order (p and q maximum of 3), the order of differencing, and the MA order. See 'Details'.
#' @param seasonal a specification of the seasonal part of the ARIMA model: the three components (P, D, Q) are the AR order (P and Q maximum of 3), the order of differencing, and the MA order.
#' @param n the size of the simulation, 10000 by default since that size makes the simulation very similar to the theorical ARIMA model, which can be seen in the ACF and PACF representations.
#' @param period the period of the seasonal part of the ARIMA model.
#' @param ar the AR coefficients of the non-seasonal part of the ARIMA model in order. For example ar = c(0.5, 0.7), 0.5 coefficient of the first lag, 0.7 coefficient of the second lag. See 'Details'. 
#' @param ma the MA coefficients of the non-seasonal part of the ARIMA model in order. Works the same way as \code{ar} parameter.
#' @param sar the AR coefficients of the seasonal part of the ARIMA model in order. For example sar = c(-0.3, 0.6), -0.3 coefficient of the \code{period}-th lag, 0.6 coefficient of the \code{period*2}-th lag. See 'Details'.
#' @param sma the MA coefficients of the seasonal part of the ARIMA model in order. Works the same way as \code{sar} parameter.
#' @param lag.max maximum number of lags at which to calculate the ACF. Default is \code{period*round((length(sample)/4)/period)}.
#' @param plot if \code{TRUE}, it will show all the plots in execution.
#' @details The function uses \code{\link[stats]{arima.sim}} to simulate the ARIMA model. The order of differencing specified in \code{order} and \code{seasonal} parameters will affect only the sample, which will be differenced if those orders are greater than 0.
#' 
#' If a multiplicative ARIMA model is chosen the coefficients product of the non-seasonal and seasonal part will be calculated automatically with the coefficients of the main lags.
#' @return A list containing the following components: 
#'    \item{simplot}{A plot of the simulated model as specified.}
#'    \item{sampplot}{A plot of the original sample time series if given.}
#'    \item{diffsampplot}{A plot of the differenced sample time series if specified so.}
#'    \item{acfPlot}{An ACF plot of the simulated model and (if specified) the time series given differenced as specified.}
#'    \item{pacfPlot}{A PACF plot of the simulated model and (if specified) the time series given differenced as specified.}
#' @seealso \code{\link[gridExtra]{grid.arrange}}, \code{\link[ggplot2]{ggplot}}.
#' @examples 
#' #Easy example with AirPassengers 
#' dates <- dateSeq(from = time(AirPassengers)[1], to = time(AirPassengers)[length(AirPassengers)], 
#'                  by="month")
#'
#' arimaSim <- arimaSimComp(dates, sample = log(AirPassengers), order = c(0, 1, 1), seasonal = c(0, 1, 1),
#'              ma = c(-0.4), period = 12, sma = c(-0.6), plot = F)
#' arimaSim 
#' 
#' #Example without a sample
#' arimaSim2 <- arimaSimComp(order = c(0, 0, 2), seasonal = c(0, 0, 1), period = 12, ma = c(-0.5, 0.4),
#'              sma = c(-0.6), plot = F)
#' arimaSim2
#' @import forecast
#' @import gridExtra
#' @export
arimaSimComp <- function(dates, sample, order = c(0, 0, 0), seasonal = c(0, 0, 0), 
                         n = 10000, period = 1, ar, ma, sar, sma, lag.max, plot = T){
  result <- list()
  if(!missing(sample)){
    ogSample <- sample
    ogDates <- dates
  }
  dif <- F
  
  if(!missing(sample)){
    #Diff of the sample
    if(order[2] != 0){
      dif <- T
      sample <- diff(sample, differences = order[2])
      dates <- dates[(order[2]+1):length(dates)]
    }
    if(seasonal[2] != 0){
      dif <- T
      sample <- diff(sample, differences = seasonal[2], lag = period)
      dates <- dates[((seasonal[2]*period)+1):length(dates)]
    }
    if(missing(lag.max))
      lag.max <- period*round((length(sample)/4)/period)
    
    acfSamp<-Acf(sample, lag.max=lag.max, plot=F)
    pacfSamp<-Pacf(sample, lag.max=lag.max, plot=F)
  }
  else{
    if(missing(lag.max))
      lag.max <- 40
  }
  
  if(!missing(ar) & !missing(sar)){
    sar <- append(sar, rep(0, 3-length(sar)))
    arCoef <- c(ar, rep(0, period-order[1]-1), sar[1], ar*sar[1],
                rep(0, (2*period)-order[1]-2), sar[2], ar*sar[2],
                rep(0, (3*period)-order[1]-2), sar[3], ar*sar[3])
  }
  else if(missing(sar) & !missing(ar))
    arCoef <- ar
  else if(missing(ar) & !missing(sar)){
    sar <- append(sar, rep(0, 3-length(sar)))
    arCoef <- c(rep(0, period-1), sar[1], rep(0, (2*period)-1), sar[2], 
                rep(0, (3*period)-1), sar[3])
  }
  else
    arCoef <- c(0)
  
  if(!missing(ma) & !missing(sma)){
    sma <- append(sma, rep(0, 3-length(sma)))
    maCoef <- c(ma, rep(0, period-order[3]-1), sma[1], ma*sma[1],
                rep(0, (2*period)-order[3]-2), sma[2], ma*sma[2],
                rep(0, (3*period)-order[3]-2), sma[3], ma*sma[3])
  }
  else if(missing(sma) & !missing(ma))
    maCoef <- ma
  else if(missing(ma) & !missing(sma)){
    sma <- append(sma, rep(0, 3-length(sma)))
    maCoef <- c(rep(0, period-1), sma[1], rep(0, (2*period)-1), sma[2], 
                rep(0, (3*period)-1), sma[3])
  }
  else
    maCoef <- c(0)
  
  
  suppressWarnings(theorical <- arima.sim(n = n, list(ma = maCoef, ar = arCoef)))
  acftheor<-Acf(theorical, lag.max=lag.max, plot=F)
  pacftheor<-Pacf(theorical, lag.max =lag.max, plot=F)
  
  #Time series plot
  if(!missing(sample)){
    datesOrd <- ordered(dates, levels=dates)
    df <- data.frame(x = factor(datesOrd), y = sample)
    df2 <- data.frame(x = factor(datesOrd), y = as.vector(theorical[1:length(dates)]))
    sampplot <- ggplot() +
      geom_line(df, mapping = aes(x = x, y = y, group = 1))+
      theme(axis.text.x=element_text(angle=60,hjust=1))+
      scale_x_discrete(breaks = levels(datesOrd)[c(T, rep(F, 9))])+
      xlab("Dates")+
      ylab("Values")
    
    simplot <- ggplot() +
      geom_line(df2, mapping = aes(x = x, y = y, group = 1))+
      theme(axis.text.x=element_text(angle=60,hjust=1))+
      scale_x_discrete(breaks = levels(datesOrd)[c(T, rep(F, 9))])+
      xlab("Dates")+
      ylab("Values")+
      ggtitle(paste("Simulation Plot of first", length(dates), "obs"))
    
    if(dif){
      diffsampplot <- sampplot + ggtitle("Differenced Sample Plot")
      
      ogDatesOrd <- ordered(ogDates, levels = ogDates)
      df3 <- data.frame(x = factor(ogDatesOrd), y = ogSample)
      sampplot <- ggplot() +
        geom_line(df3, mapping = aes(x = x, y = y, group = 1))+
        theme(axis.text.x=element_text(angle=60,hjust=1))+
        scale_x_discrete(breaks = levels(ogDatesOrd)[c(T, rep(F, 9))])+
        xlab("Dates")+
        ylab("Values")+
        ggtitle("Original Sample Plot")
      
      result$diffsampplot <- diffsampplot
    }
    else
      sampplot <- sampplot + ggtitle("Original Sample Plot")
    
    result$sampplot <- sampplot
    result$simplot <- simplot
  }
  else{
    simplot <- autoplot(theorical, color="red", main = paste("Simulation Plot of",n,"obs"))
    result$simplot <- simplot
  }
  
  if(plot){
    if(dif)
      suppressMessages(grid.arrange(simplot, diffsampplot, sampplot, nrow = 3))
    else if(!missing(sample))
      suppressMessages(grid.arrange(simplot, sampplot, nrow = 2))
    else
      suppressMessages(print(simplot))
  }
  
  if(!missing(sample)){
    clim0 <- qnorm((1 + 0.95)/2)/sqrt(acfSamp$n.used)
    acf.df <- data.frame(lag=c(acftheor$lag[2:length(acftheor$lag)], acftheor$lag[2:length(acftheor$lag)]), 
                       acf=c(acfSamp$acf[2:length(acfSamp$lag)], acftheor$acf[2:length(acftheor$lag)]), 
                       type=c(rep("Sample",length(acfSamp$lag)-1), rep("Simulation", length(acfSamp$lag)-1)))
    
    acfPlot <- ggplot(acf.df, aes(x = lag, y = acf, fill = type))+
      geom_bar(stat="identity", position = 'dodge', alpha = 0.7)+
      ylim(-1,1)+
      labs(x = "Lags", y = "ACF", title = "ACF Simulation and Sample")+
      geom_hline(yintercept=0, linetype="solid", color="black")+
      geom_hline(yintercept=clim0, linetype="dashed", color="red")+
      geom_hline(yintercept=-clim0, linetype="dashed", color="red")
    
    #PACF
    pacf.df <- data.frame(lag=c(pacftheor$lag[1:length(pacftheor$lag)], pacftheor$lag[1:length(pacftheor$lag)]), 
                        acf=c(pacfSamp$acf[1:length(pacfSamp$lag)],pacftheor$acf[1:length(pacftheor$lag)]), 
                        type=c(rep("Sample",length(pacfSamp$lag)), rep("Simulation", length(pacftheor$lag))))
    
    pacfPlot <- ggplot(pacf.df, aes(x = lag, y = acf, fill = type))+
      geom_bar(stat="identity", position = 'dodge', alpha = 0.7)+
      ylim(-1,1)+
      labs(x = "Lags", y = "PACF", title = "PACF Simulation and Sample")+
      geom_hline(yintercept=0, linetype="solid", color="black")+
      geom_hline(yintercept=clim0, linetype="dashed", color="red")+
      geom_hline(yintercept=-clim0, linetype="dashed", color="red")
  }
  else{
    suppressWarnings(acfPlot <- autoplot(acftheor, ylim = c(-1,1), main="ACF Simulation"))
    suppressWarnings(pacfPlot <- autoplot(pacftheor, ylim = c(-1,1), main="PACF Simulation"))
  }
  
  result$acfPlot <- acfPlot
  result$pacfPlot <- pacfPlot
  
  if(plot)
    suppressMessages(grid.arrange(acfPlot, pacfPlot, nrow = 2))
  
  class(result) <- "arimaSimComp"
  return(result)
}

#Modified version of print function for the arimaSimComp function output
#'@export
print.arimaSimComp<- function(x){
  if(!is.null(x$diffsampplot))
    suppressMessages(grid.arrange(x$simplot, x$diffsampplot, x$sampplot, nrow = 3))
  else if(!is.null(x$sampplot))
    suppressMessages(grid.arrange(x$simplot, x$sampplot, nrow = 2))
  else
    suppressMessages(print(x$simplot))
  
  suppressMessages(grid.arrange(x$acfPlot, x$pacfPlot, nrow = 2))
}

#' @importFrom forecast Acf
#' @export
forecast::Acf

#' @importFrom forecast Pacf
#' @export
forecast::Pacf

