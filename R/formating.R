#' @title Generate Sequences of Dates for Time Series
#' 
#' @description Creates a sequence of dates or times.
#' @param from starting date. 
#' @param to end date. Can be given instead of \code{length.out} parameter
#' @param lenght.out desired length of the sequence. A non-negative integer. Can be given instead of \code{to} parameter.
#' @param by character string, increment of the sequence. One of "\code{min}", "\code{hour}", "\code{day}", "\code{week}", "\code{month}", "\code{quarter}", "\code{year}". See 'Details'
#' @details Depending on the value of \code{by}, there are some requirements for the other arguments.
#' \itemize{
#'    \item For "\code{min}" and "\code{hour}" the format of \code{from} and \code{to} should be a correct format for the function \code{\link[base]{as.POSIXct}}.
#'    \item For "\code{day}" and "\code{week}" the format of \code{from} and \code{to} should be a correct format for the function \code{\link[base]{as.Date}}.
#'    \item For "\code{month}" the format of \code{from} and \code{to} should be a correct format for the function \code{\link[zoo]{yearmon}}.
#'    \item For "\code{quarter}" the format of \code{from} and \code{to} should be a correct format for the function \code{\link[zoo]{yearqtr}}.
#' }
#' @return A vector of dates.
#' @seealso \code{\link[base]{as.POSIXct}}, \code{\link[base]{as.Date}}, \code{\link[zoo]{yearmon}}, \code{\link[zoo]{yearqtr}}.
#' @examples 
#' dateSeq(from = "2022-01", length.out = 24, by = "month")
#' dateSeq(from = "2022-01-01", to = "2023-12", by = "month")
#' dateSeq(from = "2022-01-03", length.out = 12, by = "week")
#' dateSeq(from = "2022-01-03", to = "2022-03-21", by = "week")
#' 
#' #An example extracting the times of AirPassengers
#' dateSeq(from = time(AirPassengers)[1], length.out = length(AirPassengers), 
#'         by = "month")
#' @import zoo
#' @export
dateSeq <- function(from, to, length.out, by){
  if(by == "min"){
    if(missing(length.out)){
      dates <- as.character(strftime(seq(from=as.POSIXct(from),to=as.POSIXct(to),by=by)))
    }
    else{
      dates <- as.character(strftime(seq(from=as.POSIXct(from),length.out=length.out,by=by)))
    }
  } 
  else if(by == "hour"){
    if(missing(length.out)){
      dates <- as.character(strftime(seq(from=as.POSIXct(from),to=as.POSIXct(to),by=by)))
    }
    else{
      dates <- as.character(strftime(seq(from=as.POSIXct(from),length.out=length.out,by=by)))
    }
  }
  else if(by == "day"){
    if(missing(length.out)){
      dates <- as.character(seq(from=as.Date(from), to=as.Date(to), by=by))
    }
    else{
      dates <- as.character(seq(from=as.Date(from), length.out=length.out, by=by))
    }
  }
  else if(by == "week"){
    if(missing(length.out)){
      dates <- as.character(seq(from=as.Date(from), to=as.Date(to), by=by))
    }
    else{
      dates <- as.character(seq(from=as.Date(from), length.out=length.out, by=by))
    }
  }
  else if(by == "month"){
    if(missing(length.out)){
      dates <- as.character(as.yearmon(seq(from=as.numeric(as.yearmon(from)),to=as.numeric(as.yearmon(to)), by=1/12)))
    }
    else{
      dates <- as.character(as.yearmon(seq(from=as.numeric(as.yearmon(from)),length.out=length.out, by=1/12)))
    }
  }
  else if(by == "quarter"){
    if(missing(length.out)){
      dates <- as.character(as.yearqtr(seq(from=as.numeric(as.yearqtr(from)),to=as.numeric(as.yearqtr(to)), by=1/4)))
    }
    else{
      dates <- as.character(as.yearqtr(seq(from=as.numeric(as.yearqtr(from)), length.out=length.out, by=1/4)))
    }
  }
  else if(by == "year"){
    if(missing(length.out)){
      dates <- as.character(seq(from=as.numeric("2000"),to=as.numeric("2005"),by=1))
    }
    else{
      dates <- as.character(seq(from=as.numeric("2000"),length.out=length.out,by=1))
    }
  }
  return(dates)
}
