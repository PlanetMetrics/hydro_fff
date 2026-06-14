#' Summary information about spring climate
#'
#' computes summary information about spring temperature
#' @param clim_data  data frame with columns tavg (C)
#'	 year, month (integer), day
#' @param spring_months (vector of integers) to include in spring; default 4,5,6
#' @param spring_out (default FALSE) set to TRUE to output average spring  temperature for all years
#' @return returns list of extremes
#' \describe{
#'  \item{Tavg_spring}{mean_springT mean spring temperature (C)}
#'  \item{Tmax_spring}{  highest spring temperature  (C)}
#'  \item{Tmin_spring}{  lowest spring temperature  (C)}
#'  \item{warmest_spring}{ year with highest spring temperature  (year)}
#'  \item{coldest_spring}{ year with lowest spring temperature  (year)}
#'  and array of spring temperature means for all years if spring_out set to TRUE
#'  }

spring_summary = function(clim_data, spring_months = c(4:6), spring_out=FALSE) {

  library(tidyverse)
  
  # check to make sure all colunns available
  stopifnot(all(c("month", "year", "tavg") %in% names(clim_data)))
  
  spring = clim_data %>% subset(month %in% spring_months)

  S_means_all = spring %>% group_by(year) %>% select(year, tavg) %>% summarize_all(list(mean=mean))


  S_extremes = spring %>% summarize( Tavg_spring = mean(tavg),
                                          Tmax_spring = max(tavg),
                                          Tmin_spring=min(tavg),

                                          warmest_spring = year[which.max(tavg)],
                                          coldest_spring = year[which.min(tavg)]
                                          )


  if (spring_out)
  return(list(as.list(S_extremes), S_means_all[,c("year","mean")]))
  else
  return(as.list(S_extremes))

}
