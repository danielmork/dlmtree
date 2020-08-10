#' PM2.5 Exposure data
#' @description Data.frame containing a sample of weekly average PM2.5 exposures
#' across a range of states/counties. The PM2.5 data was downloaded from US EPA
#' (https://aqs.epa.gov/aqsweb/airdata/download_files.html) daily data summaries
#' and averaged by week. Forty-week ranges were assess for non-missingness
#' and grouped for this dataset.
#'
#' @docType data
#'
#' @usage data(pm25Exposures)
#'
#' @format data.frame; columns: S = state, C = city, 1-40 = weekly exposure data 
#'
#' @keywords datasets
#'
#' @references \url{https://www.epa.gov/outdoor-air-quality-data}
#'
#' @source \url{https://aqs.epa.gov/aqsweb/airdata/download_files.html}
#'
"pm25Exposures"
