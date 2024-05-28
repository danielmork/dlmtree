


#' Download simulated data for dlmtree articles
#'
#' @return A data frame with 10000 rows (observations) and 202 variables. All data is simulated. The variables are:
#'   \item{bwgaz}{Outcome to be used. Simulated birth weight for gestational age z-score.}
#'   \item{ChildSex}{Binary sex of child.}
#'   \item{MomAge}{Continuous age in years.}
#'   \item{GestAge}{Continuous estimated gestational age at birth in weeks.}
#'   \item{MomHeightIn}{Continuous maternal height in inches.}
#'   \item{MomPriorWeightLbs}{Continuous mothers pre-pregnancy weight in pounds.}
#'   \item{MomPriorBMI}{Continuous mothers pre-pregnancy BMI.}
#'   \item{race}{Categorical race.}
#'   \item{Hispanic}{Binary indicator of Hispanic.}
#'   \item{MomEdu}{Categorical maternal highest educational attainment.}
#'   \item{SmkAny}{Binary indicator of any smoking during pregnancy.}
#'   \item{Marital}{Categorical maternal marital status.}
#'   \item{Income}{Categorical income.}
#'   \item{EstDateConcept}{Estimated date of conception.}
#'   \item{EstMonthConcept}{Estimated month of conception.}
#'   \item{EstYearConcept}{Estimated year of conception.}
#'   \item{pm25_1 - pm25_37}{Weekly average exposure to PM2.5 for weeks 1 to 37. The columns are already scaled by the exposure IQR of 0.35.}
#'   \item{no2_1 - no2_37}{Weekly average exposure to NO2 for weeks 1 to 37. The columns are already scaled by the exposure IQR of 9.13.}
#'   \item{so2_1 - so2_37}{Weekly average exposure to SO2 for weeks 1 to 37. The columns are already scaled by the exposure IQR of 0.96.}
#'   \item{co2_1 - co2_37}{Weekly average exposure to CO for weeks 1 to 37. The columns are already scaled by the exposure IQR of 0.15.}
#'   \item{temp_1 - temp_37}{Weekly average exposure to temperature for weeks 1 to 37. The columns are already scaled by the exposure IQR of 27.93}
#'   \item{source}{Variable indicating that the data came from the bdlim package.}
#' @export
#'
#' @examples 
#' sbd_dlmtree <- get_sbd_dlmtree()
#' 
#' 
get_sbd_dlmtree <- function(){
  
  sbd_dlmtree <- NULL
  
  temp <- tempfile()
  download.file("https://github.com/danielmork/dlmtree/raw/master/vignettes/articles/sbd_dlmtree.rda", temp )
  load(temp)
  
  return(sbd_dlmtree)
}

