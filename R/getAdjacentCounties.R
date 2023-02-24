#' get.adjacent.counties
#' @title Finding an adjacent county
#' @param ... NA
#'
#' @return An adjacent county to a provided county with FIPS coordinate
#' @export
#'

getAdjacentCounties <- function(fips, state_abbv = "CO") {
  
  # Load the counties shapefile
  counties <- counties(state = state_abbv, cb = TRUE)
  
  # Get the row for the county with the specified FIPS code
  county_row <- counties %>% filter(COUNTYFP == fips)
  
  # Get the polygon for the county
  county_polygon <- county_row$geometry
  
  # Get the rows for all counties that intersect the county polygon
  adjacent_counties <- counties %>%
    st_intersects(county_polygon) %>%
    as.data.frame()
  
  # Return the FIPS codes for the adjacent counties
  adj_county_fips = fips
  
  while(adj_county_fips == fips){
    adj_county_idx <- sample(adjacent_counties$row.id, 1)
    adj_county_fips <- counties$COUNTYFP[adj_county_idx]
  }
  
  return(adj_county_fips)
}