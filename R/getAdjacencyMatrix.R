#' get.adjacency.matrix
#' @title Get an adjacency matrix
#' @param ... FIPS coordinates of counties
#'
#' @return An adjacency matrix consisting of 1s and 0s
#' @export
#'
getAdjacencyMatrix <- function(fips, state_abbv = "CO") {
  
  # Colorado county data
  county_data <- counties(state = state_abbv, cb = TRUE)  # Get adjacency information
  
  # Merge the county data with your original data
  counties_with_coordinates <- county_data %>% filter(COUNTYFP %in% fips)
  
  counties_sf <- st_as_sf(counties_with_coordinates)
  counties_list <- st_geometry(counties_sf)
  
  nb <- poly2nb(counties_list)
  
  W <- nb2mat(nb, style = "B")      # Convert the spatial weights matrix to a binary matrix:
  
  return(W)
}
