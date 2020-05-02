#' @title Geographical representativeness score estimation (Ex-situ conservation)
#' @name GRSex
#' @description This function performs an estimation of the geographical representativeness score for ex-situ gap analysis (GRSex)
#' using Ramirez-Villegas et al., (2010) methodology. GRS ex-situ score is calculated as:
#'
#' \deqn{GRSex = min(100,(Masked Area of Buffered G Occurrences / Total Area of Predicted Habitat)*100)}
#'
#' @param occurrenceData A data frame object with the species name, geographical coordinates, and type of records (G or H) for a given species
#' @param species_list An species list to calculate the GRSex metrics.
#' @param raster_list A list representing the species distribution models for the species list provided loaded in raster format. This list must match the same order of the species list.
#' @param bufferDistance Geographical distance used to create circular buffers around germplasm. Default: 50000 that is 50 km around germplasm accessions (CA50)
#'
#' @return This function returns a data frame with two columns:
#'
#' \tabular{lcc}{
#' species \tab Species name \cr
#' GRSex \tab GRSex value calculated\cr
#' }
#'
#' @examples
#' ##Obtaining occurrences from example
#' data(CucurbitaData)
#' speciesList <- unique(CucurbitaData$taxon)
#' ## Obtaining rasterList objet. ##
#' data(CucurbitaRasters)
#' #Calculating GRSex value
#' GRSex_df <- GRSex(species_list = speciesList,
#'                     occurrenceData = CucurbitaData,
#'                     raster_list = CucurbitaRasters,
#'                     bufferDistance = 50000)
#'
#' @references
#' Ramirez-Villegas, J., Khoury, C., Jarvis, A., Debouck, D. G., & Guarino, L. (2010).
#' A Gap Analysis Methodology for Collecting Crop Genepools: A Case Study with Phaseolus Beans.
#' PLOS ONE, 5(10), e13497. Retrieved from https://doi.org/10.1371/journal.pone.0013497
#'
#'Ramirez-Villegas, J., Khoury, C., Jarvis, A., Debouck, D. G., & Guarino, L. (2010).
#'A Gap Analysis Methodology for Collecting Crop Genepools: A Case Study with Phaseolus Beans.
#'PLOS ONE, 5(10), e13497. Retrieved from https://doi.org/10.1371/journal.pone.0013497
#'
#' Khoury, C. K., Amariles, D., Soto, J. S., Diaz, M. V., Sotelo, S., Sosa, C. C., … Jarvis, A. (2019).
#' Comprehensiveness of conservation of useful wild plants: An operational indicator for biodiversity
#' and sustainable development targets. Ecological Indicators. https://doi.org/10.1016/j.ecolind.2018.11.016
#'
#'
#'Hijmans, R.J. and Spooner, D.M. (2001).
#'Geographic distribution of wild potato species. Am. J. Bot., 88: 2101-2112. doi:10.2307/3558435
#'
#' @export
#' @import fasterize sp
#' @importFrom tidyr drop_na
#' @importFrom magrittr "%>%"
#' @importFrom dplyr filter select
#' @importFrom fasterize fasterize
#' @importFrom stats median


GRSex <- function(occurrenceData, species_list, raster_list, bufferDistance) {

  longitude <- NULL
  taxon <- NULL
  type <- NULL
  latitude <-NULL
  # suppressMessages(require(rgdal))
  # suppressMessages(require(raster))

  #importFrom("methods", "as")
  #importFrom("stats", "complete.cases", "filter", "median")
  #importFrom("utils", "data", "memory.limit", "read.csv", "write.csv")
  if(missing(bufferDistance)){
    bufferDistance <- 50000
  }
  # create a dataframe to hold the components
  df <- data.frame(matrix(ncol = 2, nrow = length(species_list)))
  colnames(df) <- c("species", "GRSex")

  for(i in 1:length(sort(species_list))){
    # select species G occurrences
    occData <- occurrenceData %>%
      tidyr::drop_na(longitude)%>%
      dplyr::filter(taxon == species_list[i]) %>%
      dplyr::filter(type == "G")%>%
      dplyr::select(longitude,latitude)

    sp::coordinates(occData) <- ~longitude+latitude
    sp::proj4string(occData) <- sp::CRS("+proj=longlat +datum=WGS84")
    # select raster with species name
    for(j in 1:length(raster_list)){
      if(grepl(j, i, ignore.case = TRUE)){
        sdm <- raster_list[[j]]
      }
    }
    # convert SDM from binary to 1-NA for mask and area
    sdmMask <- sdm
    sdmMask[sdmMask == 0] <- NA
    # buffer G points
#    buffer <- geobuffer::geobuffer_pts(xy = occData,
    buffer <- gapAnalysisR::geobuffer_pts(xy = occData,
                                       dist_m = bufferDistance,
                                       output = 'sf')

    # rasterizing and making it into a mask
    buffer_rs <- fasterize::fasterize(buffer, sdm)
    buffer_rs[!is.na(buffer_rs[])] <- 1
    buffer_rs <- buffer_rs * sdmMask
    # calculate area of buffer
    cell_size<-raster::area(buffer_rs, na.rm=TRUE, weights=FALSE)
    cell_size<-cell_size[!is.na(cell_size)]
    gBufferRas_area<-length(cell_size)*median(cell_size)

    # calculate area of the threshold model
    cell_size<- raster::area(sdmMask, na.rm=TRUE, weights=FALSE)
    cell_size<- cell_size[!is.na(cell_size)]
    pa_spp_area <- length(cell_size)*median(cell_size)
    # calculate GRSex
    grs <- min(c(100, gBufferRas_area/pa_spp_area*100))

    df$species[i] <- as.character(species_list[i])
    df$GRSex[i] <- grs
  }

  return(df)
}