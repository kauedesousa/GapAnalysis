#' @title Geographical representativeness score in situ
#' @name GRSin
#' @description The GRSin process provides a geographic measurement of the proportion of a species’ range that can be considered
#' to be conserved in protected areas. The GRSin compares the area of the distribution model located within protected areas versus
#' the total area of the model, considering comprehensive conservation to have been accomplished only when the entire distribution
#' occurs within protected areas.
#'
#' @inheritParams GRSex
#' @param Pro_areas A raster object of protected areas, with values 0 and 1, where 1 indicates that a given grid-cell is located within
#'  a protected area. If \code{Pro_areas = NULL} the function will use a protected area raster
#'  file provided for your use after run \code{GetDatasets()}
#' @param Gap_Map logical, if \code{TRUE} the function will calculate gap maps for each species analyzed and will return a list
#'  with two slots ERSin and gap_maps
#'
#' @return This function returns a data frame with two columns:
#'
#' \tabular{lcc}{
#' species \tab Species name \cr
#' GRSin \tab GRSin value calculated\cr
#' }
#'
#' @examples
#' ##Obtaining occurrences from example
#' data(CucurbitaData)
#' ##Obtaining species names from the data
#' Cucurbita_splist <- unique(CucurbitaData$taxon)
#' ##Obtaining Raster_list
#' data(CucurbitaRasters)
#' CucurbitaRasters <- raster::unstack(CucurbitaRasters)
#' ##Obtaining protected areas raster
#' data(ProtectedAreas)
#' #Running GRSin
#' GRSin_df <- GRSin(Species_list = Cucurbita_splist,
#'                     Occurrence_data = CucurbitaData,
#'                     Raster_list = CucurbitaRasters,
#'                     Pro_areas=ProtectedAreas)
#'
#'
#'@references
#'
#' Khoury et al. (2019) Ecological Indicators 98:420-429. doi: 10.1016/j.ecolind.2018.11.016
#'
#' @export
#' @importFrom stats median
#' @importFrom raster raster crop area
GRSin <- function(Species_list,
                  Occurrence_data,
                  Raster_list,
                  Pro_areas=NULL,
                  Gap_Map=FALSE){

# suppressMessages(require(rgdal))
# suppressMessages(require(raster))
# suppressMessages(require(tmap))
# suppressMessages(require(fasterize))
# suppressMessages(require(sf))

#importFrom("methods", "as")
#importFrom("stats", "complete.cases", "filter", "median")
#importFrom("utils", "data", "memory.limit", "read.csv", "write.csv")

  #Checking Occurrence_data format
  par_names <- c("taxon","latitude","longitude","type")


  if(missing(Occurrence_data)){
    stop("Please add a valid data frame with columns: taxon,latitude,longitude,type")
  }


  if(identical(names(Occurrence_data),par_names)==FALSE){
    stop("Please format the column names in your dataframe as taxon,latitude,longitude,type")
  }

  #Checking if user is using a raster list or a raster stack
  if (isTRUE("RasterStack" %in% class(Raster_list))) {
    Raster_list <- raster::unstack(Raster_list)
  } else {
    Raster_list <- Raster_list
  }

  df <- data.frame(matrix(ncol=2, nrow = length(Species_list)))

  colnames(df) <- c("species", "GRSin")
  # load in protected areas raster
  if (is.null(Pro_areas)) {
    if(file.exists(system.file("data/preloaded_data/protectedArea/wdpa_reclass.tif",package = "GapAnalysis"))){
      Pro_areas <- raster::raster(system.file("data/preloaded_data/protectedArea/wdpa_reclass.tif",package = "GapAnalysis"))
    } else {
      stop("Protected areas file is not available yet. Please run the function GetDatasets()  and try again")
    }
  } else{
    Pro_areas <- Pro_areas
  }

  if(Gap_Map==TRUE){
    GapMapIn_list <- list()
  }


  # loop over species list
  for(i in seq_along(Species_list)){

    # select raster with species name
    # this assumes that the user provided the rasters in the same order as the Species_list
    # as stated in the documentation
    sdm <- Raster_list[[i]]

    # determine the area of predicted presence of a species based on the species distribution map
    sdm1 <- sdm
    Pro_areas1 <- raster::crop(x = Pro_areas,y = sdm1)
    sdm1[sdm1[] != 1] <- NA

    if(raster::res(Pro_areas1)[1] != raster::res(sdm)[1]){
      Pro_areas1 <- raster::resample(x = Pro_areas1, y = sdm)
    }
    cell_size <- raster::area(sdm1, na.rm=TRUE, weights=FALSE)
    cell_size <- cell_size[!is.na(cell_size)]
    thrshold_area <- length(cell_size)*median(cell_size)

    # mask the protected area raster to the species distribution map and calculate area
    Pro_areas1[Pro_areas1[] == 0] <-NA
    Pro_areas1 <- Pro_areas1 * sdm1
    # calculate area
    cell_size <- raster::area(Pro_areas1, na.rm=TRUE, weights=FALSE)
    cell_size <- cell_size[!is.na(cell_size)]
    protected_area <- length(cell_size)*stats::median(cell_size)
    if(!is.na(protected_area)){
      # calculate GRSin
      GRSin <- min(c(100, protected_area/thrshold_area*100))
      df$species[i] <- as.character(Species_list[i])
      df$GRSin[i] <- GRSin
    }else{
      df$species[i] <- as.character(Species_list[i])
      df$GRSin[i] <- 0
    }
    #GRSin gap map

    if(Gap_Map==TRUE){
      message(paste0("Calculating GRSin gap map for ",as.character(Species_list[i])),"\n")
      Pro_areas1[is.na(Pro_areas1),] <-  0
      gap_map <- sdm - Pro_areas1
      gap_map[gap_map == 0,] <- NA
      GapMapIn_list[[i]] <- gap_map
      names(GapMapIn_list[[i]] ) <- Species_list[[i]]
    }
  }
  if(Gap_Map==TRUE){
    df <- list(GRSin= df,gap_maps=GapMapIn_list)
  } else {
    df <- df
  }
  return(df)
}
