#' @title Final conservation score ex situ
#' @name FCSex
#' @description This function calculates the average of the three ex situ conservation metrics
#'   returning a final conservation score summary table. It also assigns conservation priority
#'   categories
#' @inheritParams GRSex
#' @param Gap_Map logical, if \code{TRUE} the function will calculate gap maps for each species
#'  analyzed and will return a list with three FCSex, GRSex_maps, and ERSex_maps
#' @return A data frame summarizing the ex-situ gap analysis scores:
#'
#' \tabular{lcc}{
#' species \tab Species name \cr
#' SRSex \tab Sampling representativeness score ex situ \cr
#' GRSex \tab Geographical representativeness score ex situ \cr
#' ERSex \tab Ecological representativeness score ex situ \cr
#' FCSex \tab Final conservation score ex situ  \cr
#' FCSex_class \tab The conservation priority category \cr
#' }
#'
#' @examples
#' ##Obtaining occurrences from example
#' data(CucurbitaData)
#' ##Obtaining species names from the data
#' Cucurbita_splist <- unique(CucurbitaData$taxon)
#' ##Obtaining raster_list
#' data(CucurbitaRasters)
#' CucurbitaRasters <- raster::unstack(CucurbitaRasters)
#' ##Obtaining ecoregions shapefile
#' data(ecoregions)
#' # Running all three Ex-situ gap analysis steps using a unique function
#' FCSex_df <- FCSex(Species_list=Cucurbita_splist,
#'                   Occurrence_data=CucurbitaData,
#'                   Raster_list=CucurbitaRasters,
#'                   Buffer_distance=50000,
#'                   Ecoregions_shp=ecoregions)
#'
#'@references
#'
#' Khoury et al. (2019) Ecological Indicators 98:420-429.
#' \url{https://doi.org/10.1016/j.ecolind.2018.11.016}
#'
#' @export
FCSex <- function(Species_list,
                  Occurrence_data,
                  Raster_list,
                  Buffer_distance = 50000,
                  Ecoregions_shp = NULL,
                  Gap_Map = FALSE){

  SRSex_df <- NULL
  GRSex_df <- NULL
  ERSex_df <- NULL
  FCSex_df <- NULL

  #Checking Occurrence_data format
  par_names <- c("taxon","latitude","longitude","type")

  if(missing(Occurrence_data)){
    stop("Please add a valid data frame with columns: taxon, latitude, longitude, type")
  }

  if(identical(names(Occurrence_data),par_names)==FALSE){
    stop("Please format the column names in your dataframe as taxon,latitude,longitude,type")
  }

  # Load in ecoregions shp
  if(is.null(Ecoregions_shp) | missing(Ecoregions_shp)){
    if(file.exists(system.file("data/preloaded_data/ecoRegion/tnc_terr_ecoregions.shp",
                               package = "GapAnalysis"))){
      Ecoregions_shp <- raster::shapefile(system.file("data/preloaded_data/ecoRegion/tnc_terr_ecoregions.shp",
                                                      package = "GapAnalysis"),
                                          encoding = "UTF-8")
    } else {
      stop("Ecoregions file is not available yet. Please run the function GetDatasets() and try again")
    }
  } else{
    Ecoregions_shp <- Ecoregions_shp
  }

  # call SRSex
  SRSex_df <- SRSex(Species_list = Species_list,
                    Occurrence_data = Occurrence_data)
  # call GRSex
  GRSex_df <- GRSex(Occurrence_data = Occurrence_data,
                    Species_list = Species_list,
                    Raster_list = Raster_list,
                    Buffer_distance = Buffer_distance,
                    Gap_Map = Gap_Map)
  # call ERSex
  ERSex_df <- ERSex(Species_list = Species_list,
                    Occurrence_data = Occurrence_data,
                    Raster_list = Raster_list,
                    Buffer_distance = Buffer_distance,
                    Ecoregions_shp=Ecoregions_shp,
                    Gap_Map = Gap_Map)

  # join the dataframes based on species
  if (is.list(GRSex_df)) {
    FCSex_df <- merge(SRSex_df, GRSex_df, by = "species", all.x = TRUE)
  } else {
    FCSex_df <- merge(SRSex_df, GRSex_df$GRSex, by ="species", all.x = TRUE)
  }

  FCSex_df <- merge(FCSex_df, ERSex_df$ERSex, by = "species", all.x = TRUE)

  FCSex_df[is.na(FCSex_df)] <- 0

  # calculate the mean value for each row to determine fcs per species
  FCSex_df$FCSex <- rowMeans(FCSex_df[, c("SRSex", "GRSex", "ERSex")])

  #assign classes (exsitu)
  FCSex_df$FCSex_class <- with(FCSex_df, ifelse(FCSex < 25, "HP",
                                                ifelse(FCSex >= 25 & FCSex < 50, "MP",
                                                       ifelse(FCSex >= 50 & FCSex < 75, "LP",
                                                              "SC"))))

  if (isTRUE(Gap_Map)) {
    FCSex_df <- list(FCSex = FCSex_df,
                     GRSex_maps = GRSex_df$gap_maps,
                     ERSex_maps = ERSex_df$gap_maps)
  }

  return(FCSex_df)

}
