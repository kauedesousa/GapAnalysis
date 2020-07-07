#' @title Final conservation score in situ
#' @name FCSin
#' @description This function calculates the average of the three in situ conservation metrics and
#' assigns a priority category based on the results
#' @inheritParams GRSex
#' @inheritParams GRSin
#' @param Gap_Map logical, if \code{TRUE} the function will calculate gap maps for each species
#'  analyzed and will return a list with four slots FCSin, SRSin_maps,GRSin_maps,and ERSin_maps
#'
#' @return A data frame summarizing the in situ gap analysis scores:
#' \tabular{lcc}{
#' species \tab Species name \cr
#' SRSin \tab Sampling representativeness score in situ  \cr
#' GRSin \tab Geographical representativeness score in situ \cr
#' ERSin \tab Ecological representativeness score in situ \cr
#' FCSin \tab Final conservation score in situ  \cr
#' }
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
#' ##Obtaining Ecoregions_shpions shapefile
#' data(ecoregions)
#' #Running all three In-situ gap analysis steps using a unique function
#' FCSin_df <- FCSin(Species_list=Cucurbita_splist,
#'                                       Occurrence_data=CucurbitaData,
#'                                       Raster_list=CucurbitaRasters,
#'                                       Ecoregions_shp=ecoregions,
#'                                       Pro_areas=ProtectedAreas)
#'
#' @references
#'
#' Khoury et al. (2019) Diversity and Distributions 26(2):209-225. doi: 10.1111/DDI.13008
#'
#' @export

#' @importFrom raster overlay crop raster extent

FCSin <- function(Species_list,
                  Occurrence_data,
                  Raster_list,
                  Ecoregions_shp = NULL,
                  Pro_areas = NULL,
                  Gap_Map = FALSE) {

  SRSin_df <- NULL
  GRSin_df <- NULL
  ERSin_df <- NULL
  FCSIn_df <- NULL

  #Checking Occurrence_data format
  par_names <- c("taxon","latitude","longitude","type")


  if(missing(Occurrence_data)){
    stop("Please add a valid data frame with columns: taxon,latitude,longitude,type")
  }

  if(identical(names(Occurrence_data),par_names)==FALSE){
    stop("Please format the column names in your dataframe as taxon,latitude,longitude,type")
  }

  # load in protected area raster
  if(is.null(Pro_areas)){
    if(file.exists(system.file("data/preloaded_data/protectedArea/wdpa_reclass.tif",
                               package = "GapAnalysis"))){
      Pro_areas <- raster::raster(system.file("data/preloaded_data/protectedArea/wdpa_reclass.tif",
                                              package = "GapAnalysis"))
    } else {
      stop("Protected areas file is not available yet. Please run the function GetDatasets()  and try again")
    }
  }

  # Load in ecoregions shp
  if (is.null(Ecoregions_shp)) {
    if(file.exists(system.file("data/preloaded_data/ecoRegion/tnc_terr_ecoregions.shp",
                               package = "GapAnalysis"))){
      Ecoregions_shp <- raster::shapefile(system.file("data/preloaded_data/ecoRegion/tnc_terr_ecoregions.shp",
                                                      package = "GapAnalysis"),encoding = "UTF-8")
    } else {
      stop("Ecoregions file is not available yet. Please run the function GetDatasets() and try again")
    }
  }

  # call SRSin
  SRSin_df <- SRSin(Species_list = Species_list,
                    Occurrence_data = Occurrence_data,
                    Raster_list = Raster_list,
                    Pro_areas = Pro_areas,
                    Gap_Map = Gap_Map)

  GRSin_df <- GRSin(Species_list = Species_list,
                    Occurrence_data = Occurrence_data,
                    Raster_list = Raster_list,
                    Pro_areas = Pro_areas,
                    Gap_Map = Gap_Map)

  ERSin_df <- ERSin(Species_list = Species_list,
                    Occurrence_data =Occurrence_data,
                    Raster_list = Raster_list,
                    Pro_areas = Pro_areas,
                    Ecoregions_shp = Ecoregions_shp,
                    Gap_Map = Gap_Map)


  if(isFALSE(Gap_Map)){
    # join the dataframes based on species
    FCSin_df <- merge(SRSin_df, GRSin_df, by ="species")
    FCSin_df <- merge(FCSin_df, ERSin_df, by = "species")
  }

  if(isTRUE(Gap_Map)){

    FCSin_df <- merge(SRSin_df$SRSin, GRSin_df$GRSin, by ="species")
    FCSin_df <- merge(FCSin_df, ERSin_df$ERSin, by = "species")

  }

  # calculate the mean value for each row to determine fcs per species
  FCSin_df$FCSin <- rowMeans(FCSin_df[,c("SRSin", "GRSin", "ERSin")])

  FCSin_df$FCSin_class <- with(FCSin_df, ifelse(FCSin < 25, "HP",
                                                ifelse(FCSin >= 25 & FCSin < 50, "MP",
                                                       ifelse(FCSin >= 50 & FCSin < 75, "LP",
                                                              "SC"))))

  if (isTRUE(Gap_Map)) {
    FCSin_df <- list(FCSin = FCSin_df,
                     SRSin_maps = SRSin_df$gap_maps,
                     GRSin_maps = GRSin_df$GapMapIn_list,
                     ERSin_maps = ERSin_df$gap_maps)
  }

  return(FCSin_df)

}
