
#' Delineate a watershed from Canadian Digital Elevation Data within British Columbia
#'
#'This function requires a pour point and extent to delineate a watershed
#'from digital elevation data derived from CDED, obtained using \link[bcmaps]{cded}.
#'Flow tracing is done using \link[RSAGA] and multiple flow directions (MFD). Other
#'transformations are done using gdal (gdalwarp,gdal_traslate & gdal_polygonize) and the \link[raster] package.
#'
#' @param pour_point Geo-referenced sf point object, make sure pour point lies exactly within channel of interest.
#' @param xmin Estimated minimum X coordinate in map units defining extent of basin to delineate.
#' @param xmax Estimated maximum X coordinate in map units defining extent of basin to delineate.
#' @param ymin Estimated minimum Y coordinate in map units defining extent of basin to delineate.
#' @param ymax Estimated maximum Y coordinate in map units defining extent of basin to delineate.
#' @param convergence Convergence parameter for SAGA MFD algorithm.
#' @param minslope Minimum slope parameter for SAGA fill DEM algorithm.
#' @return A list with the delineated basin as a sf polygon, and the CDED raster object for the extent provided.
#' @export
bc_delin_basin<-function(pour_point,xmin,xmax,ymin,ymax,convergence=1.1,minslope=0.1)
{
  #Get extent as sf rectangle
  extent<-sf::st_sfc(sf::st_polygon(x = list(rbind(c(xmin,ymin),c(xmax,ymin),c(xmax,ymax),c(xmin,ymax),c(xmin,ymin)))),crs=sf::st_crs(pour_point))

  #Get cded dem intersecting extent aoi
  cded_dem <- bcmaps::cded(aoi=extent)

  #Project cded to CRS of pour point using gdal (raster too slow)
  gdal_call <- paste0('gdalwarp -t_srs ',sf::st_crs(pour_point)[[1]],' ',cded_dem,' ',file.path(tempdir(),'cded_dem_warp.vrt -r cubic'))
  system(gdal_call)

  #Import the warped cded dem into SAGA-GIS format
  RSAGA::rsaga.import.gdal(in.grid = file.path(tempdir(),'cded_dem_warp.vrt'), out.grid = file.path(tempdir(),'cded_dem'))

  #Fill the cded dem before flow tracing with MFD
  RSAGA::rsaga.geoprocessor(lib='ta_preprocessor',
                            module = 4,
                            param = list(ELEV=file.path(tempdir(),'cded_dem.sgrd'),
                                         FILLED=file.path(tempdir(),'cded_dem_filled.sgrd'),
                                         MINSLOPE=minslope))
  #Compute flow tracing with MFD
  RSAGA::rsaga.geoprocessor(lib='ta_hydrology',
                            module = 4,
                            param = list(TARGET_PT_X=sf::st_coordinates(pour_point)[,'X'],
                                         TARGET_PT_Y=sf::st_coordinates(pour_point)[,'Y'],
                                         ELEVATION=file.path(tempdir(),'cded_dem_filled.sgrd'),
                                         AREA=file.path(tempdir(),'upslope_area.sgrd'),
                                         CONVERGE=convergence))


  #Translate to GeoTIFF format using gdal
  gdal_call<-paste0('gdal_translate ',file.path(tempdir(),'upslope_area.sdat'),' ',file.path(tempdir(),'upslope_area.tif'),sep='')
  system(gdal_call)

  #Read in the MFD raster
  mfd<-raster::raster(file.path(tempdir(),'upslope_area.tif'))

  #Get only cells contributing more than 50% flow, convert to binary
  mfd[mfd<50]<-NA
  mfd<-mfd>=50

  #write out so gdal can access data
  raster::writeRaster(mfd,file.path(tempdir(),'basin_binary.tif'))

  #Convert to polygons using gdal (raster too slow)
  gdal_call<-paste('gdal_polygonize.py ',file.path(tempdir(),'basin_binary.tif'),' ',file.path(tempdir(),'basin_ply.gpkg'))
  system(gdal_call)

  #Read in the basin polygon and union
  basin<-sf::read_sf(file.path(tempdir(),'basin_ply.gpkg'))
  basin<-sf::st_union(basin)

  #CRS is lost somewhere along the line
  sf::st_crs(basin)<-sf::st_crs(pour_point)

  #Return basin polygon, and path to warped cded dem as list
  return(list(basin,file.path(tempdir(),'cded_dem_warp.vrt')))
}
