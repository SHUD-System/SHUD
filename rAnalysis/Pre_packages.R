#' ===================================================
#' Required R package installation.
#' Author: Lele Shu <shulele@lzb.ac.cn>
#' Website: www.shud.xyz
#' Date: 2022-10-10
#' Copyright: MIT-Licence
#' ===================================================

libs=c('raster', 'sp', 'rgdal', 'rgeos', 'rasterVis', 'rgl',
       'gstat', 'zoo', 'xts', 'hydroGOF', 'hydroTSM',
       'ggplot2', 'lattice', 'reshape2',
       'abind', 'spam', 'fields',
       'doParallel', 'foreach', 'doMC',
       'RTriangle', 'geometry',
       'Rcpp', 'sf',
       'deldir', 'proj4', 'lubridate', 'interp', 'ncdf4', 'GGally',
       'blogdown', 'roxygen2',
       'mapproj', 'coord_map', 'akima', 'graticule',
       'magick', 'pdftools'
)
# libs=c('lavaan', 'multcompView', 'pander', 'psych', 'pwr', 'rmdpartials', 'ufs', 'kableExtra') # ROSETTA

fun.inst <- function(libs){
  nx=length(libs)
  i=1
  for(i in 1:nx){
    x=libs[i]
    message(i, '/', nx, '\t', x)
    if(!require(x, character.only = TRUE, )){
      message('Installing ', x)
      install.packages(x)
    }else{
      message('Existing ', x)
    }
  }
}
fun.inst(libs)


if(!require(devtools)){
  install.packages("devtools")
}
devtools::install_github("davidcsterratt/RTriangle", subdir="pkg")
