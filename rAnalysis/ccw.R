#' ===================================================
#' Analysis the output results from SHUD model
#' Author: Lele Shu <shulele@lzb.ac.cn>
#' Website: www.shud.xyz
#' Date: 2022-10-10
#' Copyright: MIT-Licence
#' ===================================================
#' Function:
#'
#'


rm(list=ls())
library(rSHUD)
library(xts)
library(raster)
library(rgdal)
library(sp)
library(ggplot2)


shud.env(prjname = 'ccw', inpath = '../input/ccw/', outpath = '../output/ccw.out/')
AA = getArea()
oid = getOutlets()
qo = ts2Daily(read.tsd('../input/ccw/ccw.tsd.obs')[[1]])
qdown = readout('rivqdown')
qs = ts2Daily(qdown[, oid])
ts = time(qs)
to = time(qo)
ct= ts[to %in% ts]
qq = cbind('SIM' = qs, 'OBS' = qo/86400)[ts,]

autoplot(qq[,1])

xl = loaddata()
wb = wb.all(xl=xl)
