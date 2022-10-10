#' ===================================================
#' Analysis the output results from SHUD model
#' Author: Lele Shu <shulele@lzb.ac.cn>
#' Website: www.shud.xyz
#' Date: 2022-10-10
#' Copyright: MIT-Licence
#' ===================================================
#' Function:
#'  The analysis for simulation of Heihe River, China
#'  Simulation period: 19890101 - 



rm(list=ls())
library(rSHUD)
library(xts)
library(raster)
library(rgdal)
library(sp)
library(ggplot2)


shud.env(prjname = 'heihe', inpath = '../input/heihe/', outpath = '../output/heihe.out/')
AA = getArea()
oid = getOutlets()
qo = ts2Daily(read.tsd('../input/heihe/heihe.tsd.obs')[[1]])
qdown = readout('rivqdown')
qs = ts2Daily(qdown[, oid])
ts = time(qs)
to = time(qo)
ct= ts[to %in% ts]
qq = cbind('SIM' = qs, 'OBS' = qo/86400)[ts,]

autoplot(qq)

loaddata()
