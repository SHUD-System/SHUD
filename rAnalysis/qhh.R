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

#' ===========CACHE CREEK ====================
fn = '/Volumes/home/Project_backup/SAC/forcing/csv2000-2017/x19y114.csv'
x=readLines(fn)
tb=read.table(text=x[-1], header = TRUE, sep = ' ')
tsd = as.xts(tb[, -1], order.by = as.POSIXct('2000-01-01 00:00')+tb[, 1]*86400)
tsd = cbind(tsd[, 1]*1000, tsd[, 2],
            tsd[, 3],tsd[,4]/86400, tsd[,5]/86400)
tt = time(tsd)
x = tsd[tt >= as.POSIXct('2000-01-01') & tt <= as.POSIXct('2010-01-01')]
x = round(x, 4)
fn = '/Users/leleshu/CloudDrive/SHUD/github/SHUD/input/ccw/forcing.csv'
write.tsd(x, fn)



#' ===========CACHE CREEK ====================
fn = '/Users/leleshu/CloudDrive/SHUD/github/SHUD/input/ccw/forcing.csv'
x=read.tsd(file =fn)[[1]]
x=round(x, 3)
plot.zoo(x)
write.tsd(x, fn)
