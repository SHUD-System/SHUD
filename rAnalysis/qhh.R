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
library(gridExtra)

# 配置rSHUD工作环境。
shud.env(prjname = 'qhh', inpath = '../input/qhh/', outpath = '../output/qhh.out/')

ia = getArea()
AA = sum(ia)
hist(ia)
fdc(ia/1e6)
ncell =  length(ia)

oid = getOutlets()  # 获得河流出口的序号(index)
oid
riv = rSHUD::readriv()  # 读取河流文件
?rSHUD::readriv
View(riv)
View(riv@river)

lakid = rSHUD::getLakeEleID() #读：湖泊ID

# sp.mesh = rSHUD::sp.mesh2Shape() # gis:读取三角形网
# sp.riv = rSHUD::sp.riv2shp() # gis： 读取河流
# raster::plot(sp.mesh, axes=TRUE, lwd=0.2, border='gray70') #绘图：三角形网
# raster::plot(sp.riv, add=TRUE, col='blue')
# raster::plot(sp.mesh[lakid, ], col='lightblue', add=TRUE)
# grid()


# 河川径流
qdown = rSHUD::readout('rivqdown')
dim(qdown)
xts::plot.xts(qdown[, oid])

# 地下水
igw = readout('eleygw')
dim(igw)
plot(apply(X = igw, MARGIN = 1, FUN = mean))

# 潜在蒸散发
ietp = readout('elevetp')  # 读取：潜在蒸散发
head(ietp[, 1:5]) # 查看数据前5列的头部数据
dim(ietp)  # 数据维度
etp.daily = apply.daily(ietp, sum) /ncell  #时间序列：日均值
etp.monthly = xts::apply.monthly(ietp, sum) /ncell # 时间序列：月均值
plot(etp.daily)
g1 = ggplot2::autoplot(etp.daily*1000) +
  xlab('Time') +ylab('ETP (mm/day)') + theme_bw()
g2 = ggplot2::autoplot(etp.monthly*1000) +
  xlab('Time') +ylab('ETP (mm/month)') + theme_bw()
gg = gridExtra::arrangeGrob(grobs = list(g1, g2))
ggsave(filename = 'etp.png', plot=gg, width = 11, height = 9, units = 'in',dpi = 600)

# 实际蒸散发
ieta = readout('eleveta')  # 读取：实际蒸散发
eta.daily = apply.daily(ieta, sum) /ncell  #时间序列：日均值
eta.monthly = xts::apply.monthly(ieta, sum) /ncell # 时间序列：月均值
plot(etp.daily)


pa = cbind(etp.monthly, eta.monthly)
class(pa)
colnames(pa) = c('ETP', 'ETA')
autoplot(pa*1000, facets = NULL)+
  xlab('Time') +ylab('Evapotranspiration (mm/month)') + theme_bw()

# 对比陆地与湖泊实际蒸散发量，月均值
nlake =length(lakid)
et.vs = cbind(apply.monthly(ieta[, -1 * lakid], sum)/(ncell - nlake),
              apply.monthly(ieta[, lakid], sum)/nlake )
colnames(et.vs) = c('Land', 'Lake')
autoplot(et.vs*1000, facets = NULL)+
  xlab('Time') +ylab('Evapotranspiration (mm/month)') + theme_bw()


# Hydrograph
qdown = rSHUD::readout('rivqdown')
prcp = rSHUD::readout('elevprcp')
q = apply.daily(qdown[, oid], sum)
p = apply.daily(prcp, sum)
pq = cbind(p, q)
colnames(pq) = c('Precipitation', 'Discharge')
rSHUD::hydrograph(pq[1:(730) + (6*365), ])


# 水量平衡
xl = loaddata()
wb = wb.all(xl = xl)

