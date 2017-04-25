library(tidyverse)

nao = read.fwf(
  'http://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/norm.nao.monthly.b5001.current.ascii',
  c(5, 5, 13), col.names=c('year','month', 'value'))
nao
summary(nao)

# 2009 - 2015
qplot(year, value, data=nao)

qplot(year, value, data=nao %>% filter(year %in% 2009:2015))


