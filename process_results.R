# load libraries ----
library(tidyverse) # install.packages('tidyverse')
library(stringr)
library(rgdal)
library(raster)
library(rasterVis)
library(maps)
library(mapproj)
select = dplyr::select
stack  = raster::stack

# define functions ----
# makes process_singledir a function with the parameters inside ()
process_singledir = function(dir_results, dir_simulation, do_csv=T, do_tif=T, do_png=T){
  # dir_results    = 'G:/Team_Folders/Steph/bsb_2015/2_2_15_FM_bsb_50day_results'
  # dir_simulation = 'G:/Team_Folders/Steph/bsb_2015/2_2_15_FM_bsb_50day_simulation'
  
  run = str_replace(basename(dir_results), '_results', '') # run is now equal to a tring that removes all of the path in dir_results and replaces _results with ' ' so that it returns 2_2_15_FM_bsb_50day 
  
  # # conn_lns (connectivity lines) reads the geodatabase called 'output.gdb' for the 'connectivity feature line class
  conn_lns = readOGR(file.path(dir_results, 'output.gdb'), 'Connectivity', verbose=F)
  
  # aggregate across all ToPatchIDs to Gray's Reef (n=4)
  conn_tbl = conn_lns@data %>%
    as_tibble() %>%    
    group_by(FromPatchID) %>%
    summarize(
      quantity = sum(Quantity)) %>%
    ungroup() %>%
    mutate(
      percent = quantity / sum(quantity) * 100) %>%
    arrange(desc(percent))
  
  # write to csv
  if(do_csv){
    write_csv(conn_tbl, sprintf('%s/connectivity.csv', dir_results))
  }

  # get patch id raster, and determine which cells are NA
  r_id = raster(sprintf('%s/PatchData/patch_ids', dir_simulation)) # plot(r_id)
  id_NA = !getValues(r_id) %in% conn_tbl$FromPatchID #look for values in r_id that were not in the conn_tbl (because they do not have connectivity)
  
  # create rasters for quantity and percent ##run through the loop twice, once for v = 'quantity' and the second v = 'percent'
  for (v in c('quantity','percent')){  
    
    # reclassify from patch id to value
    r = reclassify(r_id, conn_tbl[,c('FromPatchID', v)])
    
    # set patch ids without a value to NA
    r[id_NA] = NA
    
    # write to GeoTIFF
    if(do_tif){
      writeRaster(r, sprintf('%s/%s.tif', dir_results, v), overwrite=T)
    }
    
    
    # plot to PNG for easy preview
    if (do_png){
      png(sprintf('%s/%s.png', dir_results, v)) #sprintf = print a string in specific format 
        p = levelplot(r, par.settings=viridisTheme, main=sprintf('%s %s', run, v))
        print(p)
      dev.off()  
    }
  }
}


process_sppyr_dirs = function(dir_sppyr, ...){
  # process all model runs for given species & year

  dirs_results = list.files(dir_sppyr, '.*_results$', full.names=T)
  for (i in 1:length(dirs_results)){
    
    dir_results = dirs_results[i]
    dir_simulation = str_replace(dir_results, '_results', '_simulation') #makes a directory in the same way as before but '_simulation' insteat of '_results'
    cat(sprintf('%03d of %d: %s\n', i, length(dirs_results), basename(dir_results))) #%03d prints the number with leading zeros when the argument is less than three digits long
    
    # process from geodatabase to results csv, tifs, pngs
    process_singledir(dir_results, dir_simulation, ...)
    
  }
}


summarize_sppyr = function(dir_sppyr){

  dirs_results = list.files(dir_sppyr, '.*_results$', full.names=T)
  rasters_quantity = sprintf('%s/quantity.tif', dirs_results)
  stack_quantity = stack(rasters_quantity) #gathers all the rasters_quantity together

  r_mean = mean(stack_quantity, na.rm=T)
  r_sd = calc(stack_quantity, fun=function(x) sd(x, na.rm=T))
  r_cv = r_sd / r_mean * 100
  r_sum = sum(stack_quantity, na.rm=T) # %>% mask(r_mask)

  for (v in c('mean','cv','sum')){

    r = get(sprintf('r_%s',v))

    # write to GeoTIFF
    writeRaster(r, sprintf('%s/%s.tif', dir_sppyr, v), overwrite=T)

    # plot to PNG for easy preview
    png(sprintf('%s/%s.png', dir_sppyr, v))
    p = levelplot(r, par.settings=viridisTheme, main=sprintf('%s %s', basename(dir_sppyr), v))
    print(p)
    dev.off()

  }
}


summarize_spp = function(dir_root, sp){
  # given top-level directory and species code, eg "sp" or "rs" or "bsb",
  # summarize sp_yr/mean.tif across years as sp/mean.tif and sp/cv.tif,
  # ie average dispersal across year means and variation across year means
  # dir_root = 'G:/Team_Folders/Steph'; sp='bsb'

  # dirs_results = list.files(dir_root, sprintf('%s_[0-9]{4}$', sp), full.names=T)
  dirs_results = list.files(dir_root, sprintf('%s_[0-9]{4}$', sp), full.names=T)
  rasters_mean = sprintf('%s/mean.tif', dirs_results)
  stack_mean   = stack(rasters_mean)
  dir_sp = file.path(dir_root, sp)

  if (!file.exists(dir_sp)) dir.create(dir_sp) ##creates the folder with species name (if it doesn't already exist)

  r_mean = mean(stack_mean, na.rm=T)
  r_sd = calc(stack_mean, fun=function(x) sd(x, na.rm=T))
  r_cv = r_sd / r_mean * 100

  for (v in c('mean','cv','sd')){

    r = get(sprintf('r_%s',v))

    # write to GeoTIFF
    writeRaster(r, sprintf('%s/%s.tif', dir_sp, v), overwrite=T)

    # plot to PNG for easy preview
    png(sprintf('%s/%s.png', dir_sp, v))
    p = levelplot(r, par.settings=viridisTheme, main=sprintf('%s %s', basename(dir_sp), v))
    print(p)
    dev.off()

  }
}

##for the sum instead of mean
# summarize_spp = function(dir_root, sp){
#    # given top-level directory and species code, eg "sp" or "rs" or "bsb",
#    # summarize sp_yr/mean.tif across years as sp/mean.tif and sp/cv.tif,
#    # ie average dispersal across year means and variation across year means
#    # dir_root = 'G:/Team_Folders/Steph'; sp='bsb'
# 
#    dirs_results = list.files(dir_root, sprintf('%s_[0-9]{4}$', sp), full.names=T)
#    rasters_sum = sprintf('%s/sum.tif', dirs_results)
#    stack_sum   = stack(rasters_sum)
#    dir_sp = file.path(dir_root, sp)
# 
#    if (!file.exists(dir_sp)) dir.create(dir_sp) ##creates the folder with species name (if it doesn't already exist)
#    
# 
#    r_mean = mean(stack_sum, na.rm=T)
#    r_sd = calc(stack_sum, fun=function(x) sd(x, na.rm=T))
#    r_cv = r_sd / r_mean * 100
# 
#    for (v in c('mean','cv','sd')){
# 
#      r = get(sprintf('r_%s',v))
# 
#      # write to GeoTIFF
#      writeRaster(r, sprintf('%s/%s_sum.tif', dir_sp, v), overwrite=F)
# 
#      # plot to PNG for easy preview
#      png(sprintf('%s/%s_sum.png', dir_sp, v))
#      p = levelplot(r, par.settings=viridisTheme, main=sprintf('%s %s', basename(dir_sp), v))
#      print(p)
#      dev.off()
# 
#    }
#  }


#summarize_spp('G:/Team_Folders/Steph', sp='bsb')
for (sp in c('bsb','gg','rs','sp')){
  summarize_spp('G:/Team_Folders/Steph', sp)  
}



summarize_spp = function(dir_root='G:/Team_Folders/Steph', spp=c('bsb','gg','rs','sp')){
   # given top-level directory and species code, eg "sp" or "rs" or "bsb",
   # summarize sp_yr/mean.tif across years as sp/mean.tif and sp/cv.tif,
   # ie average dispersal across year means and variation across year means
   # dir_root = 'G:/Team_Folders/Steph'; sp='bsb'

   dirs_results = file.path(dir_root, spp)
   rasters_mean = sprintf('%s/mean.tif', dirs_results)
   stack_mean   = stack(rasters_mean)
   dir_spp = file.path(dir_root, '_allspp')

   if (!file.exists(dir_spp)) dir.create(dir_spp)

   r_mean = mean(stack_mean, na.rm=T)
   r_sd = calc(stack_mean, fun=function(x) sd(x, na.rm=T))
   r_cv = r_sd / r_mean * 100

   for (v in c('mean','cv')){

     r = get(sprintf('r_%s',v))

     # write to GeoTIFF
     writeRaster(r, sprintf('%s/%s.tif', dir_spp, v), overwrite=T)

     # plot to PNG for easy preview
     png(sprintf('%s/%s.png', dir_spp, v))
     p = levelplot(r, par.settings=viridisTheme, main=sprintf('%s %s', basename(dir_spp), v))
     print(p)
     dev.off()

   }
}

# # ##for the sum instead of mean
# summarize_spp = function(dir_root='G:/Team_Folders/Steph', spp=c('bsb','gg','rs','sp')){
#   # given top-level directory and species code, eg "sp" or "rs" or "bsb",
#   # summarize sp_yr/mean.tif across years as sp/mean.tif and sp/cv.tif,
#   # ie average dispersal across year means and variation across year means
#   # dir_root = 'G:/Team_Folders/Steph'; sp='bsb'
#   
#   dirs_results = file.path(dir_root, spp)
#   rasters_mean = sprintf('%s/mean_sum.tif', dirs_results)
#   stack_mean   = stack(rasters_mean)
#   dir_spp = file.path(dir_root, '_allspp')
#   
#   if (!file.exists(dir_spp)) dir.create(dir_spp)
#   
#   r_mean = mean(stack_mean, na.rm=T)
#   r_sd = calc(stack_mean, fun=function(x) sd(x, na.rm=T))
#   r_cv = r_sd / r_mean * 100
#   
#   for (v in c('mean','cv')){
#     
#     r = get(sprintf('r_%s',v))
#     
#     # write to GeoTIFF
#     writeRaster(r, sprintf('%s/%s_sum.tif', dir_spp, v), overwrite=T)
#     
#     # plot to PNG for easy preview
#     png(sprintf('%s/%s_sum.png', dir_spp, v))
#     p = levelplot(r, par.settings=viridisTheme, main=sprintf('%s %s', basename(dir_spp), v))
#     print(p)
#     dev.off()  
#     
#   }
# }

summarize_spp(dir_root='G:/Team_Folders/Steph', spp=c('bsb','gg','rs','sp'))




##area maps----

library(tidyverse)
library(raster)
library(plotly)

r = raster('G:/Team_Folders/Steph/_allspp/mean_sum_noseagrass.tif')

d = data_frame(
  quantity = raster::getValues(r),
  cellid   = 1:length(quantity),
  area_km2 = 64)

d2 = d %>%
  filter(!is.na(quantity)) %>%
  arrange(desc(quantity)) %>%
  mutate(
    pct_quantity     = quantity/sum(quantity)*100,
    cum_pct_quantity = cumsum(quantity/sum(quantity)*100),
    cum_area_km2     = cumsum(area_km2))
tail(d2)# 7208 km2
tail(d2$cum_area_km2, 1) # 7208 km2

d3 = d %>%
  left_join(d2, by='cellid')
summary(d3)

r2 = setValues(r, d3$cum_pct_quantity)

r_cum_pct = setValues(r, d3$cum_pct_quantity)

plot(r2) 

filled.contour()

n_cols = 5
cols = RColorBrewer::brewer.pal(n_cols, 'Spectral')
lvls = seq(0,100, length.out=n_cols+1)
filledContour(r_cum_pct, levels=lvls, col = cols)


levelplot(r_cum_pct, layers = 1, margin = list(FUN = 'median'), contour=TRUE)

x <- rasterToContour(r2, levels=c(10,30,50,80))
x
# rgdal::writeOGR(x, "G:/Team_Folders/Steph/contours", layer="contour_bsb_mean", driver="ESRI Shapefile")


plot(r2, col='Spectral')
plot(x, add=TRUE)

library(leaflet)


binpal <- colorBin("Spectral", seq(0,100), 10, pretty = FALSE, na.color = "transparent")

leaflet() %>% 
  addTiles() %>%
  addProviderTiles('Esri.OceanBasemap') %>%
  addRasterImage(r2, colors = binpal, opacity = 0.6) %>%
  addLegend(
    pal = binpal, values = seq(0,100),
    title = "cum % larvae")


d_20 = d2 %>% filter(cum_pct_quantity >= 20) %>% head(1)

d_40 = d2 %>% filter(cum_pct_quantity >=40) %>% head(1)

d_60 = d2 %>% filter(cum_pct_quantity >= 60) %>% head(1)

max_area_km2 = d2 %>%
  arrange(cum_pct_quantity) %>%
  filter(cum_pct_quantity == 100) %>%
  head(1) %>%
  .$cum_area_km2

formatter1000 <- function(x){ 
  x/1000 
}

plot(r)
p = ggplot(d2, aes(y=cum_pct_quantity, x=cum_area_km2)) +
  ggtitle("All Species 2009 - 2015")+ 
  labs(x="Cumulative Area (km2)",y="Cumulative Larvae (%)") +
  theme(plot.title = element_text(face="bold", size=28)) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title = element_text(face="bold", size=26)) +
  theme(axis.text=element_text(face="bold", size=20, colour = "black"))+
  geom_point()  +
  theme( panel.grid.minor.x = element_blank() ,
    panel.grid.minor.y = element_blank()) +
  geom_segment(x=0, xend=d_20$cum_area_km2, y=d_20$cum_pct_quantity, yend=d_20$cum_pct_quantity, size = 1.5, color = "red") +
  geom_segment(x=d_20$cum_area_km2, xend=d_20$cum_area_km2, y=0, yend=d_20$cum_pct_quantity,size = 1.5, color = "red") +
  geom_segment(x=0, xend=d_40$cum_area_km2, y=d_40$cum_pct_quantity, yend=d_40$cum_pct_quantity,size = 1.5, color = "orange") +
  geom_segment(x=d_40$cum_area_km2, xend=d_40$cum_area_km2, y=0, yend=d_40$cum_pct_quantity,size = 1.5, color = "orange") +
  # geom_segment(x=0, xend=d_60$cum_area_km2, y=d_60$cum_pct_quantity, yend=d_60$cum_pct_quantity,size = 1.5, color = "yellow") +
  # geom_segment(x=d_60$cum_area_km2, xend=d_60$cum_area_km2, y=0, yend=d_60$cum_pct_quantity,size = 1.5, color = "yellow") +
  scale_y_continuous(expand = c(0,0), breaks = c(20,40,60,80,100)) + 
  scale_x_continuous(expand = c(0,0)) +
  #scale_x_continuous() +
  coord_cartesian(xlim = c(0, max_area_km2)) 
  # coord_cartesian(xlim = c(0, tail(d$cum_area_km2, 1)), ylim = c(0, 100))
print(p)

windows()
p


ggplot2::ggsave('all_sp_area_graph.png', p)

ggplotly(p)


plot(r)



# todo ----
# - create github.com/graysreef organization
# - create R package inside github.com/graysreef/mget-conn-process repository
#     using http://ucsb-bren.github.io/env-info/wk07_package.html
# - create Dan's plot: x) cumulative percent larvel input vs y) area of included ranked patches

#Process for percent instead of quanity----
# 
# calculate_percent = function(dir_root ){
#   for (sp in c('bsb','gg','rs','sp')){
#     print(sp)
#     sp_result = sprintf ('%s.*results', sp )
#     print(sp_result)
#     sp_folders = list.files(dir_root, pattern = sp_result, include.dirs = T, full.names = T)
#     print (sp_folders)
#     for (sp_years_folder in sp_folders){
#       print(sp_years_folder)
#       sp_csv = list.files(path = sp_years_folder, pattern = 'csv', recursive = T, full.names = T)
#       print(sp_csv)
# 
#       
#             list.files('G:/Team_Folders/Steph', pattern = 'gg', recursive = T, include.dirs = T, full.names = T)
#       
#       
#     }
#   }
# }
# 
# calculate_percent('G:/Team_Folders/Steph')

library(tidyverse)
library(raster)
library(plotly)

r = raster('G:/Team_Folders/Steph/bsb/mean_sum.tif')

d = data_frame(
  quantity = raster::getValues(r),
  cellid   = 1:length(quantity),
  area_km2 = 64)

d2 = d %>%
  filter(!is.na(quantity)) %>%
  arrange(desc(quantity)) %>%
  mutate(
    pct_quantity     = quantity/sum(quantity)*100,
    cum_pct_quantity = cumsum(quantity/sum(quantity)*100),
    cum_area_km2     = cumsum(area_km2))
tail(d2) # 7208 km2
tail(d2$cum_area_km2, 1) # 7208 km2

d3 = d %>%
  left_join(d2, by='cellid')
summary(d3)

r2 = setValues(r, d3$cum_pct_quantity) #cum_pct_quantity or pct_quantity based on if you want contour or not

plot(r2) 

writeRaster(r2, "percent_sum.tif", format = "GTiff")

#sensitivities ---- 
process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2009_diffusivity')
process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2009_mortality')
process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2009_competency')
process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2009_all')
summarize_sppyr('G:/Team_Folders/Steph/bsb_2009_all')
process_sppyr_dirs('G:/Team_Folders/Steph/rs_2009_competency')
process_sppyr_dirs('G:/Team_Folders/Steph/rs_2009_diffusivity')

# done ----
# process_geodb(
#   'G:/Team_Folders/Steph/bsb_2015/5_4_15_FM_bsb_50day_results',
#   'G:/Team_Folders/Steph/bsb_2015/5_4_15_FM_bsb_50day_simulation')
#process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2015', do_csv=F, do_tif=F, do_png=T)
#summarize_sppyr('G:/Team_Folders/Steph/bsb_2015')
 
# processed speices per Individual year
# process_sppyr_dirs('G:/Team_Folders/Steph/gg_2009')
# process_sppyr_dirs('G:/Team_Folders/Steph/gg_2010')
# process_sppyr_dirs('G:/Team_Folders/Steph/gg_2011')
# process_sppyr_dirs('G:/Team_Folders/Steph/gg_2012')
# process_sppyr_dirs('G:/Team_Folders/Steph/gg_2013')
# process_sppyr_dirs('G:/Team_Folders/Steph/gg_2014')
# process_sppyr_dirs('G:/Team_Folders/Steph/gg_2015')
# 
# summarize_sppyr('G:/Team_Folders/Steph/gg_2009')
# summarize_sppyr('G:/Team_Folders/Steph/gg_2010')
# summarize_sppyr('G:/Team_Folders/Steph/gg_2011')
# summarize_sppyr('G:/Team_Folders/Steph/gg_2012')
# summarize_sppyr('G:/Team_Folders/Steph/gg_2013')
# summarize_sppyr('G:/Team_Folders/Steph/gg_2014')
# summarize_sppyr('G:/Team_Folders/Steph/gg_2015')

# process_sppyr_dirs('G:/Team_Folders/Steph/sp_2009')
# process_sppyr_dirs('G:/Team_Folders/Steph/sp_2010')
# process_sppyr_dirs('G:/Team_Folders/Steph/sp_2011')
# process_sppyr_dirs('G:/Team_Folders/Steph/sp_2012')
# process_sppyr_dirs('G:/Team_Folders/Steph/sp_2013')
# process_sppyr_dirs('G:/Team_Folders/Steph/sp_2014')
# process_sppyr_dirs('G:/Team_Folders/Steph/sp_2015')
# 
# summarize_sppyr('G:/Team_Folders/Steph/sp_2009')
# summarize_sppyr('G:/Team_Folders/Steph/sp_2010')
# summarize_sppyr('G:/Team_Folders/Steph/sp_2011')
# summarize_sppyr('G:/Team_Folders/Steph/sp_2012')
# summarize_sppyr('G:/Team_Folders/Steph/sp_2013')
# summarize_sppyr('G:/Team_Folders/Steph/sp_2014')
# summarize_sppyr('G:/Team_Folders/Steph/sp_2015')

# process_sppyr_dirs('G:/Team_Folders/Steph/rs_2009')
# process_sppyr_dirs('G:/Team_Folders/Steph/rs_2010')
# process_sppyr_dirs('G:/Team_Folders/Steph/rs_2011')
# process_sppyr_dirs('G:/Team_Folders/Steph/rs_2012')
# process_sppyr_dirs('G:/Team_Folders/Steph/rs_2013')
# process_sppyr_dirs('G:/Team_Folders/Steph/rs_2014')
# process_sppyr_dirs('G:/Team_Folders/Steph/rs_2015')
# 
# summarize_sppyr('G:/Team_Folders/Steph/rs_2009')
# summarize_sppyr('G:/Team_Folders/Steph/rs_2010')
# summarize_sppyr('G:/Team_Folders/Steph/rs_2011')
# summarize_sppyr('G:/Team_Folders/Steph/rs_2012')
# summarize_sppyr('G:/Team_Folders/Steph/rs_2013')
# summarize_sppyr('G:/Team_Folders/Steph/rs_2014')
# summarize_sppyr('G:/Team_Folders/Steph/rs_2015')
# 
# process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2009')
# process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2009_all')
# process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2010')
# process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2011')
# process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2012')
# process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2012_all')
# process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2013')
# process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2014')
# process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2015')
# process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2015_all')
# 
# summarize_sppyr('G:/Team_Folders/Steph/bsb_2009')
# summarize_sppyr('G:/Team_Folders/Steph/bsb_2009_all')
# summarize_sppyr('G:/Team_Folders/Steph/bsb_2010')
# summarize_sppyr('G:/Team_Folders/Steph/bsb_2011')
# summarize_sppyr('G:/Team_Folders/Steph/bsb_2012')
# summarize_sppyr('G:/Team_Folders/Steph/bsb_2012_all')
# summarize_sppyr('G:/Team_Folders/Steph/bsb_2013')
# summarize_sppyr('G:/Team_Folders/Steph/bsb_2014')
# summarize_sppyr('G:/Team_Folders/Steph/bsb_2015')
# summarize_sppyr('G:/Team_Folders/Steph/bsb_2015_all')
# 
# summarize_sppyr('G:/Team_Folders/Steph/bsb_2009')
# summarize_sppyr('G:/Team_Folders/Steph/bsb_2010')
# summarize_sppyr('G:/Team_Folders/Steph/bsb_2011')
# summarize_sppyr('G:/Team_Folders/Steph/bsb_2012')
# summarize_sppyr('G:/Team_Folders/Steph/bsb_2013')
# summarize_sppyr('G:/Team_Folders/Steph/bsb_2014')
# summarize_sppyr('G:/Team_Folders/Steph/bsb_2015')
# 
# summarize_sppyr('G:/Team_Folders/Steph/rs_2009')
# summarize_sppyr('G:/Team_Folders/Steph/rs_2010')
# summarize_sppyr('G:/Team_Folders/Steph/rs_2011')
# summarize_sppyr('G:/Team_Folders/Steph/rs_2012')
# summarize_sppyr('G:/Team_Folders/Steph/rs_2013')
# summarize_sppyr('G:/Team_Folders/Steph/rs_2014')
# summarize_sppyr('G:/Team_Folders/Steph/rs_2015')
# 
# summarize_sppyr('G:/Team_Folders/Steph/gg_2009')
# summarize_sppyr('G:/Team_Folders/Steph/gg_2010')
# summarize_sppyr('G:/Team_Folders/Steph/gg_2011')
# summarize_sppyr('G:/Team_Folders/Steph/gg_2012')
# summarize_sppyr('G:/Team_Folders/Steph/gg_2013')
# summarize_sppyr('G:/Team_Folders/Steph/gg_2014')
# summarize_sppyr('G:/Team_Folders/Steph/gg_2015')
# 
# summarize_sppyr('G:/Team_Folders/Steph/sp_2009')
# summarize_sppyr('G:/Team_Folders/Steph/sp_2010')
# summarize_sppyr('G:/Team_Folders/Steph/sp_2011')
# summarize_sppyr('G:/Team_Folders/Steph/sp_2012')
# summarize_sppyr('G:/Team_Folders/Steph/sp_2013')
# summarize_sppyr('G:/Team_Folders/Steph/sp_2014')
# summarize_sppyr('G:/Team_Folders/Steph/sp_2015')
# 
# process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2009_expanded')
# process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2010_expanded')
# process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2011_expanded')
# process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2012_expanded')
# process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2013_expanded')
# process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2014_expanded')
# process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2015_expanded')
# 
# summarize_sppyr('G:/Team_Folders/Steph/bsb_2009_expanded')
# summarize_sppyr('G:/Team_Folders/Steph/bsb_2010_expanded')
# summarize_sppyr('G:/Team_Folders/Steph/bsb_2011_expanded')
# summarize_sppyr('G:/Team_Folders/Steph/bsb_2012_expanded')
# summarize_sppyr('G:/Team_Folders/Steph/bsb_2013_expanded')
# summarize_sppyr('G:/Team_Folders/Steph/bsb_2014_expanded')
# summarize_sppyr('G:/Team_Folders/Steph/bsb_2015_expanded')
# 
# process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2009_from_gr')
# process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2010_from_gr')
# process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2011_from_gr')
# process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2012_from_gr')
# process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2013_from_gr')
# process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2014_from_gr')
# process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2015_from_gr')
# 
# summarize_sppyr('G:/Team_Folders/Steph/bsb_2009_from_gr')
# summarize_sppyr('G:/Team_Folders/Steph/bsb_2010_from_gr')
# summarize_sppyr('G:/Team_Folders/Steph/bsb_2011_from_gr')
# summarize_sppyr('G:/Team_Folders/Steph/bsb_2012_from_gr')
# summarize_sppyr('G:/Team_Folders/Steph/bsb_2013_from_gr')
# summarize_sppyr('G:/Team_Folders/Steph/bsb_2014_from_gr')
# summarize_sppyr('G:/Team_Folders/Steph/bsb_2015_from_gr')

#Processing for mortality because it has a differently named geodatabase----
process_singledir = function(dir_results, dir_simulation, do_csv=T, do_tif=T, do_png=T){
  # dir_results    = 'G:/Team_Folders/Steph/bsb_2015/2_2_15_FM_bsb_50day_results'
  # dir_simulation = 'G:/Team_Folders/Steph/bsb_2015/2_2_15_FM_bsb_50day_simulation'
  
  run = str_replace(basename(dir_results), '_results', '')
  
  # read geodatabase
  conn_lns = readOGR(file.path(dir_results, 'mortality_0.1_A.gdb'), 'Connectivity', verbose=F)
  
  # aggregate across all ToPatchIDs to Gray's Reef (n=4)
  conn_tbl = conn_lns@data %>%
    as_tibble() %>%    
    group_by(FromPatchID) %>%
    summarize(
      quantity = sum(Quantity)) %>%
    ungroup() %>%
    mutate(
      percent = quantity / sum(quantity) * 100) %>%
    arrange(desc(percent))
  
  # write to csv
  if(do_csv){
    write_csv(conn_tbl, sprintf('%s/connectivity.csv', dir_results))
  }
  
  # get patch id raster, and determine which cells are NA
  r_id = raster(sprintf('%s/PatchData/patch_ids', dir_simulation)) # plot(r_id)
  id_NA = !getValues(r_id) %in% conn_tbl$FromPatchID
  
  # create rasters for quantity and percent
  for (v in c('quantity','percent')){
    
    # reclassify from patch id to value
    r = reclassify(r_id, conn_tbl[,c('FromPatchID', v)])
    
    # set patch ids without a value to NA
    r[id_NA] = NA
    
    # write to GeoTIFF
    if(do_tif){
      writeRaster(r, sprintf('%s/%s.tif', dir_results, v), overwrite=T)
    }
    
    
    # plot to PNG for easy preview
    if (do_png){
      png(sprintf('%s/%s.png', dir_results, v))
      p = levelplot(r, par.settings=viridisTheme, main=sprintf('%s %s', run, v))
      print(p)
      dev.off()  
    }
  }
}



