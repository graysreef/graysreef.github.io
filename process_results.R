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

process_singledir = function(dir_results, dir_simulation, do_csv=T, do_tif=T, do_png=T){
  # dir_results    = 'G:/Team_Folders/Steph/bsb_2015/2_2_15_FM_bsb_50day_results'
  # dir_simulation = 'G:/Team_Folders/Steph/bsb_2015/2_2_15_FM_bsb_50day_simulation'
  
  run = str_replace(basename(dir_results), '_results', '')
  
  # read geodatabase
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

process_sppyr_dirs = function(dir_sppyr, ...){
  # process all model runs for given species & year

  dirs_results = list.files(dir_sppyr, '.*_results$', full.names=T)
  for (i in 1:length(dirs_results)){
    
    dir_results = dirs_results[i]
    dir_simulation = str_replace(dir_results, '_results', '_simulation')
    cat(sprintf('%03d of %d: %s\n', i, length(dirs_results), basename(dir_results)))
    
    # process from geodatabase to results csv, tifs, pngs
    process_singledir(dir_results, dir_simulation, ...)
    
  }
}

summarize_sppyr = function(dir_sppyr){

  dirs_results = list.files(dir_sppyr, '.*_results$', full.names=T)
  rasters_quantity = sprintf('%s/quantity.tif', dirs_results)
  stack_quantity = stack(rasters_quantity)
  
  r_mean = mean(stack_quantity, na.rm=T)
  r_sd = calc(stack_quantity, fun=function(x) sd(x, na.rm=T))
  r_cv = r_sd / r_mean * 100
  
  for (v in c('mean','cv')){
    
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
  
  dirs_results = list.files(dir_root, sprintf('%s_[0-9]{4}$', sp), full.names=T)
  rasters_mean = sprintf('%s/mean.tif', dirs_results)
  stack_mean   = stack(rasters_mean)
  dir_sp = file.path(dir_root, sp)
  
  if (!file.exists(dir_sp)) dir.create(dir_sp)
  
  r_mean = mean(stack_mean, na.rm=T)
  r_sd = calc(stack_mean, fun=function(x) sd(x, na.rm=T))
  r_cv = r_sd / r_mean * 100
  
  for (v in c('mean','cv')){
    
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

summarize_spp(dir_root='G:/Team_Folders/Steph', spp=c('bsb','gg','rs','sp'))


####Processing for mortality because it has a differently named geodatabase----
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


##area maps----

library(tidyverse)
library(raster)
library(plotly)

r = raster('G:/Team_Folders/Steph/bsb/mean.tif')

d = data_frame(
  quantity = raster::getValues(r),
  cellid   = 1:length(quantity),
  area_km2 = 8)

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

r2 = setValues(r, d3$cum_pct_quantity)

plot(r2) 

x <- rasterToContour(r2, levels=c(10,30,50,80))
x
rgdal::writeOGR(x, "G:/Team_Folders/Steph/contours", layer="contour_bsb_mean", driver="ESRI Shapefile")


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


d_30 = d2 %>% filter(cum_pct_quantity >= 30) %>% head(1)

plot(r)
p = ggplot(d2, aes(y=cum_pct_quantity, x=cum_area_km2)) +
  geom_point() +
  geom_segment(x=0, xend=d_30$cum_area_km2, y=d_30$cum_pct_quantity, yend=d_30$cum_pct_quantity) +
  geom_segment(x=d_30$cum_area_km2, xend=d_30$cum_area_km2, y=0, yend=d_30$cum_pct_quantity) +
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))
  #coord_cartesian(xlim = c(0, tail(d$cum_area_km2, 1)), ylim = c(0, 100))
print(p)
ggplot2::ggsave('test.png', p)
ggplotly(p)



plot(r)



# todo ----

# for (dir in c('sp_2009','sp_2010','sp_2011','sp_2012', 'sp_2013', 'sp_2014', 'sp_2015')){
#   summarize_sppyr('G:/Team_Folders/Steph/sp_2009')
# }


# - create github.com/graysreef organization
# - create R package inside github.com/graysreef/mget-conn-process repository
#     using http://ucsb-bren.github.io/env-info/wk07_package.html
# - create Dan's plot: x) cumulative percent larvel input vs y) area of included ranked patches

#aggregate csvs ---- 
path <- 'G:/Team_Folders/Steph/rs_2015'
setwd(path)

my.dirs <- dir(pattern = "results", include.dirs = T)


for (i in 1:length(my.dirs)){
  file <- paste0("./",my.dirs[i], "/connectivity.csv")
  print(file)
  my.csv <- read.csv(file)

}

# done ----
# process_geodb(
#   'G:/Team_Folders/Steph/bsb_2015/5_4_15_FM_bsb_50day_results',
#   'G:/Team_Folders/Steph/bsb_2015/5_4_15_FM_bsb_50day_simulation')
#process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2015', do_csv=F, do_tif=F, do_png=T)
#summarize_sppyr('G:/Team_Folders/Steph/bsb_2015')


##sensitivities
process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2009_diffusivity')
process_sppyr_dirs('G:/Team_Folders/Steph/bsb_2009_mortality')


# processed speices per Individual year---- 
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

#Process for percent instead of quanity----

summarize_spp = function(dir_root, sp){
  # given top-level directory and species code, eg "sp" or "rs" or "bsb",
  # summarize sp_yr/mean.tif across years as sp/mean.tif and sp/cv.tif,
  # ie average dispersal across year means and variation across year means
  # dir_root = 'G:/Team_Folders/Steph'; sp='bsb'
  
  dirs_results = list.files(dir_root, sprintf('%s_[0-9]{4}$', sp), full.names=T)
  rasters_mean = sprintf('%s/mean.tif', dirs_results)
  stack_mean   = stack(rasters_mean)
  dir_sp = file.path(dir_root, sp)
  
  if (!file.exists(dir_sp)) dir.create(dir_sp)
  
  r_mean = mean(stack_mean, na.rm=T)
  r_sd = calc(stack_mean, fun=function(x) sd(x, na.rm=T))
  r_cv = r_sd / r_mean * 100
  r_percent =  r_mean/sum(r_mean) * 100
  
  for (v in c('mean','cv', 'percent')){
    
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

summarize_spp(dir_root='G:/Team_Folders/Steph', spp=c('bsb','gg','rs','sp'))


####Processing for mortality because it has a differently named geodatabase----
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



