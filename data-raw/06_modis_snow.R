#-----------------------------------------
#Aim: get the snow data from MODIS
#-----------------------------------------
#load the site info:
site.info<-"./data-raw/raw_data/sites_info/"
load(paste0(site.info,"Available_sites_info.RDA"))

#
write.csv(df_sites_avai,file=paste0(site.info,"sites_info.csv"))
