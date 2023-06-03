#-----------------------------------------
#Aim: tidy the snow data from MODIS 
#-----------------------------------------
library(dplyr)
library(tidyverse)
#---------------------
#1) load the site info:
#---------------------
site.info<-"./data-raw/raw_data/sites_info/"
load(paste0(site.info,"Available_sites_info.RDA"))
#
#write.csv(df_sites_avai,file=paste0(site.info,"sites_info.csv"))
#then download the snow product from MODIS through the following link:
# https://appeears.earthdatacloud.nasa.gov/task/point

#------------------------
#2) after download the data, tidy the snow data
#------------------------
#data are download from MODIS sensor onboard on Terra satellite 
load.path<-"D:/data/MODIS_Snow_product/MODIS_Snow_MOD/download_site_for_photocold_project/"
snow_MOD<-read.csv(paste0(load.path,"Fluxnet-sites-MOD10A1-061-results.csv"))
snow_MOD_sel<-snow_MOD %>%
  select(ID, Latitude,Longitude,Date, 
         MOD10A1_061_NDSI_Snow_Cover,MOD10A1_061_NDSI_Snow_Cover_Basic_QA)
#rename the variables
#also refer the user guide from "D:/data/MODIS_Snow_product/MODIS_Snow_MOD/download_site_for_photocold_project/"
names(snow_MOD_sel)<-c("sitename","latitude","longtitude","date",
                       "snowval","snowval_QA")

#-----------------------------
#3)to further tidy the dataset
#-----------------------------
#the meaning of the snowval:
# 0-100=fractional 
# snow, 
# 200=missing data, 
# 201=no decision, 
# 211=night, 
# 225=land, 
# 237=inland water, 
# 239=ocean, 
# 250=cloud, 
# 254=detector 
# saturated, 255=fill
# to be simple-->just keep the values when values <=100-->other set to NA
snow_MOD_sel<-snow_MOD_sel%>%
  mutate(snowval=ifelse(snowval<=100,snowval,NA))
#---------------------
#3)extrapolate the snow data
#---------------------
library(lubridate)
library(zoo)
snow_MOD_sel<-snow_MOD_sel %>%
  mutate(date=as.Date(date),
         year=year(date),
         doy=yday(date))%>%
  group_by(sitename,year)%>%
  mutate(snowval_fill=na.fill(snowval,"extend"))

##example plot:
snow_MOD_sel %>%
  filter(sitename=="BE-Bra")%>%
  ggplot()+
  geom_point(aes(x=date,y=snowval_fill),col="red")+
  geom_point(aes(x=date,y=snowval))
  
#---------------------
#4)save the data
#---------------------
save.path<-"./data-raw/raw_data/snow_data/"
save(snow_MOD_sel,file = paste0(save.path,"snow_MOD.RDA"))
