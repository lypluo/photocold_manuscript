##########################################
#Try to merge the data from Fluxnet and PhenoCam
##########################################
#------------------------------------
#(1)load the daily flux-sites data
#------------------------------------
library(lubridate)
library(dplyr)
#---------------
##a. load the processed fluxnet sites data done by me:
load.path<-"./data-raw/raw_data/processed_data_from_FLUX2015/"
load(paste0(load.path,"Daily_data.RDA"))
#sites form Fluxnet2015-->correspounding to Beni's datasets
df_YP_daily<-df_all_sel_daily
df_YP_daily$Date<-ymd(df_YP_daily$Date)
#rename the variables:
#"ppfd_fluxnet2015" corresponds to "ppfd_in_fluxnet2015"--> to combine with Beni's datasets, I used the first name
#adopt the "temp_day_fluxnet" name from Beni's dataset to represent mean daily temperature
names(df_YP_daily)<-c("sitename","date","temp_day_fluxnet2015","temp_min_fluxnet2015","temp_max_fluxnet2015",
                      "SW_IN_fluxnet2015","patm_fluxnet2015","prec_fluxnet2015","ws_fluxnet2015",
                      "ppfd_fluxnet2015","ppfd_out_fluxnet2015","vpd_day_fluxnet2015",
                      "nee","gpp_nt","gpp_dt",
                      "TS_1_fluxnet2015","TS_2_fluxnet2015","TS_3_fluxnet2015","TS_4_fluxnet2015","TS_5_fluxnet2015","TS_6_fluxnet2015",
                      "TS_7_fluxnet2015",
                      "SWC_1_fluxnet2015","SWC_2_fluxnet2015","SWC_3_fluxnet2015","SWC_4_fluxnet2015","SWC_5_fluxnet2015",
                      "SWC_6_fluxnet2015","SWC_7_fluxnet2015")

#---------------------------
##b. load Beni processed data:
load.path<-"./data-raw/raw_data/Data_sent_by_Beni/"
df_Beni<-read.csv(file=paste0(load.path,"ddf_fluxnet2015_pmodel_with_forcings_stocker19gmd.csv"))
##load the selected sites (from Fluxnet2015) info for anlysis 
load.path<-"./data-raw/raw_data/sites_info/"
load(paste0(load.path,"Pre_selected_sites_info.RDA"))
##check how many sites data available:
pos<-match(df_sites_sel$sitename,unique(df_Beni$sitename))
sel_pos<-pos[!is.na(pos)]
avai.sites<-unique(df_Beni$sitename)[sel_pos]
df_sites_avai<-df_sites_sel %>%
  filter(sitename %in% avai.sites)
#available sites data for analysis:
df_avai<-df_Beni %>%
  filter(sitename %in% avai.sites)
df_avai$date<-mdy(df_avai$date)
#save the available sites info:
save(df_sites_avai,file=paste0(load.path,"Available_sites_info.RDA"))

#-----------------------------------------
#c. keep the variables that Beni's datasets(df_avai) have, and add other variables from df_all_daily
#to check how much data are available in different datasets:
# visdat::vis_miss(df_avai, warn_large_data = FALSE)
# visdat::vis_miss(df_YP_daily,warn_large_data = FALSE)

pos<-match(names(df_avai[,-c(1:2)]),names(df_YP_daily[,-c(1:2)]))
pos_del<-pos[!is.na(pos)]+2  #real position should +2
df_YP_daily_temp<-df_YP_daily[,-pos_del]

#get the merged datasets-->merge both df_avai and df_YP_daily_temp
library(reshape2)
df_merge_daily<-merge(df_avai,df_YP_daily_temp,by = c("sitename","date"),all.x = T)

#--------------------------------------------
#d. load the updated drivers data from Koen:
load.path<-"./data-raw/raw_data/Data_prep_by_Koen/"
load(paste0(load.path,"p_model_fluxnet_drivers.rda"))
df_Koen<-p_model_fluxnet_drivers
#compare the differences:
setdiff(unique(df_merge_daily$sitename),unique(df_Koen$sitename))
sites_sel<-unique(df_merge_daily$sitename)
df_Koen_sel<-df_Koen %>%
  filter(sitename %in% sites_sel)
#convert to data.frame
df_Koen.forcing<-c()
for (i in 1:nrow(df_Koen_sel)) {
  forcing.temp<-as.data.frame(df_Koen_sel[i,]$forcing)
  df_temp<-data.frame(sitename=rep(df_Koen_sel[i,]$sitename,nrow(forcing.temp)),forcing.temp)
  df_Koen.forcing<-rbind(df_Koen.forcing,df_temp)
}
#
visdat::vis_miss(df_Koen.forcing,warn_large_data = FALSE)
visdat::vis_miss(df_merge_daily,warn_large_data = FALSE)
#comparison between df_Koen and df_merge_daily for tmin and tmax:
tt1<-df_merge_daily[,c("sitename","date","temp_day_fluxnet2015",
                       "temp_min_fluxnet2015","temp_max_fluxnet2015","vpd_day_fluxnet2015")]
tt2<-df_Koen.forcing[,c("sitename","date","temp","tmin","tmax","vpd")]
tt<-left_join(tt1,tt2,by=c("sitename","date"))
#well matched!
plot(tt$temp,tt$temp_day_fluxnet2015)
plot(tt$tmin,tt$temp_min_fluxnet2015)
abline(0,1,col="blue",lty=2)
plot(tt$tmax,tt$temp_max_fluxnet2015)
abline(0,1,col="blue",lty=2)
plot(tt$vpd,tt$vpd_day_fluxnet2015)
abline(0,1,col="blue",lty=2)
#using data from Koen to substitute the tmin and tmax in df_merge_daily:
df_Koen.forcing_sel<-df_Koen.forcing[,c("sitename","date","tmin","tmax")]
df_merge_daily<-left_join(df_merge_daily,df_Koen.forcing,by=c("sitename","date"))
df_merged_daily<-df_merge_daily %>%
  mutate(temp_min_fluxnet2015=tmin,
         temp_max_fluxnet2015=tmax,
         tmin=NULL,
         tmax=NULL)
#------------------------------------
#(2)load the daily PhenoCam data and merge to df_flux_merge
#------------------------------------
load.path<-"./data-raw/raw_data/processed_data_from_PhenoCam/"
load(paste0(load.path,"Daily_data.RDA"))
df.Phenocam_daily$date<-as.Date(df.Phenocam_daily$date)
#
df_merge<-merge(df_merge_daily,df.Phenocam_daily,by=c("sitename","date"),all.x = T)

#------------------------------------
#(3)load the daily VIs including CCI, EVI...
#------------------------------------
load.path<-"./data-raw/raw_data/processed_additional_VIs_from_MODIS/"
load(paste0(load.path,"df_modis_CCI_etal.RDA"))
#
df_merge<-left_join(df_merge,df.VIs_final,by=c("sitename","date"))

#save the merged data:
save.path<-"./data-raw/raw_data/Merged_data/"
save(df_merge,file=paste0(save.path,"Merged_Flux_and_VIs.RDA"))

