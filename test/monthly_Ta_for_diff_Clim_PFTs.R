#######################################################
##Aim: comaparing the variable differences between different Clim-PFTs
#here especially aimed for comparing the monthly Ta
#######################################################
#----------
library(dplyr)
library(tidyverse)
library(GenSA)
library(lubridate)
#-------------------------
#(1)load the data and hardening funciton
#-------------------------
#####
#load the data updated p-model data uploaded by Koen
df_recent <- readRDS(paste0("./data-raw/raw_data/P_model_output/model_data.rds")) %>%
  mutate(
    year = format(date, "%Y")
  ) %>%
  na.omit()

sites<-unique(df_recent$sitename)
for (i in 1:length(sites)) {
  df_temp<-df_recent %>%
    filter(sitename==sites[i])
  text.Date<-min(df_temp$date)+c(max(df_temp$date)-min(df_temp$date))*0.1
  #
  df_plot<-df_temp %>%
    ggplot()+
    geom_point(aes(x=date,y=gpp))+
    geom_point(aes(x=date,y=gpp_mod),col="red")+
    annotate(geom = "text",x=text.Date,y=15,label=sites[i])
  #
  k=quantile(df_temp$gpp,probs = c(0.01,0.05,seq(0.1,0.9,0.1)),0.95,0.99)
  df_plot+
    geom_hline(yintercept = k,col="blue")
}

# #filter the observational gpp data:
# df_recent_new<-c()
# for (i in 1:length(sites)) {
#   df_temp<-df_recent %>%
#     filter(sitename==sites[i])
#   #filter the gpp observation data(remove the gpp that below 5 percentile and negative):
#   k=quantile(df_temp$gpp,probs = c(0.01,0.05,seq(0.1,0.9,0.1)),0.95,0.99)
#   df_temp$gpp[df_temp$gpp<as.numeric(k[2])& df_temp$gpp<0]<-NA
#   # df_temp %>%
#   #   ggplot()+
#   #   geom_point(aes(x=date,y=gpp))+
#   #   geom_point(aes(x=date,y=gpp_mod),col="red")+
#   #   annotate(geom = "text",x=text.Date,y=15,label=sites[i])
#   df_recent_new<-rbind(df_recent_new,df_temp)
# }

#do not filter the observation gpp in study:
df_recent_new<-df_recent

#load the data Beni sent me before:
df_old<-read.csv(file=paste0("./data-raw/raw_data/Data_sent_by_Beni/","ddf_fluxnet2015_pmodel_with_forcings_stocker19gmd.csv"))
df_old<-df_old %>%
  mutate(date=lubridate::mdy(date),
         year=lubridate::year(date)) %>%
  na.omit(gpp_obs)

#--------------------------------------------------------------
#(3)merge datasets
#--------------------------------------------------------------
#first load the PFTs information:
#load the modis data-->tidy from Beni
load(paste0("./data-raw/raw_data/sites_info/","Pre_selected_sites_info.RDA"))
sites.info<-df_sites_sel
#
df_merge<-df_recent_new %>%
left_join(
    sites.info,
    by = "sitename"
  )
#day of the year
library(sirad)  #calculate the day of the year
df_merge$doy<-dayOfYear(df_merge$date)

#main Clim-PFTs
df_merge$Clim_PFTs<-paste0(df_merge$koeppen_code,"-",df_merge$classid)
Clim.PFTs<-sort(unique(df_merge$Clim_PFTs))

#--------------------------------------------------------------
#(4) compare environmental variables between specific Clim-PFTs types:
#here for monthly Ta
#--------------------------------------------------------------
df_final_new<-df_merge
# need to remove the sites that do not used in this analysis:
rm.sites<-c("BE-Bra","CA-SF1","CA-SF2","FI-Sod","US-Wi4")
df_final_new<-df_final_new %>%
  filter(sitename!=rm.sites[1] & sitename!=rm.sites[2]&sitename!=rm.sites[3]&sitename!=rm.sites[4]&sitename!=rm.sites[5])
##------------
#calculate the monthly temperature:
##------------
df_final_new %>%
  mutate(Month=month(date))%>%
  select(sitename,date,temp,tmin,tmax,Month,Clim_PFTs)%>%
  group_by(Clim_PFTs,Month)%>%
  summarise(Tmean=mean(temp,na.rm=T)) %>%
  ggplot(aes(x=Month,y=Tmean))+
  geom_point()+
  geom_line()+
  geom_vline(xintercept = c(3,4,5),col="blue",lty=2)+
  geom_hline(yintercept = c(0,5),col="red")+
  facet_wrap(.~Clim_PFTs,nrow = 3)

##from the results-->all the Clim-PFTs will above 0 degree after April-->
##paramaterize the model from April...
