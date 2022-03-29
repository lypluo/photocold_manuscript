##########################################
#To check the MODIS VIs(CCI and PRI et al.,) for Fluxnet sites
##########################################
library(lubridate)
library(dplyr)
library(tidyr)
#------------------------------------
#(1)load the daily flux-sites VI data
#------------------------------------
data.path<-"D:/data/FLUXNET_MODIS_VIs_extended/"
load(paste0(data.path,"MOD09GA_MODOCGA_filter_indices.RDA"))

#
df.VIs<-filter_state_data[,c("YY","MM","DD","sites_id","ndvi","evi","NIRv","cci","pri")]
df.VIs$date<-ymd(paste0(df.VIs$YY,"-",df.VIs$MM,"-",df.VIs$DD))
#re-tidy the names to keep the consistence with other datasets:
df.VIs_new<-df.VIs[,c("sites_id","date","ndvi","evi","NIRv","cci","pri")]
names(df.VIs_new)<-c("sitename","date","ndvi","evi","NIRv","cci","pri")
df.VIs_new$NIRv<-df.VIs_new$NIRv*0.0001##needs to multiple scaling factor 0.0001
#------------------------------------
#(2)sites examples
#------------------------------------
library(ggplot2)
norm_ts<-function(x){
  mm<-max(x,na.rm = T)
  mn<-min(x,na.rm = T)
  norm_x<-(x-mn)/c(mm-mn)
  return(norm_x)
}
df_sel_site<-df.VIs_new %>%
  filter(sitename=="FI-Hyy")
df_sel_site %>%
  ggplot(aes(x=date,y=cci))+
  geom_point()

# fast processing
df_sel_site$pri[df_sel_site$pri<c(-0.3)]<-NA
df_sel_site$cci[df_sel_site$cci<c(-0.6)]<-NA
df_sel_site<-cbind(df_sel_site[,c("sitename","date")],
                   as.data.frame(apply(df_sel_site[,-c(1:2)],2,norm_ts)))
df_sel_site$doy<-yday(df_sel_site$date)


df_sel_site %>%
  group_by(doy)%>%
  summarise(ndvi=mean(ndvi,na.rm=T),
            evi=mean(evi,na.rm=T),
            NIRv=mean(NIRv,na.rm=T),
            cci=mean(cci,na.rm=T),
            pri=mean(pri,na.rm=T)) %>%
  pivot_longer(c(ndvi,evi,NIRv,cci),names_to = "Sources",values_to = "VI") %>%
  ggplot(aes(x=doy,y=VI,col=Sources))+
  geom_point()+
  geom_line()+
  geom_vline(xintercept = c(110,130),lty=2)
  # facet_wrap(~Sources,scales = "free_y")

#------------------------------------
#(3)needs to filter the VIs data for each site (to remove the extreme values):
#------------------------------------
#I. first have an overview the VIs variation range in differnt sites
apply(df.VIs_new[,-c(1:2)],2,function(x){range(x,na.rm = T)})
#
hist(df.VIs_new$ndvi)
hist(df.VIs_new$evi)
hist(df.VIs_new$NIRv)
hist(df.VIs_new$cci)
hist(df.VIs_new$pri) #mainly focus on the VIs beyond PRI
#extreme values
apply(df.VIs_new[,-c(1:2)],2,
function(x){quantile(x,probs=c(0,0.005,0.01,0.05,seq(0.1,0.9,0.1),0.95,0.99,0.995,1),na.rm=T)})

#II.from the quantile-->we alreaday see VIs between[0.5%,99.5%] are looks normal
#hence remove the outlier by remove the data beyond quantile 0.5% and 99.5%
var_names<-names(df.VIs_new)
df_temp<-c()
for (i in 3:length(var_names)) {
  temp<-df.VIs_new[,i]
  q_values<-quantile(temp,probs = c(0.005,0.995),na.rm = T)
  #assign the values beyond 0.5 and 99.5 percentile are NAs
  temp[temp<q_values[1]|temp>q_values[2]]<-NA
  df_temp<-cbind(df_temp,temp)
}
#
df.temp<-as.data.frame(df_temp)
df.VIs_final<-cbind(df.VIs_new[,c(1:2)],df.temp)
names(df.VIs_final)<-var_names
#
apply(df.VIs_final[,-c(1:2)],2,function(x){range(x,na.rm = T)})

##save the data:
save.path<-"./data-raw/raw_data/processed_additional_VIs_from_MODIS/"
save(df.VIs_final,file=paste0(save.path,"df_modis_CCI_etal.RDA"))

#########################additional code(: add comments by YP in March, 2022#################
##III. make the plots for each sites(for sites in analyzed boreal and tempearate forest):
#load the event information
load.path<-"D:/data/photocold_project/event_length/Using_sites_in_Fluxnet2015/"
load(paste0(load.path,"df_events_length.RDA"))

df.VIs_final$doy<-yday(df.VIs_final$date)

#add "event" informaiton:
df.VIs_final$Year<-year(df.VIs_final$date)
df.VIs_end<-left_join(df.VIs_final,
                      df_events_all[,c("sitename","Year","sos","peak","Over_days_length")],
                      by=c("sitename","Year"))

##load the PFT information:
library(tidyverse) #-->function read_rds
PFT.path<-"D:/Github/photocold/data/"
df_sites_modis_era <- read_rds(paste0(PFT.path,"df_sites_modis_era.csv"))
#
df.VIs_end<-df.VIs_end %>%
  left_join(
    df_sites_modis_era[,c("sitename","classid","koeppen_code")],
    by = "sitename"
  )

#only keep the sits that analyzed in this study:
df.VIs_end<-df.VIs_end[!is.na(df.VIs_end$sos),]

#normalize the VIs in each site
df_norm<-df.VIs_end %>%
  group_by(sitename) %>%
  select(classid,date,doy,Year,sos,peak,Over_days_length,ndvi,evi,NIRv,cci) %>%
  mutate(ndvi_norm=norm_ts(ndvi),
         evi_norm=norm_ts(evi),
         NIRv_norm=norm_ts(NIRv),
         cci_norm=norm_ts(cci))

#plotting:
#-------------------
#only for DBF
#------------------
PFT_name<-"DBF"
thres_overdays<-50
#(1) For over_day_length>20
###############
# df.events_eachsite<-df_norm %>% 
#   filter(Over_days_length >30) %>%
#   group_by(sitename,doy) %>%
#   summarise(NDVI=mean(ndvi_norm,na.rm=T),
#             EVI=mean(evi_norm,na.rm=T),
#             NIRv=mean(NIRv_norm,na.rm=T),
#             CCI=mean(cci_norm,na.rm=T)) %>%
#   pivot_longer(c(NDVI,EVI,NIRv,CCI),names_to = "Sources",values_to = "VIs")
# 
# df.events_eachsite %>%
#   # filter(sitename=="FI-Hyy") %>%
#   ggplot(aes(x=doy,y=VIs,col=Sources))+
#   geom_point()+
#   geom_line()+
#   facet_wrap(~sitename)
##################
df.events_allsite<-df_norm %>% 
  filter(classid==PFT_name)%>%
  filter(Over_days_length >thres_overdays) %>%
  group_by(doy) %>%
  summarise(NDVI=mean(ndvi_norm,na.rm=T),
            EVI=mean(evi_norm,na.rm=T),
            NIRv=mean(NIRv_norm,na.rm=T),
            CCI=mean(cci_norm,na.rm=T)) %>%
  pivot_longer(c(NDVI,EVI,NIRv,CCI),names_to = "Sources",values_to = "VIs")

df_plot1<-df.events_allsite %>%
  group_by(doy) %>%
  ggplot(aes(x=doy,y=VIs,col=Sources))+
  geom_point()+
  geom_line()+
  annotate(geom="text",x=50,y=0.95,label=paste0(PFT_name,"-event"))+
  annotate(geom = "text",x=180,y=0.1,label=paste0("events-thresholds:",thres_overdays))
  
  

#(2) For over_day_length<20
###############
# df.nonevents_eachsite<-df_norm %>% 
#   filter(Over_days_length <30) %>%
#   group_by(sitename,doy) %>%
#   summarise(NDVI=mean(ndvi_norm,na.rm=T),
#             EVI=mean(evi_norm,na.rm=T),
#             NIRv=mean(NIRv_norm,na.rm=T),
#             CCI=mean(cci_norm,na.rm=T)) %>%
#   pivot_longer(c(NDVI,EVI,NIRv,CCI),names_to = "Sources",values_to = "VIs")
# 
# df.nonevents_eachsite %>%
#   # filter(sitename=="FI-Hyy") %>%
#   ggplot(aes(x=doy,y=VIs,col=Sources))+
#   geom_point()+
#   geom_line()+
#   facet_wrap(~sitename)

##################
df.nonevents_allsite<-df_norm %>% 
  filter(classid==PFT_name)%>%
  filter(Over_days_length <thres_overdays) %>%
  group_by(doy) %>%
  summarise(NDVI=mean(ndvi_norm,na.rm=T),
            EVI=mean(evi_norm,na.rm=T),
            NIRv=mean(NIRv_norm,na.rm=T),
            CCI=mean(cci_norm,na.rm=T)) %>%
  pivot_longer(c(NDVI,EVI,NIRv,CCI),names_to = "Sources",values_to = "VIs")

df_plot2<-df.nonevents_allsite %>%
  group_by(doy) %>%
  ggplot(aes(x=doy,y=VIs,col=Sources))+
  geom_point()+
  geom_line()+
  annotate(geom = "text",x=50,y=0.95,label=paste0(PFT_name,"-nonevent"))+
  annotate(geom = "text",x=180,y=0.1,label=paste0("events-thresholds:",thres_overdays))
#
library(cowplot)
df_plot1<-df_plot1+
  xlim(0,365)+
  ylim(0,1)
df_plot2<-df_plot2+
  xlim(0,365)+
  ylim(0,1)

plot_grid(df_plot1,df_plot2,
          labels = "auto",nrow=1,label_size = 18,align = "hv")

