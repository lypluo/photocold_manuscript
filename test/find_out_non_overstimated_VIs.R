#--------------------------------------------
#Aim: to find out the sites(especially GPP non-overestimated sites) that have PhenoCam VIs
#--------------------------------------------
library(lubridate)
#---------------------
#(I).Calculate the mean T in the first half-year:(Jan-Jun) 
#for GPP overestimated and non-overestimatd sites
#---------------------
#---------------------
#A.load the event_length data
#---------------------
load.path<-"./data/event_length/"
load(paste0(load.path,"df_events_length.RDA"))
# a function to separate out the site-year that event_length higher than some thresholds(e.g. 30 days)
sep_siteyears<-function(df,sep_var,sep_thrs){
  # df<-df_events_all
  # sep_var<-"Over_days_length"
  # sep_thrs<-30
  #
  df_sep<-df[df[,sep_var]>sep_thrs & !is.na(df[,sep_var]),]
  pos_sep<-as.vector(which(df[,sep_var]>sep_thrs))
  #as some site the sep_var==NA, hence using this way to separate the data
  pos_left<-setdiff(as.vector(1:length(df[,sep_var])),pos_sep)
  df_left<-df[pos_left,]
  df_new<-list(df_sep,df_left)
  names(df_new)<-c("event_siteyears","noevent_siteyears")
  return(df_new)
}
##separate the site-year when the over_days_length > 20 days
df.sep20<-sep_siteyears(df_events_all,"Over_days_length",20)

#---------------------
#B.load the data 
#---------------------
load.path<-"./data/data_used/"
#from new method:
load(paste0(load.path,"ddf_labeled_norm_trs_newmethod_all_overestimation_Fluxnet2015_sites.RDA"))
df_all_sites<-ddf_labeled;rm(ddf_labeled) 
df_norm_all<-df_all_sites

#-------------------------------------------------------------------------
#C.start to align the data according to Beni's functions of "align_events" and "get_consecutive"
#-------------------------------------------------------------------------
#---------------------------------
#source the functions to detect the consecutive events
#---------------------------------
fun.path<-"./R/functions_from_beni/"
#source get_consecutive.R
source(paste0(fun.path,"get_consecutive.R"))
#source get_consecutive_greenup.R
source(paste0(fun.path,"get_consecutive_greenup.R"))

#source align_event.R -->for event site years
source(paste0(fun.path,"align_events_df.R"))
#source align_nonevent.R -->for non-event site years
source(paste0(fun.path,"align_nonevents_df.R"))
#-------------------------------------------------
#add additional code for PFTs: add PFTs information for each sites-->2021-12-13
#------------------------------------------------
#load the modis data for classifying the vegetation types:
library(tidyverse) #-->function read_rds
PFT.path<-"./data-raw/raw_data/sites_info/"
load(paste0(PFT.path,"Available_sites_info.RDA"))
#
df_norm_all<-df_norm_all%>%
  left_join(
    df_sites_avai,
    by = "sitename"
  )

sep_siteyears_data<-function(df.data,dovars,df.sep,leng_threshold,before,after,nbins,do_norm){
  # df.data<-df_norm_all
  # dovars<-c("gpp_obs")
  # df.sep<-df.sep20
  # leng_threshold<-5
  # before=30
  # after=0
  # nbins=10
  # do_norm=TRUE
  
  #-------------------------------------------------------------------------
  #A.separating the datasets:"event" and "nonevent" site years
  #-------------------------------------------------------------------------
  N<-nrow(df.sep$event_siteyears)
  df_events.years<-c()    #target
  df_nonevents.years<-c() #target
  #for is.event years
  pos.isevent<-c()
  for(i in 1:N){
    df.temp<-subset(df.data,sitename==df.sep$event_siteyears$sitename[i] & Year==df.sep$event_siteyears$Year[i])
    pos.temp<-c(which(df.data$sitename==df.temp$sitename & df.data$Year==df.temp$Year))
    
    df_events.years<-rbind(df_events.years,df.temp)
    pos.isevent<-c(pos.isevent,pos.temp)
  }
  #for no is.event years
  pos.non_isevent<-setdiff(c(1:nrow(df.data)),pos.isevent)
  df_nonevents.years<-df.data[pos.non_isevent,] #target
  
  #---------------------------------
  #B.start to align different events(for event site years) and green-up(for non-event site years) in each site-year
  #---------------------------------
  # df<-df_norm_all
  # leng_threshold<-20
  # before=30
  # after=30
  # nbins=10
  # do_norm=FALSE
  
  #function format
  # align_events <- function( df, leng_threshold, before, after, nbins, do_norm=FALSE )
  #at this stage-->bins do not used(now set default value as 10)
  #-------------
  #set different length_threshold-->5 days(consecutive 5 days overestimation will be an event)
  #-------------
  #and select the 30 days before the events 
  df_len_events<-align_events(df_events.years,dovars,leng_threshold = leng_threshold,before = before,after = after,
                              nbins = nbins,do_norm = do_norm)
  print("ok")
  df_len_nonevents<-align_nonevents(df_nonevents.years,dovars,leng_threshold = leng_threshold,before = before,after = after,
                                    nbins = nbins,do_norm = do_norm)
  #---------------------------------
  #C.separate the df_event and df_nonevent;
  #---------------------------------
  df_all<-c()
  #for event site years
  df_dday<-df_len_events$df_dday
  #for non_event site years, take the doy belongs to green-up period:
  df_noevent_dday<-df_len_nonevents$df_dday
  #
  df_all<-list(df_dday=df_dday,df_noevent_dday=df_noevent_dday)
  return(df_all)
}
#!!important step:leng_threshold=5-->merge the events(consecutive days are over 5 days)
#do_vars-->the variables that are going to be processed(normalized)-->actually do not processed in this process 
# names(df_norm_all)
# do_vars<-c("gpp_obs","fapar_itpl","fapar_spl",paste0(c("ppfd","temp_day","temp_min","temp_max",
#                                                        "vpd_day","prec","patm","SW_IN","ws",paste0("TS_",1:7),paste0("SWC_",1:5)),"_fluxnet2015"),
#            "gcc_90","rcc_90")
#set the before events days from 30 days to 60 days
df_len5_nonnorm<-sep_siteyears_data(df_norm_all,do_vars,df.sep20,5,60,0,10,FALSE)
#update in March, 2022:additional-->calculate the doy range when the overestimate happens in different PFTs:
#for sites have overesimation:
df_overestimated_range<-df_len5_nonnorm$df_dday %>%
  filter(is_event=="yes")%>%
  group_by(sitename,classid)%>%
  dplyr::summarise(doy_min=range(doy)[1],doy_max=range(doy)[2])%>%
  group_by(classid)%>%
  dplyr::summarise(doy_min=floor(mean(doy_min)),doy_max=ceiling(mean(doy_max)))
df_overestimated_range
#---------------------------------
#D.Calculating the Ta in first half-year(Jan-June) in overestimated and non-overestimated site-years
#for different PFTs
#---------------------------------
Cal_1half_Year_Ta<-function(df){
  # df<-df_len5_nonnorm$df_dday
  
  #
  df$month<-month(df$date)
  #
  each_site_Tmean<-df %>%
    group_by(classid,sitename) %>%
    filter(month<=6) %>%
    dplyr::summarise(Ta=mean(temp,na.rm=T))
  multi_sites_Tmean<-each_site_Tmean %>%
    group_by(classid)%>%
    dplyr::summarise(Tmean=mean(Ta,na.rm = T),
                     Tsd=sd(Ta,na.rm = T))
  Tmean_all<-list(each_site_Tmean=each_site_Tmean,
                  multi_sites_Tmean=multi_sites_Tmean)
  return(Tmean_all)
}
#
Tmean_over<-Cal_1half_Year_Ta(df_len5_nonnorm$df_dday)
Tmean_nonover<-Cal_1half_Year_Ta(df_len5_nonnorm$df_noevent_dday)

#---------------------
#(II).Calculate the mean T in the first half-year:(Jan-Jun) 
#for PhenoCam sites
#---------------------
#A.load the PhenoCam info:
load(paste0("./data-raw/raw_data/processed_data_from_PhenoCam/Daily_data.RDA"))
PhenoCam.sites<-unique(df.Phenocam_daily$sitename) #23 sites

#B.load the Meteological data tidy by Jiangong:
load.path<-"D:/Github/FallPheno/data/"
load(paste0(load.path,"flux_meteo_VIs_fromJiangong.rda"))
df.meteo<-df;rm(df)

#C.load the site info tidy by Jiangong:
site.info<-read.csv(paste0(load.path,"SiteInfo.csv"))
site.info<-site.info[,c(1,4:6)]
names(site.info)<-c("sitename","PFT","Lat.","Long.")
#
Tmean_Camsites<-df.meteo %>%
  filter(sitename %in% PhenoCam.sites)%>%
  filter(Month<=6) %>%
  group_by(sitename)%>%
  dplyr::summarise(Tmean=mean(TA_F))
names(site.info)
Tmean_Camsites<-left_join(Tmean_Camsites,site.info) #18 sites available

##################################
#according to the above analysis and compare between the mean T between "non-overestimation"
#and "overestimation" sites in different PFTs-->summarize as below:
#for DBF:put US-UMB, US-UMd, US-WCr, CA-Gro as sites under low-Ta stress; while
#        put DE-Hai, DK-Sor, and DE-Lnf as sites with moderate or no Ta stress

#for ENF:put CA-Qfo,US-NR1,FI-Hyy,CA-Obs as sites under low-Ta strees, while
#        put DE-Tha,CA-TP1,CA-TP3, CA-TP4 as sites with moderate and no Ta stress

#the GPP overestimation period are as follows for differnt PFTs:
# DBF: 70-103
# ENF: 81-143
# MF: 69-147
# for all PFTs range: 69-147
apply(df_overestimated_range[,2:3],2,mean)
