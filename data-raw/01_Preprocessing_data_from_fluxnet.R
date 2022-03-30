##########################################
#Try to preprocess the fluxnet data from FLUXNET2015
##########################################
library(lubridate)
library(zoo)
library(dplyr)
library(plyr)
#------------------------------------
#(1)for fluxnet data:
#------------------------------------
#----------------------
#a. unzip the zip file:for all the available sites from FLUXNET2015
#----------------------
unzip_file<-function(ori.zip.path,sel_site,exdir.path){
  # ori.zip.path<-"C:/Users/yluo/Desktop/CES/Data_for_use/Fluxnet_Data/Download_data/"
  # sel_site<-"CH-Oe1"
  # exdir.path<-'C:/Users/yluo/Desktop/CES/Data_for_use/Fluxnet_Data/Preprocessed_data/Unzip_ori_HH_data'
  
  #1).set the working directory to ori.zip.path
  setwd(ori.zip.path)
  all.files.name<-list.files()
  
  #2).find the zip file going to be unzipped
  pos_site_zip<-grep(sel_site,all.files.name)
  zip.file.name<-all.files.name[pos_site_zip]
  
  #3). only unzip the half-hourly data in the zip file
  files.name.inzip<-unzip(zip.file.name,list=TRUE)
  #select half-hourly(HH) "FULLSET" or halfly(HR) csv
  pos_FULL_HH<-c(grep("FULLSET_HH",files.name.inzip$Name),grep("FULLSET_HR",files.name.inzip$Name))
  #only unzip the "FULLSET" csv file-->need to very careful with write the right directory of"zipfile" that going to be unzipped
  #also here exdir should end as".." rather than "../"
  unzip(zipfile=paste0(zip.file.name),files = files.name.inzip$Name[pos_FULL_HH],exdir=exdir.path)  
}
#
ori.zip.path<-"D:/CES/Data_for_use/Fluxnet_Data/Download_data/"
exdir.path<-'D:/CES/Data_for_use/Fluxnet_Data/Preprocessed_data/Unzip_ori_HH_data'

### for the selected sites for the analysis:
sites.path<-"./data-raw/raw_data/sites_info/"
load(paste0(sites.path,"Pre_selected_sites_info.RDA"))
#unzip the files-->have already done this before
# sites_goingto_processed<-df_sites_sel$sitename
# for(i in 1:length(sites_goingto_processed)){
#   unzip_file(ori.zip.path,sites_goingto_processed[i],exdir.path)
# }

#----------------------------------------------------
#b. selected the variables we are interested in:
#----------------------------------------------------
library(lubridate)
library(tidyverse)
#source the function based on Beni:
fun.path<-"./R/"
source(paste0(fun.path,"remove_outliers.R"))
source(paste0(fun.path,"clean_fluxnet_gpp.R"))
#----------------
#b1. write a function to subset the interested variables
#----------------
subset_interest_vars<-function(proc.path,sel_site){
  # proc.path<-'D:/CES/Data_for_use/Fluxnet_Data/Preprocessed_data/Unzip_ori_HH_data'
  # sel_site<-"DK-Sor"
  
  setwd(proc.path)
  csv.names<-list.files()
  #
  pos_sel_site<-grep(sel_site,csv.names)
  sel_site_csv<-csv.names[pos_sel_site]
  #read the sel csv file
  df<-read.csv2(file=sel_site_csv,sep = ",")
  #selecting the variables that we are interested in:
  #refer the fluxnet variables description: https://fluxnet.org/data/fluxnet2015-dataset/fullset-data-product/
  #select following variables:
  #TA_F(TA_F_QC) -->degC
  #SW_IN_F(SW_IN_F_QC)-->W m-2
  #PA_F(PA_F_QC)-->atmospheric pressure:kPa
  #P_F(P_F_QC)-->precipitaiton:mm
  #WS_F(WS_F_QC)-->wind speed:m s-1
  #PPFD_IN/PPFD_OUT(PPFD_IN_QC)-->photosynthetic photon flux denstiy: umol m-2 s-1
  #NEE_VUT_REF(NEE_VUT_REF_QC)-->umolCO2 m-2 s-1
  #GPP_NT/DT_VUT_REF-->VUT:variable u*threshold REF-->reference GPP, using model efficiency approach
  #TS_F_MDS_#(TS_F_MDS_#_QC)-->soil temperature, # stands for soil depth
  #SWC_F_MDS_#(SWC_F_MDS_#_QC)-->soil water content (%)
  
  #
  vars_names<-names(df)
  #
  pos_TIMESTAMP<-grep("TIMESTAMP",vars_names)
  pos_TA<-c(match("TA_F",vars_names),match("TA_F_QC",vars_names))
  pos_SW_IN<-c(match("SW_IN_F",vars_names),match("SW_IN_F_QC",vars_names))
  pos_SW_OUT<-c(grep("SW_OUT",vars_names))
  pos_PA_F<-grep("PA_F",vars_names)
  pos_P_F<-grep("P_F",vars_names)
  pos_WS_F<-grep("WS_F",vars_names)
  pos_PPFD<-c(grep("PPFD_IN",vars_names),grep("PPFD_OUT",vars_names))
  pos_VPD<-grep("VPD_F",vars_names)
  pos_NEE_VUT_REF<-c(match("NEE_VUT_REF",vars_names),match("NEE_VUT_REF_QC",vars_names))
  pos_GPP_VUT_REF<-c(grep("GPP_NT_VUT_REF",vars_names),grep("GPP_DT_VUT_REF",vars_names))
  pos_TS_F_MDS<-grep("TS_F",vars_names)
  pos_SWC_F_MDS<-grep("SWC_F",vars_names)
  #sel vars position:
  pos_all<-c(pos_TIMESTAMP,pos_TA,pos_SW_IN,pos_SW_OUT,pos_PA_F,pos_P_F,pos_WS_F,pos_PPFD,pos_VPD,
             pos_NEE_VUT_REF,pos_GPP_VUT_REF,pos_TS_F_MDS,pos_SWC_F_MDS)
  df_sel<-df[,pos_all]
  #tidy the format for the half hourly data:
  #nrow(df_sel[df_sel$NEE_VUT_REF_QC==0|df_sel$NEE_VUT_REF_QC==1,])/nrow(df_sel)
  #assign -9999=NA
  df_sel[df_sel==c(-9999)]<-NA
  #time format:
  df_sel$TIMESTAMP_START<-lubridate::ymd_hm(df_sel$TIMESTAMP_START)
  df_sel$TIMESTAMP_END<-lubridate::ymd_hm(df_sel$TIMESTAMP_END)
  #format of other variables-->other variables are numeric 
  df_sel[,3:ncol(df_sel)]<-apply(df_sel[,3:ncol(df_sel)],2,as.numeric)
  #only keep the GPP and NEE value when NEE_VUT_REF_QC==0 or ==1, others set to NA
  N_NA<-nrow(df_sel[df_sel$NEE_VUT_REF_QC>1 & !is.na(df_sel$NEE_VUT_REF_QC),])
  df_sel[df_sel$NEE_VUT_REF_QC>1&!is.na(df_sel$NEE_VUT_REF_QC),]$GPP_DT_VUT_REF<-rep(NA,N_NA)
  df_sel[df_sel$NEE_VUT_REF_QC>1&!is.na(df_sel$NEE_VUT_REF_QC),]$GPP_NT_VUT_REF<-rep(NA,N_NA)
  df_sel[df_sel$NEE_VUT_REF_QC>1&!is.na(df_sel$NEE_VUT_REF_QC),]$NEE_VUT_REF<-rep(NA,N_NA)
  #further removing the flux data accorrding to Stocker et al., 2020(https://gmd.copernicus.org/articles/13/1545/2020/#bib1.bibx178)
  #and Beni's code (https://github.com/stineb/ingestr/blob/master/R/get_obs_bysite_fluxnet.R)
  # using clean_fluxnet_gpp <- function(df, nam_gpp_nt, nam_gpp_dt, nam_nt_qc, nam_dt_qc, threshold, remove_neg = FALSE, filter_ntdt){
  # for "filter_ntdt" parameter -->Remove data points where the two flux decompositions are inconsistent,
  ## i.e. where the residual of their regression is above the 97.5% or below the 2.5% quantile.
  df_sel_new<-clean_fluxnet_gpp(df_sel,remove_neg = FALSE,filter_ntdt = T)
  #return the results:
  return(df_sel_new)
}
#----------------
#b2. put the selected variables from different sites into a data.frame
#----------------
library(dplyr)
proc.path<-'D:/CES/Data_for_use/Fluxnet_Data/Preprocessed_data/Unzip_ori_HH_data'
###II.for the available analyzed sites in FLUXNET2015:
df_all<-c()
for(i in 1:nrow(df_sites_sel)){
  df.temp<-subset_interest_vars(proc.path = proc.path,df_sites_sel$sitename[i])  
  df.temp<-data.frame(sitename=rep(df_sites_sel$sitename[i],nrow(df.temp)),df.temp)
  #using dplyr::bind_rows(df1, df2) to bind two sites that might have different variable numbers
  if(i==1){
    df_all<-df.temp
  }
  if(i>1){
    df_all<-bind_rows(df_all,df.temp)
  }
  rm(df.temp)
  print(i)
}

#----------------
#b3.save the preprocessed selected HH data
#----------------
#some unit conversion: convert VPD (hPa)-->to VPD (Pa)
df_all$VPD_F<-df_all$VPD_F*100
df_all$VPD_F_MDS<-df_all$VPD_F_MDS*100
#
#the sites according Beni' datasets-->from Fluxnet2015
setwd("D:/Github/photocold_manuscript/")
save.path<-"./data-raw/raw_data/processed_data_from_FLUX2015/"
#for the other available sites in FLUXNET2015:
save(df_all,file=paste0(save.path,"HH_data.RDA"))

#----------------------------------------------------
#c.summary the HH data to daily data:
#----------------------------------------------------
library(plyr)
#----------------
#c1.summary half-hourly to daily data for each site
#----------------
#I.for the sites accoring to the Fluxnet2015 tidy by Beni:
#first to select the variables interested:
sel_variables<-c("sitename","TIMESTAMP_START","TA_F","SW_IN_F","SW_OUT",
                 "PA_F","P_F","WS_F","PPFD_IN","PPFD_OUT","VPD_F",  
                 "NEE_VUT_REF",
                 "GPP_NT_VUT_REF","GPP_DT_VUT_REF",
                 paste0("TS_F_MDS_",c(1:7)),paste0("SWC_F_MDS_",c(1:7))) #from the datasets-->most layer of SWC =7, and TS=7
df_all_sel<-df_all[,sel_variables]
##Working here!
#merge to daily
df_all_sel$Date<-format(df_all_sel$TIMESTAMP_START,format = "%Y-%m-%d")
df_all_sel$HH<-hour(df_all_sel$TIMESTAMP_START)
#summarize the data to daily 
#-for VPD,only using the data in the day time (HH>=6 <=18)-->set the VPD value == NA for non-day 
df_all_sel[df_all_sel$HH<6 | df_all_sel$HH>18,]$VPD_F<-NA
#-for SW_IN, SW_OUT, ppfd_IN, and ppfd_OUT-->calculated daily mean, daily midday(10-14) mean, and daily midday max
df_all_sel_Rg<-df_all_sel[,c("SW_IN_F","SW_OUT","PPFD_IN","PPFD_OUT","HH")]
df_all_sel_Rg[df_all_sel$HH<10 | df_all_sel$HH>14,]<-NA
names(df_all_sel_Rg)<-c(paste0(c("SW_IN_F","SW_OUT","PPFD_IN","PPFD_OUT"),"_midday"),"HH")
#
df_all_sel<-cbind(df_all_sel,df_all_sel_Rg[,c(1:4)]) #do not add "HH"


df_all_sel_daily<-plyr::ddply(df_all_sel,.(sitename,Date),summarise,
  Ta_mean=mean(TA_F,na.rm = T),TA_min=min(TA_F,na.rm = T),TA_max=max(TA_F,na.rm = T),
  SW_IN_fullday_mean=mean(SW_IN_F,na.rm = T),SW_IN_midday_mean=mean(SW_IN_F_midday,na.rm=T),SW_IN_midday_max=max(SW_IN_F_midday,na.rm = T),
  SW_OUT_fullday_mean=mean(SW_OUT,na.rm=T),SW_OUT_midday_mean=mean(SW_OUT_midday,na.rm=T),SW_OUT_midday_max=max(SW_OUT_midday,na.rm = T),
  PPFD_IN_fullday_mean=mean(PPFD_IN,na.rm=T),PPFD_IN_midday_mean=mean(PPFD_IN_midday,na.rm=T),PPFD_IN_midday_max=max(PPFD_IN_midday,na.rm = T),
  PPFD_OUT_fullday_mean=mean(PPFD_OUT,na.rm=T),PPFD_OUT_midday_mean=mean(PPFD_OUT_midday,na.rm=T),PPFD_OUT_midday_max=max(PPFD_OUT_midday,na.rm = T),
  PA_mean=mean(PA_F,na.rm = T),P=sum(P_F,na.rm = T),
  WS_mean=mean(WS_F,na.rm = T),
  VPD_day_mean=mean(VPD_F,na.rm=T),
  NEE_mean=mean(NEE_VUT_REF,na.rm=T),
  GPP_NT_mean=mean(GPP_NT_VUT_REF,na.rm = T),GPP_DT_mean=mean(GPP_DT_VUT_REF,na.rm = T),
  TS_1_mean=mean(TS_F_MDS_1,na.rm = T),TS_2_mean=mean(TS_F_MDS_2,na.rm = T),
  TS_3_mean=mean(TS_F_MDS_3,na.rm = T),TS_4_mean=mean(TS_F_MDS_4,na.rm = T),
  TS_5_mean=mean(TS_F_MDS_5,na.rm = T),TS_6_mean=mean(TS_F_MDS_6,na.rm = T),
  TS_7_mean=mean(TS_F_MDS_7,na.rm = T),
  SWC_1_mean=mean(SWC_F_MDS_1,na.rm = T),SWC_2_mean=mean(SWC_F_MDS_2,na.rm = T),
  SWC_3_mean=mean(SWC_F_MDS_3,na.rm = T),SWC_4_mean=mean(SWC_F_MDS_4,na.rm = T),
  SWC_5_mean=mean(SWC_F_MDS_5,na.rm = T),SWC_6_mean=mean(SWC_F_MDS_6,na.rm = T),
  SWC_7_mean=mean(SWC_F_MDS_7,na.rm = T)
)
#checke the non NAs in each variable
apply(df_all_sel_daily[,-c(1:2)],2,function(x){sum(!is.na(x))})

#----------------
#c2.save the preprocessed daily data
#----------------
#I.for the sites accoring to the data sent by Beni:
save.path<-"./data-raw/raw_data/processed_data_from_FLUX2015/"
save(df_all_sel_daily,file=paste0(save.path,"Daily_data.RDA"))

#################################additional code############################################
#----------------------------------------------------
#d.comparing the processed results with provided data by Beni
#-->this analysis is mainly for the sites according to the data sent by Beni:
#----------------------------------------------------
#load the Beni's data
load.path<-"D:/CES/Data_for_use/Data_sent_by_Beni/"
df_Beni<-read.csv(file=paste0(load.path,"ddf_fluxnet2015_pmodel_with_forcings_stocker19gmd.csv"))
#finding the corresponding sites between "df_Beni"and "df_all_sel_daily"
sel_sites<-intersect(unique(df_all_sel_daily$sitename),unique(df_Beni$sitename))
#
df_Beni_sel<-c()
for(i in 1:length(sel_sites)){
  df_temp<-subset(df_Beni,df_Beni$sitename==sel_sites[i])
  df_Beni_sel<-rbind(df_Beni_sel,df_temp)
}
df_Beni_sel$date<-mdy(df_Beni_sel$date)
#----------------
#d1.write a function to compare two datasets
#----------------
library(ggplot2)
compare_vars<-function(df_YP_daily,df_Beni_daily,var_in_YP,var_in_Beni){
  # df_YP_daily<-df_all_daily
  # df_Beni_daily<-df_Beni_sel
  # var_in_YP<-"GPP_NT_mean"
  # var_in_Beni<-"gpp_obs"
  
  #
  df_YP_daily_sel<-df_YP_daily[,c("sitename","Date",var_in_YP)];names(df_YP_daily_sel)<-c("sitename","date",var_in_YP)
  df_Beni_daily_sel<-df_Beni_daily[,c("sitename","date",var_in_Beni)]
  #remove the NA
  df_YP_daily_sel<-df_YP_daily_sel[!is.na(df_YP_daily_sel[,3]),]
  df_Beni_daily_sel<-df_Beni_daily_sel[!is.na(df_Beni_daily_sel[,3]),]
  ####merge two datasets:
  #make sure the format of variables is the same
  df_YP_daily_sel$date<-ymd(df_YP_daily_sel$date)
  df_merge<-merge(df_YP_daily_sel,df_Beni_daily_sel,by = c("sitename","date"),all.x = T)
  df_merge<-df_merge[!is.na(df_merge[,3])&!is.na(df_merge[,4]),]
  names(df_merge)<-c("sitename","date","tidy_YP","tidy_Beni")
  #plotting
  p_out<-ggplot(data=df_merge,aes(x=tidy_YP,y=tidy_Beni,col=sitename))+
    geom_point()+
    xlab(paste0("tidy_YP: ",var_in_YP))+
    ylab(paste0("tidy_Beni: ",var_in_Beni))
  # print(p_out)
  return(p_out)
}
####
#before comparing-->convert the unit:
#convert the unit of GPP from umol m-2 s-1 to gC m-2 d-1
df_all_sel_daily$GPP_NT_mean<-df_all_sel_daily$GPP_NT_mean*24*3600*12/1000000
df_all_sel_daily$GPP_DT_mean<-df_all_sel_daily$GPP_DT_mean*24*3600*12/1000000
#for GPP
p_gpp<-compare_vars(df_all_sel_daily,df_Beni_sel,"GPP_NT_mean","gpp_obs")
p_ppfd<-compare_vars(df_all_sel_daily,df_Beni_sel,"PPFD_IN_mean","ppfd_fluxnet2015")
p_temp<-compare_vars(df_all_sel_daily,df_Beni_sel,"Ta_mean","temp_day_fluxnet2015")
p_preci<-compare_vars(df_all_sel_daily,df_Beni_sel,"P","prec_fluxnet2015")
p_patm<-compare_vars(df_all_sel_daily,df_Beni_sel,"PA_mean","patm_fluxnet2015")
p_vpd<-compare_vars(df_all_sel_daily,df_Beni_sel,"VPD_day_mean","vpd_day_fluxnet2015")

library(cowplot)
p_merge<-plot_grid(p_gpp,p_ppfd,p_temp,p_preci,p_patm,p_vpd,
          labels = "auto",ncol=2,label_size = 12,align = "hv")
#save the plot
# save.path<-"C:/Users/yluo/Desktop/R_testcode/PhotoCold/Second_round_of_code/plot/comp_Beni_YP/"
# ggsave(paste0(save.path,"p_merge_update_Aug12.png"),p_merge,width = 15,height = 12)
