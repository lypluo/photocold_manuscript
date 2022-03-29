##########################################
#Try to preprocess the PhenoCam data from Phenocam network
##########################################
#----------------------------------------------------
#1. load the data
#----------------------------------------------------
#---------------
#load the data from the PhenoCam-USA
#---------------
phenocam.path<-"D:/CES/Data_for_use/PhenoCam_Data/Through_PhenoCam_USA/VIs/"
setwd(phenocam.path)
site.files<-list.files()
site.files<-site.files[grep(".csv",site.files)]
#
df.Phenocam<-c()
for(i in 1:length(site.files)){
  temp<-read.csv(site.files[i])
  site.name<-substr(site.files[i],5,10)
  temp$sitename<-rep(site.name,nrow(temp))
  df.Phenocam<-rbind(df.Phenocam,temp)
  rm(temp)
}
#---------------
#load the data from the EuroPheno 
#---------------
phenocam.path<-"D:/CES/Data_for_use/PhenoCam_Data/Through_EuroPheno/VIs/"
setwd(phenocam.path)
site.files<-list.files()
df.Phenocam_EU<-c()
#
library(dplyr)
for(i in 1:length(site.files)){
  #using dplyr::bind_rows(df1, df2) to bind two sites that might have different variable numbers
  temp<-read.csv(site.files[i])
  site.name<-substr(site.files[i],5,10)
  temp$sitename<-rep(site.name,nrow(temp))
  #
  if(i==1){
    df.Phenocam_EU<-temp
  }
  #
  if(i>1){
  df.Phenocam_EU<-bind_rows(df.Phenocam_EU,temp)
  }
  rm(temp)
}
#---------------
#merge two data sources
#---------------
co_vars<-intersect(names(df.Phenocam),names(df.Phenocam_EU))
pos_1<-match(co_vars,names(df.Phenocam))
pos_2<-match(co_vars,names(df.Phenocam_EU))
df.Phenocam_new<-rbind(df.Phenocam[,pos_1],df.Phenocam_EU[,pos_2])
#----------------------------------------------------
#2. select the variables we are interested in and calculated the new PhenoCam VIs
#----------------------------------------------------
sel_variables<-c("sitename","date","year","doy","midday_r","midday_g","midday_b","gcc_90","rcc_90","max_solar_elev","snow_flag","PFT","RoiID")
df.Phenocam_daily<-df.Phenocam_new[,sel_variables]  #23 sites
names(df.Phenocam_daily)<-c("sitename","date","year","doy","midday_r","midday_g","midday_b","gcc_90","rcc_90","max_solar_elev","snow_flag","PhenoCam_PFT","RoiID")
#calculating GRVI:(DNgreen-DNred)/(DNgreen+DNred)-->here using (gcc-rcc)/(gcc+rcc)
df.Phenocam_daily$GRVI<-c(df.Phenocam_daily$gcc_90 - df.Phenocam_daily$rcc_90)/c(df.Phenocam_daily$gcc_90 + df.Phenocam_daily$rcc_90)
#----------------
#3.save the preprocessed daily data
#----------------
setwd("D:/Github/photocold_manuscript")
save.path<-"./data-raw/raw_data/processed_data_from_PhenoCam/"
save(df.Phenocam_daily,file=paste0(save.path,"Daily_data.RDA"))

