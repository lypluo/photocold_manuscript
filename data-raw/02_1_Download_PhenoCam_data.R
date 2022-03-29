##########################################
#download the PhenoCam data from PhenoCam US and EuroPhen:
# PhenoCam US: https://phenocam.sr.unh.edu/webcam/ and http://europhen.org/
#data of European phenocam network are from Koen
##########################################
library(phenocamapi)
#----------------------------------------
#0.display the information of the PhenoCam infromation for the EC sites I am analyzing:
#----------------------------------------
library(readxl)
load.path<-"D:/CES/Data_for_use/Merge_Data/ECflux_and_PhenoCam_site_info/"
Meta.info<-read.csv(paste0(load.path,"ECflux_and_PhenoCam_site_info_add_manually.csv"))
#----------------------------------------
#1.download the data from the PhenoCam USA/ and PhenoCam data for european sites from Koen
#----------------------------------------
Meta.info<-as.data.frame(Meta.info)
#
Meta.info_PhenoCam_US<-Meta.info[!is.na(Meta.info$Cam_Source)&Meta.info$Cam_Source=="PhenoCam USA",]
#
Meta.info_PhenoCam_EU<-Meta.info[!is.na(Meta.info$Cam_Source) & Meta.info$Cam_Source=="EuroPhen",]
#---------------------
#1a.For PhenoCam-US sites:
#---------------------
#download the data through the R package phenocamapi
library(phenocamapi)
#write the function to download the PhenoCam data
ts_download<-function(Site,VegType,RoiID,Type){
  # Site="sylvania"
  # VegType="DB"
  # RoiID=c(1000)
  # Type="1day"
  
  ##download the data
  data_df<-c()
  df_names<-paste0(VegType,"_",RoiID)
  for(i in 1:length(RoiID)){
    data_temp<-get_pheno_ts(site=Site,vegType = VegType,roiID = RoiID[i],type = Type)
    #add PFT and ROI_name
    data_temp$PFT<-rep(VegType,nrow(data_temp))
    data_temp$RoiID<-RoiID[i]
    data_df<-rbind(data_df,data_temp)
  }
  return(data_df)
}
#start to downloading the phenocam data
# PhenoCam_USA_sitename<-Meta.info_PhenoCam_US$Cam_sitename
#
save.path<-"D:/CES/Data_for_use/PhenoCam_Data/Through_PhenoCam_USA/VIs/"
#download 1-day time resolution data
site_names<-Meta.info_PhenoCam_US[,c("SiteName","Cam_sitename","PFT")]
#And I have selected the most representative ROIs
#--------
#for the sites from Beni's datasets
#--------
#US-Syv
Cam_US_Syv_DB<-ts_download("sylvania","DB",c(1000),"1day")
Cam_US_Syv_EN<-ts_download("sylvania","EN",c(1001),"1day")
Cam_US_Syv<-rbind(Cam_US_Syv_DB,Cam_US_Syv_EN)
#US-UMB
Cam_US_UMB<-ts_download("umichbiological","DB",c(1000,2000),"1day")
#US-UMd
Cam_US_UMd<-ts_download("umichbiological2","DB",c(1000,2000),"1day")
#US-WCr
Cam_US_WCr<-ts_download("willowcreek","DB",c(1000),"1day")
#CA-Qfo
Cam_CA_Qfo<-ts_download("chibougamau","EN",c(1000),"1day")
#IT-Tor
Cam_IT_Tor<-ts_download("torgnon-nd","GR",c(1000),"1day")
#--------
#for the sites beyond Beni's datasets
#--------
#AT-Neu:
Cam_AT_Neu<-ts_download("innsbruck","GR",c(1000),"1day")
#CA-Gro:
Cam_CA_Gro_DB<-ts_download("groundhog","DB",c(1000,2000),"1day")  #combine two ROIs to get full period of data
Cam_CA_Gro_EN<-ts_download("groundhog","EN",c(1000,2000),"1day")  #combine two ROIs to get full period of data
Cam_CA_Gro<-rbind(Cam_CA_Gro_DB,Cam_CA_Gro_EN)
#CA-Oas
Cam_CA_Oas<-ts_download("canadaOA","DB",c(1000),"1day")
#CA-Obs
Cam_CA_Obs_DN<-ts_download("canadaOBS","DN",c(1000),"1day")
Cam_CA_Obs_EN<-ts_download("canadaOBS","EN",c(1000),"1day")
Cam_CA_Obs<-rbind(Cam_CA_Obs_DN,Cam_CA_Obs_EN)
#CA-TP1
Cam_CA_TP1<-ts_download("turkeypointenf02","EN",c(1000),"1day")
#CA-TP3
Cam_CA_TP3<-ts_download("turkeypointenf74","EN",c(1000),"1day")
#CA-TP4
Cam_CA_TP4<-ts_download("turkeypointenf39","EN",c(1000),"1day")
#IT-MBo
Cam_IT_MBo<-ts_download("montebondonegrass","GR",c(1000),"1day")
#US-Ha1
Cam_US_Ha1<-ts_download("bbc1","DB",c(1000),"1day") #select bbc1 at the end (bbc2 is similar)
#US-MMS
Cam_US_MMS<-ts_download("morganmonroe","DB",c(1000),"1day") 
#US-NR1
Cam_US_NR1<-ts_download("niwot2","EN",c(1000),"1day") 
#--------------
#save the data
#--------------
df.PhenoCam_USA<-list(Cam_AT_Neu,Cam_CA_Gro,Cam_CA_Oas,Cam_CA_Obs,Cam_CA_TP1,Cam_CA_TP3,Cam_CA_TP4,
                      Cam_IT_MBo,Cam_US_Ha1,Cam_US_MMS,Cam_US_NR1,
                      Cam_US_Syv,Cam_US_UMB,Cam_US_UMd,Cam_US_WCr,Cam_CA_Qfo,Cam_IT_Tor)
for(i in 1: length(df.PhenoCam_USA)){
  write.csv(df.PhenoCam_USA[[i]],file = paste0(save.path,"Cam_",site_names$SiteName[i],".csv"))
}

#---------------------
#1b.For PhenoCam-EU sites:
#---------------------
# first to unzip the tar.gz files-->has finished, do not need to do anymore as I receive data from Koen
# see below code
# path<-"C:/Users/yluo/Desktop/CES/Data_for_use/PhenoCam_Data/Through_EuroPheno/"
# setwd(path)
# untar("phenocam_2018.tar.gz")
# untar("validation_product_2020.tar.gz")
#refer the website:http://european-webcam-network.net/
#at the end, only find four sites data are aviable-->they are: DE-Hai, FI-Hyy, BE-Vie,and DE_Tha
load.path<-"D:/CES/Data_for_use/PhenoCam_Data/Through_EuroPheno/Ori_Data/validation_product/"
#the data from EuroPheno is less tidy as PhenoCam-Us, needs more data quality check...
ts_tidy<-function(Data.path,Site,VegType,RoiID,Type){
  # Data.path<-load.path
  # Site="DE-Hai"
  # VegType="DB"
  # RoiID=c("0001")
  # Type="1day"
  
  ##download the data
  data_df<-c()
  df_names<-paste0(VegType,"_",RoiID)

    for(i in 1:length(RoiID)){
      folder_name<-paste0(Site,"_",df_names)
      #
      if(Site=="DE-Hai"|Site=="FI-Hyy"|Site=="DE-Lei"){
      file_name<-paste0(folder_name,"_1day.csv")
      data_temp<-read.csv(file=paste0(Data.path,folder_name,"/",file_name),skip = 22)
      }
      if(Site=="BE-Vie"){
        file_name<-paste0(folder_name,"_3day.csv")
        data_temp<-read.csv(file=paste0(Data.path,folder_name,"/",file_name),skip = 24)
      }
      if(Site=="DE-Tha"|Site=="DK-Sor"){
        file_name<-paste0(folder_name,"_1day.csv")
        data_temp<-read.csv(file=paste0(Data.path,folder_name,"/",file_name),skip = 24)
      }
      #add PFT and ROI_name
      data_temp$PFT<-rep(VegType,nrow(data_temp))
      data_temp$RoiID<-RoiID[i]
      data_df<-rbind(data_df,data_temp)
    }

  return(data_df)
}
#
site_names<-Meta.info_PhenoCam_EU[,c("SiteName","Cam_sitename","PFT")]
#--------
#for the sites from Beni's datasets
#--------
#for DE-Hai
Cam_DE_Hai<-ts_tidy(load.path,"DE-Hai","DB","0001","1day")
#for FI-Hyy
Cam_FI_Hyy<-ts_tidy(load.path,"FI-Hyy","EN","1000","1day")
#for BE_Vie
Cam_BE_Vie<-ts_tidy(load.path,"BE-Vie","EN","1000","1day")
#for DE_Tha
Cam_DE_Tha<-ts_tidy(load.path,"DE-Tha","EN","0001","1day")
#--------
#for the sites beyond Beni's datasets
#--------
Cam_DE_Lnf<-ts_tidy(load.path,"DE-Lei","DB","0001","1day")
Cam_DK_Sor<-ts_tidy(load.path,"DK-Sor","DB","0001","1day")

#save the data:
save.path<-"D:/CES/Data_for_use/PhenoCam_Data/Through_EuroPheno/VIs/"
df.PhenoCam_EU<-list(Cam_DE_Lnf,Cam_DK_Sor,Cam_BE_Vie,Cam_DE_Hai,Cam_DE_Tha,Cam_FI_Hyy)
site_names_EU<-c("DE-Lnf","DK-Sor","BE-Vie","DE-Hai","DE-Tha","FI-Hyy")
for(i in 1: length(df.PhenoCam_EU)){
  write.csv(df.PhenoCam_EU[[i]],file = paste0(save.path,"Cam_",site_names_EU[i],".csv"))
}

#----------------------------------------
#2.plot the selected ROIs-->update on Aug,9th:at this moment-->do not plot the ROIs-->will update in the future.
#----------------------------------------
#---------------------
#For PhenoCam-US sites:
#---------------------
library(tiff)
#information for all the rois
all_rois<-get_rois()
#a. download the selected sites ROI Masks and ROIs
ROI.path<-"D:/CES/Data_for_use/PhenoCam_Data/Through_PhenoCam_USA/ROIs/"
create_ROI_tif_folder<-function(all_rois,download_site,VegType,RoiID,save.folder){
  # download_site<-site_names$Cam_sitename[1]
  # all_rois<-all_rois
  # VegType<-"DB"
  # RoiID<-1000
  # save.folder<-ROI.path
  #
  temp_sel_rois<-all_rois[all_rois$site==download.site,]
  sel_rois<-temp_sel_rois[temp_sel_rois$roitype==VegType&temp_sel_rois$sequence_number==RoiID,]
  roi.url<-paste0(sel_rois$roi_page,sel_rois$roi_name,"_01.tif")
  #down load the tif
  # tmp <- tempfile(fileext = ".tif")
  setwd(save.folder)
  dir.create(download_site)
  #!!commented by YP: at the moment-->the downloaded tif does not in a right form-->hence only use this function to creat folder
  # and manually download the tif from websites:
  # save.path<-paste0(save.folder,download.site,"/",VegType,"_",RoiID,".tif")
  # if (class(try(download.file(roi.url, destfile = save.path, mode = "w"))) == 
  #     "try-error") {
  #   stop("file was not found on the server!")
  # }
}
#starting to downloading:
#US-Syv,#US-UMB.#US-UMd.#US-WCr,CA-Qfo,IT-Tor
#actually, do not use Vegetation type and ROIID in this function...
create_ROI_tif_folder(all_rois,site_names$Cam_sitename[1],"DB",1000,ROI.path)
create_ROI_tif_folder(all_rois,site_names$Cam_sitename[2],"DB",1000,ROI.path)
create_ROI_tif_folder(all_rois,site_names$Cam_sitename[3],"DB",1000,ROI.path)
create_ROI_tif_folder(all_rois,site_names$Cam_sitename[4],"DB",1000,ROI.path)
create_ROI_tif_folder(all_rois,site_names$Cam_sitename[5],"EN",1000,ROI.path)
create_ROI_tif_folder(all_rois,site_names$Cam_sitename[6],"GR",1000,ROI.path)
##working here-->download the tif and background manually
##in the future-->going to plot the selected ROIs in the images

#---------------------
#For PhenoCam-EU sites:
#---------------------
#-->the data obtained from Koen