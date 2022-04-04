###############################################
#tidy the PhenoCam VIs and MODIS VIs from the interested sites:
##############################################
#From previous analysis and compare between the mean T between "non-overestimation"
#and "overestimation" sites in different PFTs-->summarize as below:
#for DBF:put US-UMB, US-UMd, US-WCr, CA-Qfo as sites under low-Ta stress; while
#        put DE-Hai, DK-Sor, and DE-Lnf as sites with moderate or no Ta stress

#for ENF:put CA-Qfo,US-NR1,FI-Hyy,CA-Obs as sites under low-Ta strees, while
#        put DE-Tha,CA-TP1,CA-TP3, CA-TP4 as sites with moderate and no Ta stress
#the GPP overestimation period are as follows for differnt PFTs:
# DBF: 70-103
# ENF: 81-143
# MF: 69-147
# for all PFTs range: 69-167
#!!In this study, I would like to compare the VIs differences betweeen "overestimated"
#and "non-overestimated" site:
#---------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
sel_sites<-data.frame(sitename=c("US-UMB","US-UMd","US-WCr","CA-Gro","DE-Hai","DK-Sor","DE-Lnf",
                                 "CA-Qfo","US-NR1","FI-Hyy","CA-Obs","DE-Tha","CA-TP1","CA-TP3","CA-TP4"),
                      PFT=c(rep("DBF",7),rep("ENF",8)),
                      flag=c(rep("Ta_stress",4),rep("noTa_stress",3),rep("Ta_stress",4),rep("noTa_stress",4)))
#----------------------------------
#(1) tidy PhenoCam VIs from interested sites
#----------------------------------
Cam.path<-"./data-raw/raw_data/processed_data_from_PhenoCam/"
load(paste0(Cam.path,"Daily_data.RDA"))
df.Cam_sel<-df.Phenocam_daily %>%
  filter(sitename %in% sel_sites$sitename)
df.Cam_sel$date<-as.Date(df.Cam_sel$date)

#----------------------------------
#(2)tidy the MODIS VIs from interested sites-->using data from Walther et al., 2021
#----------------------------------
#A.manually selected the nc files from the interested sites:
#B.open the nc files 
library(ncdf4)
library(raster)
library(lubridate)
nc.files.path<-"D:/data/FallPheno_project/MODIS_VIs_from_Walther2021/nc_files/"
#
sites<-sel_sites$sitename
df.VIs<-c()
for(i in 1:length(sites)){
  fname<-list.files(path = nc.files.path, pattern=sites[i], full.names = TRUE)
  sitename<-sites[i]
  if (length(fname) > 0){
    print(paste('---> Analyzing ', sitename))
    nc <- nc_open(fname)
    
    ##a. extract the variables from the nc file
    #Date-->since 1970-01-01
    N<-nc$dim$time$len
    Date<-as.Date(nc$dim$time$vals,origin="1970-01-01")
    
    #VIs-->BLUE,GREEN,RED,NIR,EVI,NDVI,NIRv,generalized NDVI(kNDVI)
    b_blue<-ncvar_get(nc,"BLUE")
    b_green<-ncvar_get(nc,"GREEN")
    b_red<-ncvar_get(nc,"RED")
    b_nir<-ncvar_get(nc,"NIR")
    EVI<-ncvar_get(nc,"EVI")
    NDVI<-ncvar_get(nc,"NDVI")
    NIRv<-ncvar_get(nc,"NIRv")
    kNDVI<-ncvar_get(nc,"kNDVI")
    #calculating GRVI:
    GRVI_MODIS<-c(b_green - b_red)/c(b_green + b_red)
    ##b.create the df
    #produce the new datasets
    df.temp<-data.frame(sitename=rep(sitename,N),date=Date,b_blue=b_blue,b_green=b_green,b_red=b_red,
                        b_nir=b_nir,EVI=EVI,NDVI=NDVI,NIRv=NIRv,kNDVI=kNDVI,GRVI_MODIS=GRVI_MODIS)
    df.VIs<-rbind(df.VIs,df.temp)
  }
}
df.MODIS_sel<-df.VIs
##testing
# df.VIs %>% filter(sitename=="DE-Hai") %>%
#   ggplot(.,aes(x=Date,y=EVI))+
#   geom_point()

#-------------------------------------------
#(3)merge the PhenoCam VIs and MODIS VIs and compare the differences:
#-------------------------------------------
df.VIs<-merge(df.MODIS_sel,df.Cam_sel,by=c("sitename","date"),all.x=T)
#test
df.VIs$sitedate<-paste0(df.VIs$sitename,df.VIs$date)
df.VIs<-df.VIs[!duplicated(df.VIs$sitedate),] #remove the duplicated values
df.VI_new<-left_join(df.VIs,sel_sites,by="sitename")
df.VI_new<-df.VI_new %>%
  mutate(sitedate = NULL)
df.VI_new<-df.VI_new %>%
  mutate(month=month(date))
#-----------------I.exploring--------------------------
#-----------
#A. For the DBF
#----------
#check the time series:
df.VI_new %>%
  filter(PFT=="DBF") %>%
  group_by(sitename,doy)%>%
  dplyr::summarise(doy=mean(doy),mean_gcc=mean(gcc_90,na.rm=T),
                   mean_EVI=mean(EVI,na.rm=T),mean_MODIS_GRVI=mean(GRVI_MODIS,na.rm=T))%>%
  ggplot()+
    # geom_point(aes(x=doy,y=mean_gcc,col="gcc"))+
    # geom_point(aes(x=doy,y=mean_EVI,col="EVI"))+
    geom_point(aes(x=doy,y=mean_MODIS_GRVI,col="MODIS_GRVI"))+
    facet_grid(~sitename)
df.DBF<-df.VI_new %>%
  filter(PFT=="DBF")%>%
  filter(month<=6)%>%
  group_by(flag,doy)%>%
  dplyr::summarise(GCC=mean(gcc_90,na.rm=T),
                   RCC=mean(rcc_90,na.rm=T),
                   GRVI=mean(GRVI,na.rm=T),
                   NDVI=mean(NDVI,na.rm=T),
                   EVI=mean(EVI,na.rm=T),
                   NIRv=mean(NIRv,na.rm=T),
                   kNDVI=mean(kNDVI,na.rm=T))
df.DBF%>%
  ggplot(aes(x=doy,y=GCC,col=flag))+
  geom_line()

#-----------
#B. For the ENF
#----------
#check the time series:
df.VI_new %>%
  filter(PFT=="ENF") %>%
  group_by(sitename,doy)%>%
  dplyr::summarise(doy=mean(doy),mean_gcc=mean(gcc_90,na.rm=T),
                   mean_EVI=mean(EVI,na.rm=T))%>%
  ggplot()+
  geom_point(aes(x=doy,y=mean_gcc,col="gcc"))+
  # geom_point(aes(x=doy,y=mean_EVI,col="EVI"))+
  facet_grid(~sitename)

df.ENF<-df.VI_new %>%
  filter(PFT=="ENF")%>%
  filter(month<=6)%>%
  group_by(flag,doy)%>%
  dplyr::summarise(GCC=mean(gcc_90,na.rm=T),
                   RCC=mean(rcc_90,na.rm=T),
                   GRVI=mean(GRVI,na.rm=T),
                   NDVI=mean(NDVI,na.rm=T),
                   EVI=mean(EVI,na.rm=T),
                   NIRv=mean(NIRv,na.rm=T),
                   kNDVI=mean(kNDVI,na.rm=T))
df.ENF%>%
  ggplot(aes(x=doy,y=GCC,col=flag))+
  geom_line()
#-----------------II.really analysis------------------------------
plot_VIs<-function(df,PFT_name,VI_source){
  # df<-df.VI_new
  # PFT_name<-"ENF"
  # VI_source<-"MODIS"
  
  #
  df.proc<-df %>%
    dplyr::select(sitename,date,EVI,NDVI,NIRv,kNDVI,GRVI_MODIS,
              gcc_90,rcc_90,GRVI,flag,year,month,doy,PFT)%>%
    filter(PFT==PFT_name)
  #test
  df.proc %>%
    group_by(sitename,doy)%>%
    dplyr::summarise(doy=mean(doy),GCC=mean(gcc_90,na.rm=T),
                     RCC=mean(rcc_90,na.rm=T),
                     GRVI=mean(GRVI,na.rm=T))%>%
    ggplot(aes(x=doy,y=GCC))+
    geom_point()+
    facet_grid(~sitename)
  
  #For PhenoCam VIs, in order to compare among sites, first to do the following steps:
  #for each sites using the 1) smooth(mean) the data according to doy:
  #2)VIs -VImin(smooth)-->(3)then normalize the VIs by ampl=VImax-VImin:
  Cam.VIs_stats_temp<-df.proc %>%
    group_by(sitename,doy)%>%
    dplyr::summarise(doy=mean(doy),gcc=mean(gcc_90,na.rm=T),
                     rcc=mean(rcc_90,na.rm=T),
                     GRVI=mean(GRVI,na.rm=T))
  Cam.VIs_stats<-Cam.VIs_stats_temp %>%
    dplyr::summarise(gcc_mn=min(gcc,na.rm = T),gcc_q10=quantile(gcc,probs=0.1,na.rm=T),gcc_mm=max(gcc,na.rm = T),
                     rcc_mn=min(rcc,na.rm = T),rcc_q10=quantile(rcc,probs=0.1,na.rm=T),rcc_mm=max(rcc,na.rm = T),
                     GRVI_mn=min(GRVI,na.rm = T),GRVI_q10=quantile(GRVI,probs=0.1,na.rm=T),GRVI_mm=max(GRVI,na.rm = T),
                     )%>%
    mutate(gcc_ampl=gcc_mm - gcc_q10,
           rcc_ampl=rcc_mm - rcc_q10,
           GRVI_ampl=GRVI_mm - GRVI_q10)
  #adjust the original PhenoCam VIs:
  df.proc_new<-left_join(df.proc,Cam.VIs_stats,by="sitename")
  df.results<-df.proc_new %>%
    mutate(gcc_adj= c(gcc_90 - gcc_q10)/gcc_ampl,
           rcc_adj= c(rcc_90 - rcc_q10)/rcc_ampl,
           GRVI_adj= c(GRVI - GRVI_q10)/GRVI_ampl)
  #################
  df.results$doy<-yday(df.results$date)
  df.results %>%
    group_by(sitename,doy)%>%
    dplyr::summarise(doy=mean(doy),GCC=mean(gcc_adj,na.rm=T),
                     RCC=mean(rcc_adj,na.rm=T),
                     GRVI=mean(GRVI_adj,na.rm=T))%>%
    ggplot(aes(x=doy,y=GCC))+
    geom_point()+
    facet_grid(~sitename)
  
  # if(PFT=="ENF"){
  #   df_summary<-df.results %>%
  #     filter(sitename!="DE-Tha")
  # }
  df_summary<-df.results %>%
      group_by(flag,doy)%>%
      # filter(month<=6)%>%
      dplyr::summarise(doy=mean(doy),GCC=mean(gcc_adj,na.rm=T),
                       RCC=mean(rcc_adj,na.rm=T),
                       GRVI=mean(GRVI_adj,na.rm=T),
                       NDVI=mean(NDVI,na.rm=T),
                       EVI=mean(EVI,na.rm=T),
                       NIRv=mean(NIRv,na.rm=T),
                       kNDVI=mean(kNDVI,na.rm=T),
                       GRVI_MODIS=mean(GRVI_MODIS,na.rm=T)
      )%>%
      pivot_longer(c(GCC,RCC,GRVI,NDVI,EVI,NIRv,kNDVI,GRVI_MODIS),names_to = "Source",values_to = "VIs")
    
 
    #plotting:
  if(VI_source=="PhenoCam"){
    p_plot<-df_summary %>%
      filter(Source %in% c("GCC","RCC","GRVI"))%>%
      # ggplot(aes(x=doy,y=GCC,col=flag))+
      # ggplot(aes(x=doy,y=RCC,col=flag))+
      # ggplot(aes(x=doy,y=GRVI,col=flag))+
      ggplot(aes(x=doy,y=VIs,col=flag))+
      geom_point()+
      annotate("rect",xmin=69,xmax=147,ymin = -Inf,ymax = Inf,alpha=0.2)+
      scale_color_manual("",values = c("noTa_stress"="green4","Ta_stress"="tomato"),
                         labels=c("Less Ta stress","Higher Ta stress"))+
      xlab("DoY (day)")+ylab("Normalized PhenoCam VIs")+
      facet_grid(~Source)+
      theme_classic()+
      theme(legend.text = element_text(size=20),
            axis.title = element_text(size=24),
            axis.text = element_text(size = 20),
            text = element_text(size=24))+
      annotate(geom = "text",x=200,y=-0.1,label=PFT_name,col="blue",size=8)
  }
  if(VI_source=="MODIS"){
    # New facet label names for Source variable-->change "GRVI_MODIS" to "GRVI"
    Source.labs <- c("EVI","GRVI","NDVI")
    names(Source.labs) <- c("EVI", "GRVI_MODIS","NDVI")
    
    p_plot<-df_summary %>%
      filter(Source %in% c("NDVI","EVI","GRVI_MODIS"))%>%
      ggplot(aes(x=doy,y=VIs,col=flag))+
      geom_point()+
      annotate("rect",xmin=69,xmax=147,ymin = -Inf,ymax = Inf,alpha=0.2)+
      scale_color_manual("",values = c("noTa_stress"="green4","Ta_stress"="tomato"),
                         labels=c("Less Ta stress","Higher Ta stress"))+
      xlab("DoY (day)")+ylab("Normalized MODIS VIs")+
      facet_grid(~Source)+
      theme_classic()+
      theme(legend.text = element_text(size=20),
            axis.title = element_text(size=24),
            axis.text = element_text(size = 20),
            text = element_text(size=24))+
      annotate(geom = "text",x=200,y=-0.1,label=PFT_name,col="blue",size=8)+
      facet_grid(
        . ~ Source, 
        labeller = labeller(Source = Source.labs)
      )
  }
  
  return(p_plot)
  
}
#
DBF_Cam_VIs<-plot_VIs(df.VI_new,"DBF","PhenoCam")+
  xlab("")
ENF_Cam_VIs<-plot_VIs(df.VI_new,"ENF","PhenoCam")
#
DBF_MODIS_VIs<-plot_VIs(df.VI_new,"DBF","MODIS")+
  xlab("")
ENF_MODIS_VIs<-plot_VIs(df.VI_new,"ENF","MODIS")

#merge the plot
library(ggpubr)
library(grid)
library(cowplot)
plot_Cam_final<-plot_grid(DBF_Cam_VIs,ENF_Cam_VIs,nrow = 2,labels = "AUTO")
plot_MODIS_final<-plot_grid(DBF_MODIS_VIs,ENF_MODIS_VIs,nrow = 2,labels = "AUTO")
#save the plots:
save.path<-"./manuscript/figures/"
ggsave(paste0(save.path,"Figure4_PhenoCam_VIs.png"),plot_Cam_final,width = 15,height = 10)
ggsave(paste0(save.path,"FigureS_MODIS_VIs.png"),plot_MODIS_final,width = 15,height = 10)

