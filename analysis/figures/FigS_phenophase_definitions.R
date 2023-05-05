#####################################################################
#Aim:make the demonstration plot for the phenophases defintions:
#####################################################################
library(readr)
library(lubridate)
library(dplyr)
library(ggplot2)
library(tidyr)
library(cowplot)
library(phenopix)
library(zoo)

#-----------------------
#(1)load the data
#-----------------------
path<-"./data-raw/raw_data/Merged_data/"
load(paste0(path,"Merged_Flux_and_VIs.RDA"))
ddf<-df_merge;rm(df_merge)

#------------------------------
#(2)define the functions
#------------------------------
#function1: to smooth the time series
extract_smooth_y<-function(x,df){
  # x<-df_new$gpp_mod_FULL_f
  # df<-0.05
  x_sm<-SplineFit(x,df=df)
  y<-x_sm$fit$predicted
  return(y)
}

#funciton2:to determine the green-up period and basic phenophases(sos, eos)
#using the threshold method:
Deter_greenup<-function(df,pos_max){
  # df<-df_subset
  # pos_max<-pos_max
  
  df_sel_sos<-df %>%
    dplyr::filter(doy<=pos_max)
  #1) calculate the amplitude of the time series
  #based on the 10% or 25% of amplitude:
  mn_sos<-min(df_sel_sos$gpp_mod_norm_sm,na.rm = T)
  mm<-max(df_sel_sos$gpp_mod_norm_sm,na.rm = T)
  amp_sos<-mm-mn_sos
  
  trs_10_sos<-mn_sos+0.1*amp_sos
  trs_25_sos<-mn_sos+0.25*amp_sos
  
  #2)determine the start of season(sos10, and sos25)
  #set the sos beyond Feburary
  N_length<-length(df_sel_sos$gpp_mod_norm_sm)
  sos10<-which.min(abs(as.numeric(df_sel_sos$gpp_mod_norm_sm[60:N_length] - trs_10_sos)))+59
  sos25<-which.min(abs(as.numeric(df_sel_sos$gpp_mod_norm_sm[60:N_length] - trs_25_sos)))+59
  #3)the end of the green-up period
  peak<-pos_max
  #4)additional:also determine eos10 and eos25
  df_sel_eos<-df %>%
    dplyr::filter(doy>=pos_max)
  mn_eos<-min(df_sel_eos$gpp_mod_norm_sm,na.rm = T)
  amp_eos<-mm-mn_eos
  #
  trs_10_eos<-mn_eos+0.1*amp_eos
  trs_25_eos<-mn_eos+0.25*amp_eos
  #
  eos10<-which.min(abs(as.numeric(df_sel_eos$gpp_mod_norm_sm - trs_10_eos)))
  eos25<-which.min(abs(as.numeric(df_sel_eos$gpp_mod_norm_sm - trs_25_eos)))
  #adjust the eos 10 and eos25
  eos25<-eos25+peak-1
  eos10<-eos10+peak-1
  #
  pos_agg<-data.frame(sos10=sos10,sos25=sos25,peak=peak,eos25=eos25,eos10=eos10)
  #plotting demonstration:
  # plot(df_subset$date,df$gpp_mod_norm,ylab = "Norm GPP",xlab = "",pch=".")
  # lines(df_subset$date,df$gpp_mod_roll20,ylab = "Norm GPP",xlab = "",lwd=1.2,col="red")
  # abline(h=c(trs_10_sos,trs_25_sos,trs_25_eos,trs_10_eos),col="red",lty=2)
  # abline(v=df_subset$date[c(sos10,sos25,eos25,eos10)],col="red")
  # text(df_subset$date[c(sos10,sos25,eos25,eos10)],0.8,labels = c("sos10","sos25","eos25","eos10"),col = "red",pos = c(2,4,2,4))
  # abline(v=df_subset$date[c(peak)],col="red")
  # text(df_subset$date[c(peak)],0.8,labels = c("peak"),col = "red",pos = c(4))
  return(pos_agg)
}

#------------------------------------------------
#(3)start to extract the phenos 
#------------------------------------------------
 degf<-0.05  #degree of freedom for the spline smooth
 Clim<-"Cfa"
 PFT<-"DFB"
 site_name<-"US-MMS" #-->using US-MMS as the referred time series 
 avai_years<-c(2000:2014)

##remember add dplyr:: for instance as sometimes two same name functions might conflict with each other
df<-ddf %>%
  dplyr::filter(sitename==site_name) %>%
  dplyr::mutate(doy = yday(date), Year=year(date)) %>%
  dplyr::filter(Year %in% avai_years)
#
Years_unique<-unique(df$Year)

#-------------------------------------------------------------------------------
#Step1: normlization for all the years in one site
#normalized the gpp_obs and gpp_mod using the gpp_max(95 percentile of gpp)
#-------------------------------------------------------------------------------
gpp_max_obs<-quantile(df$gpp_obs,0.95,na.rm = T)
gpp_max_mod<-quantile(df$gpp_mod_FULL,0.95,na.rm = T)
df$gpp_mod_norm<-df$gpp_mod_FULL/gpp_max_mod
df$gpp_obs_norm<-df$gpp_obs/gpp_max_obs
#first to gap-fill the ts to enable the ts smoothing
#then smooth the time series of all the years uisng the spline(did not use at the moment)
#hints:do not gap-fill gpp_obs if there are long gaps in the time sereis-->for biases calculation in step 4-->need to gapfill each year
if(length(df$gpp_obs[!is.na(df$gpp_obs)])>30){
  df_new<-c()
  for(i in 1:length(Years_unique)){
    df_t<-df[df$Year==Years_unique[i],]
    df_t<-df_t %>%
      mutate(gpp_mod_norm_f=na.fill(gpp_mod_norm,c("extend")),gpp_obs_norm_f= na.fill(gpp_obs_norm,c(NA,"extend",NA))) %>%
      mutate(gpp_mod_norm_sm = extract_smooth_y(gpp_mod_norm_f,degf))
    df_new<-rbind(df_new,df_t)
  }
  #coercian gpp_obs to numeric as sometimes long NA makes the rbind wrong to convert to numeric format
  df_new$gpp_obs_norm_f<-as.numeric(df_new$gpp_obs_norm_f)
}
#
N<-length(df$gpp_mod_FULL)
#-------------------------------------------------------------------------------
#Step 2:Determine the green-up period for each year(using spline smoothed values): 
#followed analysis is based on the normlized "GPP_mod"time series(determine earlier sos)
#-------------------------------------------------------------------------------
#
library(zoo)
pos_diffYears<-c()
#also subset data beyond the green-up period
df_outgreenup<-c()
df_final<-c()
for(i in 1:length(Years_unique)){
  df_subset<-df_new[df_new$Year==Years_unique[i],]
  #spot t
  #calculate the rollmean 
  df_subset<-df_subset%>%
    dplyr::mutate(
      gpp_mod_roll5=rollapply(df_subset$gpp_mod_norm_f,5,mean,fill=NA,align="center"),
      gpp_mod_roll7=rollapply(df_subset$gpp_mod_norm_f,7,mean,fill=NA,align="center"),
      gpp_mod_roll10=rollapply(df_subset$gpp_mod_norm_f,10,mean,fill=NA,align="center"),
      gpp_mod_roll15=rollapply(df_subset$gpp_mod_norm_f,15,mean,fill=NA,align="center"),
      gpp_mod_roll20=rollapply(df_subset$gpp_mod_norm_f,20,mean,fill=NA,align="center"),
      gpp_obs_roll5=rollapply(df_subset$gpp_obs_norm_f,5,mean,fill=c(NA,"extend",NA),align="center"),
      gpp_obs_roll7=rollapply(df_subset$gpp_obs_norm_f,7,mean,fill=c(NA,"extend",NA),align="center"),
      gpp_obs_roll10=rollapply(df_subset$gpp_obs_norm_f,10,mean,fill=c(NA,"extend",NA),align="center"),
      gpp_obs_roll15=rollapply(df_subset$gpp_obs_norm_f,15,mean,fill=c(NA,"extend",NA),align="center"),
      gpp_obs_roll20=rollapply(df_subset$gpp_obs_norm_f,20,mean,fill=c(NA,"extend",NA),align="center"),
      temp_obs_roll20=rollapply(df_subset$temp_day_fluxnet2015,20,mean,fill=c(NA,"extend",NA),align="center")
    )
  # pos_sim_max<-match(max(df_subset$gpp_mod_roll20,na.rm = T),df_subset$gpp_mod_roll20)
  pos_sim_max<-match(max(df_subset$gpp_mod_norm_sm),df_subset$gpp_mod_norm_sm)
  pos_max<-pos_sim_max
  #finding out the start and end of the green-up period:
  pos_agg<-Deter_greenup(df_subset,pos_max)
  pos_diffYears<-rbind(pos_diffYears,pos_agg)
  #separate the data beyond greenup period:
  # pos_beyond<-c(1:c(pos_agg$sos10-1),c(c(pos_agg$peak+1):length(df_subset$date)))
  #pos_beyond<-c(pos_agg$peak+1):length(df_subset$date)
  pos_beyond<-c(1:c(pos_agg$sos10-1),c(c(pos_agg$eos10+1):length(df_subset$date)))
  t_outgreenup<-df_subset[pos_beyond,]
  df_outgreenup<-rbind(df_outgreenup,t_outgreenup)
  df_final<-rbind(df_final,df_subset)
}
rownames(pos_diffYears)<-Years_unique

#------------------------------------------------
#(4)making the extraction plots 
#------------------------------------------------
#using the data in 2000 as the reference:refer previous results in https://rpubs.com/yluo_BGC/843413
df_subset<-df_final[df_final$Year==2000,]
#manually to define the doy to fulfill the demonstration purpose
rect.coord_GS=data.frame(x1=df_subset$date[115],
                         x2=df_subset$date[198], 
                         y1=min(df_subset$gpp_mod_norm,na.rm = T)-0.05, 
                         y2=max(df_subset$gpp_mod_norm,na.rm = T)+0.05)

gg_plot<-ggplot() + 
  # geom_point(data=df_subset,aes(x=date,y = gpp_mod_norm),col=adjustcolor("tomato",0.5),pch=1)+
  # geom_line(data=df_subset,aes(x=date,y = gpp_mod_roll20,col="modelled"),lwd=1.2)+
  #using the observation to make the plot
  geom_point(data=df_subset,aes(x=date,y = gpp_obs_norm),col=adjustcolor("black",0.5),pch=1)+
  geom_line(data=df_subset,aes(x=date,y = gpp_obs_roll20),lwd=1.2)+
  #manually to define the doy to fulfill the demonstration purpose
  # geom_vline(xintercept = df_subset$date[115],col="forestgreen",lty=2,lwd=1.2)+
  geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),
               col="forestgreen",lty=2,lwd=1.1,
               data = data.frame(x1=df_subset$date[115],y1= -0.05,
                                 x2=df_subset$date[115],y2=0.18))+
  # geom_vline(xintercept = df_subset$date[198],col="forestgreen",lty=2,lwd=1.2)+
  geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),
               col="forestgreen",lty=2,lwd=1.1,
               data = data.frame(x1=df_subset$date[198],y1= -0.05,
                                 x2=df_subset$date[198],y2=1.05))+
  # geom_vline(xintercept = df_subset$date[300],col="forestgreen",lty=2,lwd=1.2)+
  geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),
               col="forestgreen",lty=2,lwd=1.1,
               data = data.frame(x1=df_subset$date[300],y1= -0.05,
                                 x2=df_subset$date[300],y2=0.18))+
  labs(title = "", x = "DoY", y = "Norm GPP" )+  #(g C m-2 d-1 )
  geom_rect(data=rect.coord_GS,mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="green3",alpha=0.1)+
  #
  geom_hline(yintercept = c(0,1.1),col=adjustcolor("grey",0.5),size=1.05,lty=2)+
  ##adding the amplitude-->working here-->tomorrow working!!!
  
  
  #adding text:
  annotate(geom = "text",x=df_subset$date[160],y=0,label="PRP",size=6)+
  annotate(geom = "text",x=c(df_subset$date[95],df_subset$date[215],df_subset$date[318]),
           y=rep(0.2,3),label=c("SOS","POS","EOS"),size=4)+
  ylim(-0.05,1.16)+
  theme_classic()+
  theme(axis.ticks.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.title = element_text(size=20),
    legend.position = c(0.15,0.9),
    legend.background = element_blank())
#
save.path<-"./manuscript/test_files/phenophase_definition/"
ggsave(paste0(save.path,"gpp_phenophase_define.png"),gg_plot)
