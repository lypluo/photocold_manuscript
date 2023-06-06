##---------------------------------------
#Aim: the check the variation of modifier fT along with the temperature gradient-->
#the relationship between fT and Tmin in winter or spring
##---------------------------------------
library(dplyr)
devtools::load_all("D:/Github/rbeni/")
library(rbeni) #-->make the evaluation plot
library(tidyverse)
library(cowplot)
library(grid)
library(ggpubr)
library(lubridate)
#---------------------------
#(1)load the calibrated parameters for each site
#---------------------------
load(paste0("./data/model_parameters/parameters_MAE_newfT/",
            "optim_par_run5000_eachsite_new.rds"))
#merge the parameters:
merge_pars<-c()
sites<-names(par_mutisites)
for(i in 1:length(par_mutisites)){
  temp<-t(as.data.frame(par_mutisites[i]))
  merge_pars<-rbind(merge_pars,temp)
}
pars_final<-as.data.frame(t(merge_pars))
names(pars_final)<-sites
#change the parameters name:
rownames(pars_final)<-c("tau","X0","Smax","k")

#----------------------------
#(2)load original data(meteos, gpp..)
#----------------------------
#load the data uploaded by Koen
df_recent <- readRDS(paste0("./data-raw/raw_data/P_model_output/model_data.rds")) %>%
  mutate(
    year = format(date, "%Y")
  ) %>%
  na.omit()

#------------------------------------------
#normalized the GPP-->for each site,
#normalized the "gpp" and "gpp_mod" through their 90 percentiles
#------------------------------------------
#I should use the same value to normlize the gpp and gpp_mod:
gpp_P95<-df_recent %>%
  group_by(sitename) %>%
  dplyr::summarise(gpp_norm_p95=quantile(c(gpp,gpp_mod),0.95,na.rm=T))
#
df_recent<-left_join(df_recent,gpp_P95,by="sitename")
df_recent<-df_recent %>%
  mutate(gpp=gpp/gpp_norm_p95,gpp_mod=gpp_mod/gpp_norm_p95)
##
# need to remove the sites that do not used in this analysis:
rm.sites<-c("BE-Bra","CA-SF1","CA-SF2","FI-Sod","US-Wi4")
df_recent<-df_recent %>%
  filter(sitename!=rm.sites[1] & sitename!=rm.sites[2]&sitename!=rm.sites[3]&sitename!=rm.sites[4]&sitename!=rm.sites[5])
#
sel_sites<-unique(df_recent$sitename)

#-------------------------
#to obtain the modifer fT 
#-------------------------
source(paste0("./R/functions_in_model/newly_formulated_fun/model_fT_rev.R"))
#a.get the stress factor(calibration factor) for each site
df_final<-c()
df_recent$doy<-yday(df_recent$date)
for (i in 1:length(sel_sites)) {
  df_sel<-df_recent %>%
    dplyr::filter(sitename==sel_sites[i])
  
  scaling_factors <- df_sel %>%
    # group_by(sitename, year) %>%
    do({
      scaling_factor <- f_Ts_rev(.,par_mutisites[[i]])
      data.frame(
        sitename = .$sitename,
        date = .$date,
        scaling_factor_optim = scaling_factor
      )
    })
  df_sel <- left_join(df_sel, scaling_factors)
  
  #merge different sites:
  df_final<-rbind(df_final,df_sel)
}
#-----------------------------
#need to back-convert the normalized gpp to gpp
#-----------------------------
df_final_new<-df_final %>%
  mutate(gpp=gpp*gpp_norm_p95,
         gpp_mod=gpp_mod*gpp_norm_p95)
# # need to remove the sites that do not used in this analysis:
# rm.sites<-c("BE-Bra","CA-SF1","CA-SF2","FI-Sod","US-Wi4")
# df_final_new<-df_final_new %>%
#   filter(sitename!=rm.sites[1] & sitename!=rm.sites[2]&sitename!=rm.sites[3]&sitename!=rm.sites[4]&sitename!=rm.sites[5])

#-------------------------------
#load the PFTs information:
#load the modis data-->tidy from Beni
#read_rds from tidyverse
#-------------------------------
load(paste0("./data-raw/raw_data/sites_info/","Pre_selected_sites_info.RDA"))
sites.info<-df_sites_sel
#----
#merge the data
#-----
df_merge<-left_join(df_final_new,sites.info,by="sitename")
df_merge$year<-as.numeric(df_merge$year)
#load the data Beni sent me before:
df_old<-read.csv(file=paste0("./data-raw/raw_data/Data_sent_by_Beni/","ddf_fluxnet2015_pmodel_with_forcings_stocker19gmd.csv"))
df_old<-df_old %>%
  mutate(date=lubridate::mdy(date),
         year=lubridate::year(date)) %>%
  na.omit(gpp_obs)
#----
#merge data:
#----
df_merge_new<-left_join(df_merge,df_old,by=c("sitename", "date", "year"))
#update in Nov,2022-->also add the green-up:
phenology.path<-"./data/event_length/"
load(paste0(phenology.path,"df_events_length.RDA"))
df.pheno<-df_events_all%>%
  select(sitename,Year,sos,peak)%>%
  mutate(year=Year,Year=NULL)
#
df_merge_new<-left_join(df_merge_new,df.pheno)
#-----------------------------------
#(3) start to summarize the data:
#-----------------------------------
df_merge_new<-df_merge_new%>%
  mutate(greenup=ifelse(doy>=sos & doy<=peak,"greenup","Notgreenup"))
#-----------
#A.summarize for the growing season:
#----------
df_sum_yearly_1<-df_merge_new %>%
  group_by(sitename,year,greenup) %>%
  dplyr::summarise(temp=mean(temp),
                   prec=sum(prec_fluxnet2015,na.rm = T), #do not use prec from Koen as the data seems strange
                   vpd=mean(vpd),
                   ppdf=mean(ppfd),
                   elv=mean(elv),
                   tmin_mean=mean(tmin,na.rm = T),
                   tmin_min=min(tmin,na.rm = T),
                   tmax=mean(tmax),
                   fT=mean(scaling_factor_optim),  ##modifier fT
                   fapar_itpl=mean(fapar_itpl),
                   fapar_spl=mean(fapar_spl)
  )
df_sum_yearly_2<-df_merge_new %>%
  group_by(sitename,year) %>%
  dplyr::summarise(lon=unique(lon),
                   lat=unique(lat),
                   classid=unique(classid),
                   koeppen_code=unique(koeppen_code))
df_sum_yearly<-left_join(df_sum_yearly_1,df_sum_yearly_2)
#---
#summary site-years(different seasons) for site
#---
df_sum_1<-df_sum_yearly %>%
  group_by(sitename,greenup) %>%
  summarise_at(vars(temp:fapar_spl),mean,na.rm=T)
df_sum_2<-df_sum_yearly %>%
  group_by(sitename) %>%
  dplyr::summarise(lon=unique(lon),
                   lat=unique(lat),
                   classid=unique(classid),
                   koeppen_code=unique(koeppen_code))
df_sum_greenup<-left_join(df_sum_1,df_sum_2)

#-----------
#B.summarize the meteos/fT for each site
#-----------
df_sum_yearly_1<-df_merge_new %>%
  #update in Nov,2022-->separate the temperature in different season
  mutate(month=month(date),
         season=case_when(month>=3 & month<=5 ~"spring",
                          month>5 & month<=8 ~"summer",
                          month>8 & month<=11 ~"autumn",
                          month>11 | month<3 ~"winter"))%>%
  group_by(sitename,year,season) %>%
  dplyr::summarise(temp=mean(temp),
            prec=sum(prec_fluxnet2015,na.rm = T), #do not use prec from Koen as the data seems strange
            vpd=mean(vpd),
            ppdf=mean(ppfd),
            elv=mean(elv),
            tmin_mean=mean(tmin,na.rm = T),
            tmin_min=min(tmin,na.rm = T),
            tmax=mean(tmax),
            fT=mean(scaling_factor_optim),  ##modifier fT
            fapar_itpl=mean(fapar_itpl),
            fapar_spl=mean(fapar_spl)
            )
df_sum_yearly_2<-df_merge_new %>%
  group_by(sitename,year) %>%
  dplyr::summarise(lon=unique(lon),
            lat=unique(lat),
            classid=unique(classid),
            koeppen_code=unique(koeppen_code))
df_sum_yearly<-left_join(df_sum_yearly_1,df_sum_yearly_2)

#---
#summary site-years(different seasons) for site
#---
df_sum_1<-df_sum_yearly %>%
  group_by(sitename,season) %>%
  summarise_at(vars(temp:fapar_spl),mean,na.rm=T)
df_sum_2<-df_sum_yearly %>%
  group_by(sitename) %>%
  dplyr::summarise(lon=unique(lon),
            lat=unique(lat),
            classid=unique(classid),
            koeppen_code=unique(koeppen_code))
df_sum_season<-left_join(df_sum_1,df_sum_2)

##-----------------------
#(4) compare fT difference in differnt group
##----------------------
df_sum_greenup$Clim.PFTs<-paste0(df_sum_greenup$koeppen_code,"-",df_sum_greenup$classid)
df_sum_season$Clim.PFTs<-paste0(df_sum_season$koeppen_code,"-",df_sum_season$classid)
#
#first merge the parameters with meteos:
pars_final<-as.data.frame(t(pars_final))
pars_final$sitename<-rownames(pars_final)
#merge:
df_sum_greenup_new<-left_join(df_sum_greenup,pars_final,by="sitename")
df_sum_season_new<-left_join(df_sum_season,pars_final,by="sitename")

#--------------------------
#(5)plotting:fT vs Tmin
#--------------------------
#-------
#A.for greenup
#------
library(ggrepel)
library(ggpmisc) #r package for function stat_poly_eq
plot_data<-df_sum_greenup_new %>%
  dplyr::select(sitename,greenup,temp,prec,tmin_mean,tmin_min,fT,classid,Clim.PFTs,tau:k)%>%
  mutate(PFT=classid,classid=NULL)
plot_fT_greenup<-ggplot()+
  geom_point(data = plot_data[plot_data$greenup=="greenup",],aes(x=tmin_min,y=fT,col=PFT,size=prec))+
  geom_text_repel(data = plot_data[plot_data$greenup=="greenup",],aes(x=tmin_min,y=fT,col=PFT,label=sitename),size=4)+
  scale_color_manual(values = c("DBF"="orange","MF"="cyan","ENF"="magenta"))+
  geom_smooth(data=plot_data[plot_data$greenup=="greenup",],aes(x=tmin_min,y=fT),col="blue",
              method = "lm",formula = y ~ x,lty=2,fill=adjustcolor("steelblue2",0.2))+
  stat_poly_eq(data=plot_data[plot_data$greenup=="greenup",],
               aes(x=tmin_min,y=fT,
                   label = paste(
                     after_stat(rr.label),
                     after_stat(p.value.label),
                     sep = "*\", \"*"),
               ),col="blue")+
  xlab(expression(T[min]*" (°C)"))+
  ylab(expression(f[T]))+
  theme(
    legend.text = element_text(size=20),
    legend.key.size = unit(2, 'lines'),
    axis.title = element_text(size=20),
    axis.text = element_text(size = 18),
    text = element_text(size=18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour ="grey",fill="white")
  )

plot_fT_notgreenup<-ggplot()+
  geom_point(data = plot_data[plot_data$greenup=="Notgreenup",],aes(x=tmin_min,y=fT,col=PFT,size=prec))+
  geom_text_repel(data = plot_data[plot_data$greenup=="Notgreenup",],aes(x=tmin_min,y=fT,col=PFT,label=sitename),size=4)+
  scale_color_manual(values = c("DBF"="orange","MF"="cyan","ENF"="magenta"))+
  geom_smooth(data=plot_data[plot_data$greenup=="Notgreenup",],aes(x=tmin_min,y=fT),col="blue",
              method = "lm",formula = y ~ x,lty=2,fill=adjustcolor("steelblue2",0.2))+
  stat_poly_eq(data=plot_data[plot_data$greenup=="Notgreenup",],
               aes(x=tmin_min,y=fT,
                   label = paste(
                     after_stat(rr.label),
                     after_stat(p.value.label),
                     sep = "*\", \"*"),
               ),col="blue")+
  xlab(expression(T[min]*" (°C)"))+
  ylab(expression(f[T]))+
  theme(
    legend.text = element_text(size=20),
    legend.key.size = unit(2, 'lines'),
    axis.title = element_text(size=20),
    axis.text = element_text(size = 18),
    text = element_text(size=18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour ="grey",fill="white")
  )
#-------
#B.for different season:
#------
plot_data<-df_sum_season_new %>%
  dplyr::select(sitename,season,temp,prec,tmin_mean,tmin_min,fT,classid,Clim.PFTs,tau:k)%>%
  mutate(PFT=classid,classid=NULL)
####
plot_fT_spring<-ggplot()+
  geom_point(data = plot_data[plot_data$season=="spring",],aes(x=tmin_min,y=fT,col=PFT),size=4)+
  geom_text_repel(data = plot_data[plot_data$season=="spring",],aes(x=tmin_min,y=fT,col=PFT,label=sitename),size=4)+
  scale_color_manual(values = c("DBF"="orange","MF"="cyan","ENF"="magenta"))+
  geom_smooth(data=plot_data[plot_data$season=="spring",],aes(x=tmin_min,y=fT),col="blue",
              method = "lm",formula = y ~ x,lty=2,fill=adjustcolor("steelblue2",0.2))+
  stat_poly_eq(data=plot_data[plot_data$season=="spring",],
                 aes(x=tmin_min,y=fT,
                     label = paste(
                                   after_stat(rr.label),
                                   after_stat(p.value.label),
                                   sep = "*\", \"*"),
                     ),col="blue")+
  xlab(expression("spring  "*T[min]*" (°C)"))+
  ylab(expression("spring  "*f[T]))+
  # geom_hline(yintercept = 1,lty=2)+
  theme(
    legend.position = "none",
    legend.text = element_text(size=20),
    legend.key.size = unit(2, 'lines'),
    axis.title = element_text(size=20),
    axis.text = element_text(size = 18),
    text = element_text(size=18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour ="grey",fill="white")
  )
####
plot_fT_winter<-ggplot()+
  geom_point(data = plot_data[plot_data$season=="winter",],aes(x=tmin_min,y=fT,col=PFT),size=4)+
  geom_text_repel(data = plot_data[plot_data$season=="winter",],aes(x=tmin_min,y=fT,col=PFT,label=sitename),size=4)+
  scale_color_manual(values = c("DBF"="orange","MF"="cyan","ENF"="magenta"))+
  geom_smooth(data=plot_data[plot_data$season=="winter",],aes(x=tmin_min,y=fT),col="blue",
              method = "lm",formula = y ~ x,lty=2,fill=adjustcolor("steelblue2",0.2))+
  stat_poly_eq(data=plot_data[plot_data$season=="winter",],
               aes(x=tmin_min,y=fT,
                   label = paste(
                     after_stat(rr.label),
                     after_stat(p.value.label),
                     sep = "*\", \"*"),
               ),col="blue")+
  xlab(expression("winter  "*T[min]*" (°C)"))+
  ylab(expression("winter  "*f[T]))+
  # geom_hline(yintercept = 1,lty=2)+
  theme(
    legend.position = "none",
    legend.text = element_text(size=20),
    legend.key.size = unit(2, 'lines'),
    axis.title = element_text(size=20),
    axis.text = element_text(size = 18),
    text = element_text(size=18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour ="grey",fill="white")
  )

###Additonal plot:Tmin_winter vs fT_spring:
plot_data_winter<-plot_data%>%
  filter(season=="winter")%>%
  select(sitename,temp,prec,tmin_mean,tmin_min) #select its temperature
plot_data_spring<-plot_data%>%
  filter(season=="spring")%>%
  select(sitename,fT,Clim.PFTs,PFT)
plot_data_new_test<-left_join(plot_data_winter,plot_data_spring)

plot_winterT_springfT<-ggplot()+
  geom_point(data = plot_data_new_test,aes(x=tmin_min,y=fT,col=PFT),size=4)+
  geom_text_repel(data = plot_data_new_test,aes(x=tmin_min,y=fT,col=PFT,label=sitename),size=4)+
  # scale_color_manual(values = c("DBF"="orange","MF"="cyan","ENF"="magenta"))+
  khroma::scale_color_highcontrast(aesthetics = "color")+
  geom_smooth(data=plot_data_new_test,aes(x=tmin_min,y=fT),col="blue",
              method = "lm",formula = y ~ x,lty=2,fill=adjustcolor("steelblue2",0.2))+
  stat_poly_eq(data=plot_data_new_test,
               aes(x=tmin_min,y=fT,
                   label = paste(
                     after_stat(rr.label),
                     after_stat(p.value.label),
                     sep = "*\", \"*"),
               ),col="blue")+
  xlab(expression("winter  "*T[min]*" (°C)"))+
  ylab(expression("Mean spring  "*f[T]))+
  # geom_hline(yintercept = 1,lty=2)+
  theme(
    legend.position = "none",
    legend.text = element_text(size=20),
    legend.key.size = unit(2, 'lines'),
    axis.title = element_text(size=20),
    axis.text = element_text(size = 18),
    text = element_text(size=18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour ="grey",fill="white")
  )

#-------
#C.adding the relationship: winter Tmin vs GPP bias in the early spring
#update in Nov, 2022
#update in Jan, 2023-->both before and after applying the temperature acclimation factor
#------
##C1:GPP bias calculated using original P-model
#load the GPP data(has been characterized as "overestimated" or "non-overestimated")
load.path<-"./data/data_used/"
#from new method:
load(paste0(load.path,"ddf_labeled_norm_trs_newmethod_all_overestimation_Fluxnet2015_sites.RDA"))
df_all_sites<-ddf_labeled;rm(ddf_labeled) 
#select the needed variables 
df_GPP<-df_all_sites%>%
  select(sitename,date,doy,greenup,gpp_obs,gpp_mod_FULL,gpp_res)%>%
  filter(greenup=="yes")%>% ##selecting the green-up period data
  group_by(sitename)%>%
  dplyr::summarise(gpp_obs_GP=mean(gpp_obs,na.rm=T),
      gpp_mod_GP=mean(gpp_mod_FULL,na.rm=T),
      gpp_bias=mean(c(gpp_mod_GP - gpp_obs_GP),na.rm=T),
                   #gpp_bias_rel:relative bias
      gpp_bias_rel=mean(c(gpp_mod_GP - gpp_obs_GP)/gpp_obs_GP,na.rm=T)*100
  )   #GP:for green-up(GPP resumption) period
##C2:GPP bias calculated using P-model with applying the temperature acclimation factor:
df_GPP_optim<-df_merge_new %>%
  select(sitename,date,doy,greenup,gpp_obs,gpp_mod,scaling_factor_optim)%>%
  dplyr::mutate(gpp_mod_optim=gpp_mod*scaling_factor_optim,
         gpp_mod=NULL,scaling_factor_optim=NULL)%>%
  filter(greenup=="greenup")%>% ##selecting the green-up period data
  group_by(sitename)%>%
  dplyr::summarise(gpp_obs_GP=mean(gpp_obs,na.rm=T),
                   gpp_mod_optim_GP=mean(gpp_mod_optim,na.rm=T),
                   gpp_bias=mean(c(gpp_mod_optim_GP - gpp_obs_GP),na.rm=T),
                   #gpp_bias_rel:relative bias
                   gpp_bias_rel_optim=mean(c(gpp_mod_optim_GP - gpp_obs_GP)/gpp_obs_GP,na.rm=T)*100
  )   #GP:for green-up(GPP resumption) period

###(a) for Prior GPP bias:
df_winterTmin_bias_prior<-left_join(plot_data_new_test,df_GPP)
#plotting
plot_winterT_GPPbia_prior<-ggplot()+
  geom_point(data = df_winterTmin_bias_prior,aes(x=tmin_min,y=gpp_bias_rel,col=PFT),size=4)+
  geom_text_repel(data = df_winterTmin_bias_prior,aes(x=tmin_min,y=gpp_bias_rel,col=PFT,label=sitename),size=4)+
  scale_color_manual(values = c("DBF"="orange","MF"="cyan","ENF"="magenta"))+
  geom_smooth(data=df_winterTmin_bias_prior,aes(x=tmin_min,y=gpp_bias_rel),col="blue",
              method = "lm",formula = y ~ x,lty=2,fill=adjustcolor("steelblue2",0.2))+
  khroma::scale_color_highcontrast(aesthetics = "color")+
  stat_poly_eq(data=df_winterTmin_bias_prior,
               aes(x=tmin_min,y=gpp_bias_rel,
                   label = paste(
                     after_stat(rr.label),
                     after_stat(p.value.label),
                     sep = "*\", \"*"),
               ),col="blue")+
  xlab(expression("winter  "*T[min]*" (°C)"))+
  ylab(expression("Rel. bias in GPP"[Pmodel]*" (%)"))+
  # geom_hline(yintercept = 1,lty=2)+
  ylim(-30,150)+
  theme(
    # legend.position = c(0.7,0.75),
    legend.position = "none",
    legend.background = element_blank(),
    legend.direction = "horizontal",
    legend.text = element_text(size=16),
    legend.key.size = unit(2, 'lines'),
    axis.title = element_text(size=20),
    axis.text = element_text(size = 18),
    text = element_text(size=18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour ="grey",fill="white")
  )

###(b) for Post GPP bias:
df_winterTmin_bias_post<-left_join(plot_data_new_test,df_GPP_optim)
#plotting
plot_winterT_GPPbias_post<-ggplot()+
  geom_point(data = df_winterTmin_bias_post,aes(x=tmin_min,y=gpp_bias_rel_optim,col=PFT),size=4)+
  geom_text_repel(data = df_winterTmin_bias_post,aes(x=tmin_min,y=gpp_bias_rel_optim,col=PFT,label=sitename),size=4)+
  scale_color_manual(values = c("DBF"="orange","MF"="cyan","ENF"="magenta"))+
  #update in Jan,2023
  # geom_smooth(data=df_winterTmin_bias_post,aes(x=tmin_min,y=gpp_bias_rel_optim),col="white",se=FALSE,
  #             method = "lm",formula = y ~ x,lty=2,fill=adjustcolor("steelblue2",0.2))+
  khroma::scale_color_highcontrast(aesthetics = "color")+
  stat_poly_eq(data=df_winterTmin_bias_post,
               aes(x=tmin_min,y=gpp_bias_rel_optim,
                   label = paste(
                     after_stat(rr.label),
                     after_stat(p.value.label),
                     sep = "*\", \"*"),
               ),label.y = 0.85,
               col="blue")+
  xlab(expression("winter  "*T[min]*" (°C)"))+
  ylab(expression("Rel. bias in GPP"[adj]*" (%)"))+
  ylim(-30,150)+
  # geom_hline(yintercept = 1,lty=2)+
  theme(
    legend.position = c(0.7,0.75),
    # legend.position = "none",
    legend.background = element_blank(),
    legend.direction = "horizontal",
    legend.text = element_text(size=16),
    legend.key.size = unit(2, 'lines'),
    axis.title = element_text(size=20),
    axis.text = element_text(size = 18),
    text = element_text(size=18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour ="grey",fill="white")
  )

##merge the plots:
library(cowplot)
# plot_fT<-plot_grid(plot_fT_spring,plot_fT_winterT_springfT,nrow = 2,align = 'hv')
#using ggarrange from ggpubr package:
library(ggpubr)
# ggarrange(plot_fT_spring,plot_fT_winter,nrow = 2,common.legend = TRUE,legend = "bottom")
# plot_fT<-ggarrange(plot_fT_spring,
#         plot_winterT_springfT,plot_winterT_GPPbias,align = "v",
#         nrow = 3,common.legend = TRUE,legend = "bottom")
#plot 1-->not use in Jan, 2023
# plot_fT<-plot_grid(plot_fT_spring,
#                    plot_winterT_springfT,plot_winterT_GPPbias,
#                    labels = c("a","b","c"),
#                    align = "v",
#                    nrow = 3,ncol=1)
#
plot_fT<-plot_grid(plot_winterT_springfT,
                   plot_winterT_GPPbia_prior,
                   plot_winterT_GPPbias_post,
                   labels = c("a","b","c"),
                   align = "v",
                   nrow = 3,ncol=1)


#save the plot
save.path<-"./manuscript/figures/"
ggsave(paste0(save.path,"Figure8_fT_vs_Ta_bias_new.png"),plot_fT,height = 13,width = 9)

