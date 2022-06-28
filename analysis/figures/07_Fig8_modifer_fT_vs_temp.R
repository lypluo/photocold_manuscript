##---------------------------------------
#Aim: the check the variation of modifier fT along with the temperature gradient
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
load(paste0("./data/model_parameters/parameters_MAE_newfT/","optim_par_run5000_eachsite.rds"))
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
# need to remove the sites that do not used in this analysis:
rm.sites<-c("BE-Bra","CA-SF1","CA-SF2","FI-Sod","US-Wi4")
df_final_new<-df_final_new %>%
  filter(sitename!=rm.sites[1] & sitename!=rm.sites[2]&sitename!=rm.sites[3]&sitename!=rm.sites[4]&sitename!=rm.sites[5])

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

#-----------------------------------
#(3) start to summarize the data:
#-----------------------------------
#-----------
#summarize the meteos/fT for each site
#-----------
df_sum_yearly_1<-df_merge_new %>%
  group_by(sitename,year) %>%
  dplyr::summarise(temp=mean(temp),
            prec=sum(prec_fluxnet2015,na.rm = T), #do not use prec from Koen as the data seems strange
            vpd=mean(vpd),
            ppdf=mean(ppfd),
            elv=mean(elv),
            tmin=mean(tmin),
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
#summary site-years for site
#---
df_sum_1<-df_sum_yearly %>%
  group_by(sitename) %>%
  summarise_at(vars(temp:fapar_spl),mean,na.rm=T)
df_sum_2<-df_sum_yearly %>%
  group_by(sitename) %>%
  dplyr::summarise(lon=unique(lon),
            lat=unique(lat),
            classid=unique(classid),
            koeppen_code=unique(koeppen_code))
df_sum<-left_join(df_sum_1,df_sum_2)

##-----------------------
#(4) compare fT difference in differnt group
##----------------------
df_sum$Clim.PFTs<-paste0(df_sum$koeppen_code,"-",df_sum$classid)
#
#first merge the parameters with meteos:
pars_final<-as.data.frame(t(pars_final))
pars_final$sitename<-rownames(pars_final)
#merge:
df_sum_new<-left_join(df_sum,pars_final,by="sitename")

#--------------------------
#(5)plotting:fT vs Tmean
#--------------------------
library(ggrepel)
plot_data<-df_sum_new %>%
  select(sitename,temp,prec,tmin,fT,classid,Clim.PFTs,tau:k)%>%
  mutate(PFT=classid,classid=NULL)
plot_fT<-ggplot()+
  geom_point(data = plot_data,aes(x=temp,y=fT,col=PFT,size=prec))+
  geom_text_repel(data = plot_data,aes(x=temp,y=fT,col=PFT,label=sitename),size=4)+
  scale_color_manual(values = c("DBF"="orange","MF"="cyan","ENF"="magenta"))+
  geom_smooth(data=plot_data,aes(x=temp,y=fT),col="blue",
              method = "lm",formula = y ~ x,lty=2,fill=adjustcolor("steelblue2",0.2))+
  xlab(expression(T[mean]*" (°C)"))+
  ylab(expression(f[T]))+
  theme(
    legend.text = element_text(size=22),
    legend.key.size = unit(2, 'lines'),
    axis.title = element_text(size=26),
    axis.text = element_text(size = 22),
    text = element_text(size=24),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour ="grey",fill="white")
  )

#save the plot
save.path<-"./manuscript/figures/"
ggsave(paste0(save.path,"Figure8_fT_vs_Ta.png"),plot_fT)

