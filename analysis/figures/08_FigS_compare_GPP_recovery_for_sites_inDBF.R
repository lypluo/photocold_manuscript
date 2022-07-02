##---------------------------------------
#Aim: To compare the GPP annual variation pattern for the sites from temperature gradient 
##---------------------------------------
library(dplyr)
library(devtools)
# devtools::load_all("D:/Github/rbeni/")
# install_github("stineb/rbeni")
library(rbeni) #-->make the evaluation plot
library(tidyverse)
library(cowplot)
library(grid)
library(ggpubr)
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

#load the PFTs information:
#load the modis data-->tidy from Beni
#read_rds from tidyverse
load(paste0("./data-raw/raw_data/sites_info/","Pre_selected_sites_info.RDA"))
sites.info<-df_sites_sel

#----
#merge the data
#-----
df_merge<-left_join(df_recent,sites.info,by="sitename")
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

#-----------
#summarize the meteos for each site
#-----------
df_sum_yearly_1<-df_merge_new %>%
  group_by(sitename,year) %>%
  dplyr::summarise(temp=mean(temp),
            prec=sum(prec),
            vpd=mean(vpd),
            ppdf=mean(ppfd),
            elv=mean(elv),
            tmin=mean(tmin),
            tmax=mean(tmax),
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
#(3) compare the parameter difference in differnt group
##----------------------
df_sum$Clim.PFTs<-paste0(df_sum$koeppen_code,"-",df_sum$classid)
#only target the sites we used for the analysis:
##---------------------
#A.load the event_length data-->the sites we were used
#---------------------
load.path<-"./data/event_length/"
load(paste0(load.path,"df_events_length.RDA"))
#
used_sites<-unique(df_events_all$sitename)

#-----select the data for thoses used sites----
#delete the sites do not used:
df_final<-df_sum %>%
  filter(sitename %in% used_sites)

#--------------------------
#(4)plotting:
#--------------------------
#first merge the parameters with meteos:
pars_final<-as.data.frame(t(pars_final))
pars_final$sitename<-rownames(pars_final)
#merge:
df_final_new<-left_join(df_final,pars_final,by="sitename")

#-------------
#plot the GPP annual cycle for DBF sites
#-------------
library(lubridate)
library(colorspace)
GPP_merge<-df_merge_new %>%
  dplyr::filter(classid=="DBF")%>%
  mutate(doy=yday(date))%>%
  dplyr::select(sitename,doy,gpp_obs)%>%
  group_by(sitename)%>%
  mutate(gpp_obs_95=quantile(gpp_obs,probs=0.95,na.rm=T))%>%
  mutate(gpp_norm=gpp_obs/gpp_obs_95)%>%
  dplyr::select(sitename,doy,gpp_norm)%>%
  group_by(sitename,doy)%>%
  dplyr::summarise(obs=mean(gpp_norm,na.rm=T))
#merge GPP_merge with general Climates condition in different sites:
df_plot_final<-left_join(GPP_merge,df_final_new)

#-------
#plotting
#-------
df_plot<-df_plot_final%>%
   ggplot(aes(doy,obs,col=temp))+
   geom_point()+
   geom_smooth(aes(group=temp),method = "loess",se=F,n=30)+
   scale_color_continuous_sequential(palette = "OrYel",rev = T)+
   xlab("DoY")+
   ylab(expression("GPP"[EC]*" (g C m"^-2, " d"^-1, ")" ) )+
   guides(color=guide_legend(title=expression(T[mean])))+##change the legend title
  theme(
    legend.text = element_text(size=20),
    legend.key.size = unit(2, 'lines'),
    axis.title = element_text(size=24),
    axis.text = element_text(size = 20),
    text = element_text(size=24),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour ="grey",fill="white")
  )
print(df_plot)
#save the plot
save.path<-"./manuscript/figures/"
ggsave(paste0(save.path,"FigureS_GPP_cycle_for_diffT.png"),df_plot)

