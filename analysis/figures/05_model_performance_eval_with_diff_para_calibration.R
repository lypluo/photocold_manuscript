#################################################
#Aim:to evaluate the model performances with the parameters from the model 
#calibrated on the different scale 
#(each site(site-level), PFT scale(PFT-level), pool all the sites together(All sites-level))
#################################################
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(lme4)
library(tidyverse)
# remotes::install_github("computationales/ingestr") #install the package
library(ingestr)
devtools::load_all("D:/Github/rbeni/")
library(rbeni)

#-------------------------
#(1)load the data and hardening funciton
#-------------------------
#####
#load the data uploaded by Koen
df_recent <- readRDS(paste0("./data-raw/raw_data/P_model_output/model_data.rds")) %>%
  mutate(
    year = format(date, "%Y")
  ) %>%
  na.omit()

#load the data Beni sent me before:
df_old<-read.csv(file=paste0("./data-raw/raw_data/Data_sent_by_Beni/","ddf_fluxnet2015_pmodel_with_forcings_stocker19gmd.csv"))
df_old<-df_old %>%
  mutate(date=lubridate::mdy(date),
         year=lubridate::year(date)) %>%
  na.omit(gpp_obs)
#load the temperature acclimation factor fT:
# source(paste0("./R/functions_in_model/model_hardening_byBeni_addbaseGDD_rev.R"))
source(paste0("./R/functions_in_model/newly_formulated_fun/model_fT_rev.R"))

#--------------------------------------------------------------
#(2) retreive the optimized parameter for the sites
#--------------------------------------------------------------
# set initial value
par <- c("tau"=5,"X0"=-10,"Smax"=5,"k"=1)
#
lower=c(1,-10,5,0)
upper=c(25,10,25,2)

# run model and compare to true values
# returns the RMSE
cost <- function(
  data,
  par
) {
  
  scaling_factor <- data %>%
    # group_by(sitename) %>%
    do({
      scaling_factor <- f_Ts_rev(
        .,
        par
      )
      
      data.frame(
        sitename = .$sitename,
        date = .$date,
        scaling_factor = scaling_factor
      )
    })
  
  df <- left_join(data, scaling_factor)
  
  #rmse
  # rmse <- sqrt(
  #   sum(
  #     (df$gpp - df$gpp_mod * df$scaling_factor)^2)
  #   )/nrow(df)
  #mse:mean square error
  mse<-mean((df$gpp - df$gpp_mod * df$scaling_factor)^2,na.rm=T)
  
  return(mse)
}

#--------------------------------------------------------------
#(3) parameters optimization (for site-level, PFT-level, and All-sites level)
#--------------------------------------------------------------
#a.adding the PFTs information for site-->load the modis data-->tidy from Beni
load(paste0("./data-raw/raw_data/sites_info/","Pre_selected_sites_info.RDA"))
sites.info<-df_sites_sel
#
df_merge<-df_recent %>%
  left_join(
    sites.info,
    by = "sitename"
  )
##adding day of the year
df_merge$doy<-yday(df_merge$date)

#main Clim-PFTs
df_merge$Clim_PFTs<-paste0(df_merge$koeppen_code,"-",df_merge$classid)
Clim.PFTs<-sort(unique(df_merge$Clim_PFTs))

#------------------------------------------
#b.normalized the GPP-->for each site,
#normalized the "gpp" and "gpp_mod" through their 90 percentiles
#------------------------------------------
#I should use the same value to normlize the gpp and gpp_mod:
gpp_P95<-df_merge %>%
  group_by(sitename) %>%
  dplyr::summarise(gpp_norm_p95=quantile(c(gpp,gpp_mod),0.95,na.rm=T))
#
df_merge.new<-left_join(df_merge,gpp_P95,by="sitename")
df_merge.new<-df_merge.new %>%
  mutate(gpp=gpp/gpp_norm_p95,gpp_mod=gpp_mod/gpp_norm_p95)

# need to remove the sites that do not used in this analysis:
rm.sites<-c("BE-Bra","CA-SF1","CA-SF2","FI-Sod","US-Wi4")
df_merge_new<-df_merge.new %>%
  filter(sitename!=rm.sites[1] & sitename!=rm.sites[2]&sitename!=rm.sites[3]&sitename!=rm.sites[4]&sitename!=rm.sites[5])

#------------------------
#c.optimization:
#-----------------------
#this step is done in other script, save the parameterized parameters in the folder..
#....

#--------------------------------------------------------------
#(4) compare the gpp_obs, ori modelled gpp, and gpp modelled 
# using optimated parameters calibrated on different scale(for site-level, PFT-level, and All-sites level)
#--------------------------------------------------------------
#---a.load the optimized parameters from different scales
# site-level
load(paste0("./data/model_parameters/parameters_MAE_newfT/","optim_par_run5000_eachsite.rds"))
#PFT-level
# load(paste0("./data/model_parameters/parameters_MAE_newfT/","optim_par_run5000_PFTs.rds"))
#with updated parameter
load(paste0("./data/model_parameters/parameters_MAE_newfT/","optim_par_run5000_PFTs_with_newMF_paras.rds"))
#All-sites level
load(paste0("./data/model_parameters/parameters_MAE_newfT/","optim_par_run5000_allsites.rds"))

#---b.retrieve the stress factor(calibration factor) for each scale-level
#--------
#for optimized paramater from site-level calibration:
#--------
sel_sites<-unique(df_merge_new$sitename)
df_final_sitelevel<-c()
for (i in 1:length(sel_sites)) {
  df_sel<-df_merge_new %>%
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
  df_final_sitelevel<-rbind(df_final_sitelevel,df_sel)
}
#--------
#for optimized paramater from PFT-level calibration:
#--------
PFTs<-unique(df_merge_new$classid)
df_final_PFTlevel<-c()
for (i in 1:length(PFTs)) {
  df_sel<-df_merge_new %>%
    dplyr::filter(classid==PFTs[i])
  
  scaling_factors <- df_sel %>%
    # group_by(sitename, year) %>%
    do({
      scaling_factor <- f_Ts_rev(.,par_PFTs[[i]])
      data.frame(
        sitename = .$sitename,
        date = .$date,
        scaling_factor_optim = scaling_factor
      )
    })
  df_sel_new <- left_join(df_sel, scaling_factors)
  
  #merge different sites:
  df_final_PFTlevel<-rbind(df_final_PFTlevel,df_sel_new)
}
#--------
#for optimized paramater from All sites-level calibration:
#--------
df_final_Allsiteslevel<-c()
#
df_sel<-df_merge_new
scaling_factors <- df_sel %>%
    # group_by(sitename, year) %>%
    do({
      scaling_factor <- f_Ts_rev(.,par_allsites)
      data.frame(
        sitename = .$sitename,
        date = .$date,
        scaling_factor_optim = scaling_factor
      )
    })
df_sel_new <- left_join(df_sel, scaling_factors)
  
#merge different sites:
df_final_Allsiteslevel<-df_sel_new

#-----------------------------
#---c.need to back-convert the normalized gpp to gpp
#-----------------------------
#for site-level:
df_final_sitelevel_new<-df_final_sitelevel %>%
  mutate(gpp=gpp*gpp_norm_p95,
         gpp_mod=gpp_mod*gpp_norm_p95,
         year=year(date))
df_merge_sitelevel<-left_join(df_final_sitelevel_new,df_old,by = c("sitename", "date", "year")) %>%
  mutate(gpp_obs_recent=gpp,
         gpp_obs_old=gpp_obs,
         gpp_mod_FULL_ori=gpp_mod_FULL,
         gpp_mod_recent_ori=gpp_mod,
         gpp_mod_recent_optim=gpp_mod*scaling_factor_optim,
         gpp=NULL,
         gpp_obs=NULL,
         gpp_mod=NULL)

#for PFT-level:
df_final_PFTlevel_new<-df_final_PFTlevel %>%
  mutate(gpp=gpp*gpp_norm_p95,
         gpp_mod=gpp_mod*gpp_norm_p95,
         year=year(date))
df_merge_PFTlevel<-left_join(df_final_PFTlevel_new,df_old,by = c("sitename", "date", "year")) %>%
  mutate(gpp_obs_recent=gpp,
         gpp_obs_old=gpp_obs,
         gpp_mod_FULL_ori=gpp_mod_FULL,
         gpp_mod_recent_ori=gpp_mod,
         gpp_mod_recent_optim=gpp_mod*scaling_factor_optim,
         gpp=NULL,
         gpp_obs=NULL,
         gpp_mod=NULL)

#for Allsites-level:
df_final_Allsiteslevel_new<-df_final_Allsiteslevel %>%
  mutate(gpp=gpp*gpp_norm_p95,
         gpp_mod=gpp_mod*gpp_norm_p95,
         year=year(date))
df_merge_Allsiteslevel<-left_join(df_final_Allsiteslevel_new,df_old,by = c("sitename", "date", "year")) %>%
  mutate(gpp_obs_recent=gpp,
         gpp_obs_old=gpp_obs,
         gpp_mod_FULL_ori=gpp_mod_FULL,
         gpp_mod_recent_ori=gpp_mod,
         gpp_mod_recent_optim=gpp_mod*scaling_factor_optim,
         gpp=NULL,
         gpp_obs=NULL,
         gpp_mod=NULL)
#--------------------------
#(5).making the evaluation plots:
#using modelled and observed gpp to make the plots
#-------------------------
library(cowplot)

##general assessment: how much bias was reduced(during the green-up period) for 
##original p-model and adjusted p-model(using PFT-specific parameters):
#load the green-up data:
pheno.path<-"./data/event_length/"
load(paste0(pheno.path,"df_events_length.RDA"))
df_pheno<-df_events_all %>%
  select(sitename,Year,sos,peak)%>%
  mutate(year=Year,Year=NULL)
df_merge_PFTlevel<-left_join(df_merge_PFTlevel,df_pheno)

## update in Nov,2022
df_modobs_comp<-df_merge_PFTlevel%>%
  filter(doy>=sos & doy<=peak) %>% ##only select the data during the growing season 
  select(sitename,date,year,gpp_obs_recent,gpp_mod_FULL_ori,gpp_mod_recent_ori,gpp_mod_recent_optim)%>%
  mutate(gpp_obs=gpp_obs_recent,
         gpp_mod_old_ori=gpp_mod_FULL_ori,             #gpp_mod_old_ori-->corrsponds to Stocker et al., 2022
         gpp_mod_recent_ori=gpp_mod_recent_ori,        #gpp_mod_recent_ori-->updated GPP from Beni
         gpp_mod_recent_optim=gpp_mod_recent_optim) %>%#gpp_mod_recent_optim-->updated GPP calibrated with paras
  mutate(gpp_obs_recent=NULL,
         gpp_mod_FULL_ori=NULL)
##how much bias is reduced for the model with optimized parameter compaerd to old original p-model:
library(sirad)
MAE_Pmodel_ori<-unlist(modeval(df_modobs_comp$gpp_mod_old_ori,df_modobs_comp$gpp_obs,
                        stat = c("MAE","RMAE","RMSE","RRMSE")))
MAE_Pmodel_optim<-unlist(modeval(df_modobs_comp$gpp_mod_recent_optim,df_modobs_comp$gpp_obs,
                               stat = c("MAE","RMAE","RMSE","RRMSE")))
#reduced bias:
c(MAE_Pmodel_ori[1]-MAE_Pmodel_optim[1])/MAE_Pmodel_optim[1]

###############################
#---a.making scatter plots:evaluation for all sites (GPP_obs vs GPP_adj(with optimized paras))
###############################
#for the site-level comparision(pooled all sites for this scatter plot):
df_modobs_sitelevel<-df_merge_sitelevel%>%
  select(sitename,date,gpp_obs_recent,gpp_mod_FULL_ori,gpp_mod_recent_ori,gpp_mod_recent_optim)%>%
  mutate(gpp_obs=gpp_obs_recent,
       gpp_mod_old_ori=gpp_mod_FULL_ori,             #gpp_mod_old_ori-->corrsponds to Stocker et al., 2022
       gpp_mod_recent_ori=gpp_mod_recent_ori,        #gpp_mod_recent_ori-->updated GPP from Beni
       gpp_mod_recent_optim=gpp_mod_recent_optim) %>%#gpp_mod_recent_optim-->updated GPP calibrated with paras
  mutate(gpp_obs_recent=NULL,
         gpp_mod_FULL_ori=NULL)
#
plot_gpp_modobs_sitelevel<-df_modobs_sitelevel %>%
     analyse_modobs2("gpp_mod_recent_optim", "gpp_obs", type = "heat")

#for the PFT-level comparision(pooled all sites for this scatter plot):
df_modobs_PFTlevel<-df_merge_PFTlevel%>%
  select(sitename,date,gpp_obs_recent,gpp_mod_FULL_ori,gpp_mod_recent_ori,gpp_mod_recent_optim)%>%
  mutate(gpp_obs=gpp_obs_recent,
         gpp_mod_old_ori=gpp_mod_FULL_ori,             #gpp_mod_old_ori-->corrsponds to Stocker et al., 2022
         gpp_mod_recent_ori=gpp_mod_recent_ori,        #gpp_mod_recent_ori-->updated GPP from Beni
         gpp_mod_recent_optim=gpp_mod_recent_optim) %>%#gpp_mod_recent_optim-->updated GPP calibrated with paras
  mutate(gpp_obs_recent=NULL,
         gpp_mod_FULL_ori=NULL)
#
plot_gpp_modobs_PFTlevel<-df_modobs_PFTlevel %>%
  analyse_modobs2("gpp_mod_recent_optim", "gpp_obs", type = "heat")

#for the Allsites-level comparision(pooled all sites for this scatter plot):
df_modobs_Allsiteslevel<-df_merge_Allsiteslevel%>%
  select(sitename,date,gpp_obs_recent,gpp_mod_FULL_ori,gpp_mod_recent_ori,gpp_mod_recent_optim)%>%
  mutate(gpp_obs=gpp_obs_recent,
         gpp_mod_old_ori=gpp_mod_FULL_ori,             #gpp_mod_old_ori-->corrsponds to Stocker et al., 2022
         gpp_mod_recent_ori=gpp_mod_recent_ori,        #gpp_mod_recent_ori-->updated GPP from Beni
         gpp_mod_recent_optim=gpp_mod_recent_optim) %>%#gpp_mod_recent_optim-->updated GPP calibrated with paras
  mutate(gpp_obs_recent=NULL,
         gpp_mod_FULL_ori=NULL)
#
plot_gpp_modobs_Allsiteslevel<-df_modobs_Allsiteslevel %>%
  analyse_modobs2("gpp_mod_recent_optim", "gpp_obs", type = "heat")


##additional-->calculate the stats
library(sirad)
stats_sitelevel_allsitespooled<-round(unlist(modeval(df_merge_sitelevel$gpp_mod_recent_optim,
      df_merge_sitelevel$gpp_obs_recent,stat = c("MAE","RMSE","R2"))),2)
stats_PFTlevel_allsitespooled<-round(unlist(modeval(df_merge_PFTlevel$gpp_mod_recent_optim,
      df_merge_PFTlevel$gpp_obs_recent,stat = c("MAE","RMSE","R2"))),2)
stats_Allsitelevel_allsitespooled<-round(unlist(modeval(df_merge_Allsiteslevel$gpp_mod_recent_optim,
      df_merge_Allsiteslevel$gpp_obs_recent,stat = c("MAE","RMSE","R2"))),2)
#update in Nov,2022:also calculate the stats for the original p-model gpp (before the calibration):
stats_priorCalibration_allsitespooled<-round(unlist(modeval(df_merge_sitelevel$gpp_mod_FULL_ori,
      df_merge_sitelevel$gpp_obs_recent,stat = c("MAE","RMSE","R2"))),2)


#---
#change x,y axis labels
#---
plot.theme<-theme(
  legend.text = element_text(size=20),
  axis.title = element_text(size=24),
  axis.text = element_text(size = 20),
  text = element_text(size=24)
)
plot_gpp_modobs_sitelevel$gg<-plot_gpp_modobs_sitelevel$gg+
  # xlab("")+
  ylab("")+
  xlab(expression("GPP"[adj]*" (g C m"^-2*"d"^-1*")"))+
  ylab(expression("GPP"[obs]*" (g C m"^-2*"d"^-1*")"))+
  xlim(0,25)+ylim(-5,25)+
  annotate(geom="text",x=20,y=0,label="site-specific paras",size=6)+
  plot.theme
plot_gpp_modobs_PFTlevel$gg<-plot_gpp_modobs_PFTlevel$gg+
  # xlab("")+
  xlab(expression("GPP"[adj]*" (g C m"^-2*"d"^-1*")"))+
  ylab(expression("GPP"[obs]*" (g C m"^-2*"d"^-1*")"))+
  xlim(0,25)+ylim(-5,25)+
  annotate(geom="text",x=20,y=0,label="PFT-specific paras",size=6)+
  plot.theme
plot_gpp_modobs_Allsiteslevel$gg<-plot_gpp_modobs_Allsiteslevel$gg+
  xlab("")+
  ylab("")+
  xlab(expression("GPP"[adj]*" (g C m"^-2*"d"^-1*")"))+
  ylab(expression("GPP"[obs]*" (g C m"^-2*"d"^-1*")"))+
  xlim(0,25)+ylim(-5,25)+
  annotate(geom="text",x=20,y=0,label="One general paras",size=6)+
  plot.theme

#merge the plots
evaulation_merge_plot<-plot_grid(plot_gpp_modobs_sitelevel$gg,
       plot_gpp_modobs_PFTlevel$gg,plot_gpp_modobs_Allsiteslevel$gg,
       widths=15,heights=4,
       labels = "auto",ncol =3,nrow = 1,label_size = 12,align = "hv")
evaulation_merge_plot

###############################
#---b.RMSE and R2 for eval performance (display for different scales)
###############################
library(sirad)
#function#
stat_fun<-function(df,paras_level){
  # df<-df_merge_sitelevel
  # paras_level<-"site"
  # 
  df.stat<-df %>%
    select(sitename,date,Clim_PFTs,gpp_mod_recent_optim,gpp_obs_recent)%>%
    dplyr::mutate(gpp_adj=gpp_mod_recent_optim,gpp_mod_recent_optim=NULL,
                  gpp_obs=gpp_obs_recent,gpp_obs_recent=NULL)%>%
    group_by(Clim_PFTs)%>%
    dplyr::summarise(N=as.numeric(unlist(modeval(gpp_adj,gpp_obs,stat = "N"))),
                     Rsquare=as.numeric(unlist(modeval(gpp_adj,gpp_obs,stat = "R2"))),
                     MAE=as.numeric(unlist(modeval(gpp_adj,gpp_obs,stat="MAE"))),
                     RMSE=as.numeric(unlist(modeval(gpp_adj,gpp_obs,stat="RMSE"))))%>%
    mutate(para_flag=paras_level)
  return(df.stat)
}

#for the parameters is calibrated for site-level(RMSE,R2 for each sites)
#evaluation on the PFT-Clim categories
stats_para_f_site<-stat_fun(df_merge_sitelevel,"site-specific")
#for the parameters is calibrated for PFT-level
#evaluation on the PFT-Clim categories
stats_para_f_PFT<-stat_fun(df_merge_PFTlevel,"PFT-specific")
#for the parameters is calibrated for Allsites-level(general)
#evaluation on the PFT-Clim categories
stats_para_f_Allsites<-stat_fun(df_merge_Allsiteslevel,"general")
##update in Nov,2022-->also calculate the stats between original p-model gpp and gpp_obs:
stats_prior_cali<-df_merge_Allsiteslevel%>%
  select(sitename,date,Clim_PFTs,gpp_mod_FULL_ori,gpp_obs_recent)%>%
  dplyr::mutate(gpp_pmodel=gpp_mod_FULL_ori,gpp_mod_FULL_ori=NULL,
                gpp_obs=gpp_obs_recent,gpp_obs_recent=NULL)%>%
  group_by(Clim_PFTs)%>%
  dplyr::summarise(N=as.numeric(unlist(modeval(gpp_pmodel,gpp_obs,stat = "N"))),
                   Rsquare=as.numeric(unlist(modeval(gpp_pmodel,gpp_obs,stat = "R2"))),
                   MAE=as.numeric(unlist(modeval(gpp_pmodel,gpp_obs,stat="MAE"))),
                   RMSE=as.numeric(unlist(modeval(gpp_pmodel,gpp_obs,stat="RMSE"))))

####summary the stats:
stats_all<-rbind(rbind(stats_para_f_site,stats_para_f_PFT),stats_para_f_Allsites)
stats_all_tidy<-stats_all %>%
  select(Clim_PFTs,Rsquare,MAE,RMSE,para_flag)%>%
  pivot_longer(c(Rsquare,MAE,RMSE),names_to = "eval_metrics",values_to = "metrics")
#set the order:
stats_all_tidy$para_flag<-factor(stats_all_tidy$para_flag,
      levels = c("site-specific","PFT-specific","general"))
#sum for all Clim_PFTs groups
stats_all_tidy_sum<-stats_all_tidy %>%
  group_by(para_flag,eval_metrics) %>%
  dplyr::summarise(metrics_mean=mean(metrics),SD=sd(metrics))
#update Nov, 2022-->also for the statsbetween original p-model gpp and gpp_obs:
stats_prior_cali%>%
  select(Clim_PFTs,Rsquare,MAE,RMSE)%>%
  pivot_longer(c(Rsquare,MAE,RMSE),names_to = "eval_metrics",values_to = "metrics")%>%
  group_by(eval_metrics) %>%
  dplyr::summarise(metrics_mean=mean(metrics),SD=sd(metrics))

###making the plots:
df_plot<-stats_all_tidy %>%
  ggplot(aes(x=eval_metrics,y=metrics))+
  geom_point(size=3,col=adjustcolor("gray"))+
  geom_boxplot(size=1.1,col=adjustcolor("steelblue",0.9),fill=adjustcolor("white",0.1))+
  facet_wrap(.~para_flag)+
  # facet_wrap(para_flag~.)+
  xlab("")+ylab("")+
  theme_classic()+
  plot.theme
#add the point_range
# df_stats_final<-df_plot+
#   geom_pointrange(data = stats_all_tidy_sum,
#      aes(x=eval_metrics,y=metrics_mean,
#          ymin=metrics_mean-SD,ymax=metrics_mean+SD),
#      col=adjustcolor("steelblue",0.9),size=1.1)
#change the panel labels:
df_plot<-plot_grid(df_plot,nrow=1,labels = "d")

###put plots a and b together 
#merge the plots
plot_final<-plot_grid(evaulation_merge_plot,df_plot,
    widths=15,heights=8,
    labels = "",
    rel_heights = c(1.5,1),
    ncol =1,nrow = 2,label_size = 12,align = "hv")
#save the plot
save.path<-"./manuscript/figures/"
ggsave(paste0(save.path,"FigureS_pmodel_evaluation_with_diff_cali_paras.png"),
       plot_final,width = 20,height = 15)

