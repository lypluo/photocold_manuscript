#######################################################
##Aim: improve p-model performance
##both for the early spring and peak season
#######################################################
#calibrated each site separately to improve the model
#-->after set the iteration to 5000, the model improved substantially
#----------
library(tidyverse)
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
#
# sites<-unique(df_recent$sitename)
# for (i in 1:length(sites)) {
#   df_temp<-df_recent %>%
#     filter(sitename==sites[i])
#   text.Date<-min(df_temp$date)+c(max(df_temp$date)-min(df_temp$date))*0.1
#   df_temp %>%
#     ggplot()+
#     geom_point(aes(x=date,y=gpp))+
#     geom_point(aes(x=date,y=gpp_mod),col="red")+
#     annotate(geom = "text",x=text.Date,y=15,label=sites[i])
# }

#load the data Beni sent me before:
df_old<-read.csv(file=paste0("./data-raw/raw_data/Data_sent_by_Beni/","ddf_fluxnet2015_pmodel_with_forcings_stocker19gmd.csv"))
df_old<-df_old %>%
  mutate(date=lubridate::mdy(date),
         year=lubridate::year(date)) %>%
  na.omit(gpp_obs)
#####
source(paste0("./R/functions_in_model/model_hardening_byBeni_addbaseGDD_rev.R"))
#--------------------------------------------------------------
#(2) retreive the optimized parameter for the selected sites
#--------------------------------------------------------------
# set initial value
par <- c("a" = 0, "b" = 0.5, "c" = 50, "d" = 0.1, "e" = 1,"f"=1,"k"=5)
lower=c(-50,0,0,0,0,0,-10)
upper=c(50,20,200,20,2,2,10)

# run model and compare to true values
# returns the RMSE
cost <- function(
  data,
  par
) {

  scaling_factor <- data %>%
    # group_by(sitename) %>%
    do({
      scaling_factor <- model_hardening_2par(
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
  #mae:mean absolute error:
  # mae<-sum(abs(df$gpp - df$gpp_mod * df$scaling_factor))/nrow(df)
  # This visualizes the process,
  # comment out when running for real
  # plot(df$gpp, type = 'p',ylim=c(0,12))
  # lines(df$gpp_mod, col = "red")
  # lines(df$gpp_mod * df$scaling_factor, col = "blue",cex=1.2)
  # Sys.sleep(0.1)

  return(mse)
}

#--------------------------------------------------------------
#(3) optimize for all sites
#--------------------------------------------------------------
#first load the PFTs information:
#load the modis data-->tidy from Beni
load(paste0("./data-raw/raw_data/sites_info/","Pre_selected_sites_info.RDA"))
sites.info<-df_sites_sel
#
df_recent<-df_recent %>%
  left_join(
    sites.info,
    by = "sitename"
  )

#adding day of the year
library(sirad)  #calculate the day of the year
df_recent$doy<-dayOfYear(df_recent$date)

#main Clim-PFTs
df_recent$Clim_PFTs<-paste0(df_recent$koeppen_code,"-",df_recent$classid)
Clim.PFTs<-sort(unique(df_recent$Clim_PFTs))
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

# optimize for all sites
# library(tictoc)#-->record the parameterization time
# tic("start to parameterize")
# par_allsites<-c()
# # for(i in 1:length(sel_sites)){
#   df_sel<-df_recent
#     # dplyr::filter(sitename==sel_sites[i])
# 
#   optim_par <- GenSA::GenSA(
#   par = par,
#   fn = cost,
#   data = df_sel,
#   lower = lower,
#   upper = upper,
#   control = list(max.call=5000))$par
# 
#   # print(i)
#   par_allsites<-optim_par
# # }
# print("finish parameterization")
# toc()
# #
# # names(par_allsites)<-all_sites
# print(par_allsites)
# # save the optimized data
# save(par_allsites,file = paste0("./data/model_parameters/parameters_MSE_add_baseGDD/","optim_par_run5000_beni_allsite_updated.rds"))

#--------------------------------------------------------------
#(4) compare the gpp_obs, ori modelled gpp, and gpp modelled using optimated parameters
#--------------------------------------------------------------
load(paste0("./data/model_parameters/parameters_MSE_add_baseGDD/","optim_par_run5000_beni_allsite_updated.rds"))
#a.get the stress factor(calibration factor) for each Clim-PFT:using the same parameters
par_Clim_PFTs<-list(a1=par_allsites,a2=par_allsites,a3=par_allsites,a4=par_allsites,
                    a5=par_allsites,a6=par_allsites,a7=par_allsites,a8=par_allsites,
                    a9=par_allsites,a10=par_allsites)
names(par_Clim_PFTs)<-Clim.PFTs

df_final<-c()
for (i in 1:length(Clim.PFTs)) {
  df_sel<-df_recent %>%
    dplyr::filter(Clim_PFTs==Clim.PFTs[i])

  scaling_factors <- df_sel %>%
    # group_by(sitename, year) %>%
    do({
      scaling_factor <- model_hardening_2par(.,par_Clim_PFTs[[i]])
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

#b.make evaluation plots
#!!first need to merge the modelled gpp from different sources:
df_final_new$year<-lubridate::year(df_final_new$date)
df_final_new<-left_join(df_final_new,df_old,by = c("sitename", "date", "year")) %>%
  mutate(gpp_obs_recent=gpp,
         gpp_obs_old=gpp_obs,
         gpp_mod_FULL_ori=gpp_mod_FULL,
         gpp_mod_recent_ori=gpp_mod,
         gpp_mod_recent_optim=gpp_mod*scaling_factor_optim,
         gpp=NULL,
         gpp_obs=NULL,
         gpp_mod=NULL)

#--------------------------
#5.modelled and observed gpp:scatter plots
#-------------------------
#--------
#5a.plot for site
#--------
plot_modobs_general<-c()
df_modobs<-c()
for(i in 1:length(sel_sites)){

  df_modobs_each<-df_final_new %>%
    select(sitename,date,gpp_obs_recent,gpp_mod_FULL_ori,gpp_mod_recent_ori,gpp_mod_recent_optim) %>%
    mutate(gpp_obs=gpp_obs_recent,
           gpp_mod_old_ori=gpp_mod_FULL_ori,
           gpp_mod_recent_ori=gpp_mod_recent_ori,
           gpp_mod_recent_optim=gpp_mod_recent_optim) %>%
    mutate(gpp_obs_recent=NULL,
           gpp_mod_FULL_ori=NULL)
  #
  df_modobs<-rbind(df_modobs,df_modobs_each)

  # #scatter plots to compare the model and observation gpp
  # gpp_modobs_comp1<-df_modobs_each %>%
  #   analyse_modobs2("gpp_mod_old_ori", "gpp_obs", type = "heat")
  # gpp_modobs_comp2<-df_modobs_each %>%
  #   analyse_modobs2("gpp_mod_recent_ori", "gpp_obs", type = "heat")
  # gpp_modobs_comp3<-df_modobs_each %>%
  #   analyse_modobs2("gpp_mod_recent_optim", "gpp_obs", type = "heat")
  # # add the site-name:
  # gpp_modobs_comp1$gg<-gpp_modobs_comp1$gg+
  #   annotate(geom="text",x=15,y=0,label=sel_sites[i])
  # gpp_modobs_comp2$gg<-gpp_modobs_comp2$gg+
  #   annotate(geom="text",x=15,y=0,label=sel_sites[i])
  # gpp_modobs_comp3$gg<-gpp_modobs_comp3$gg+
  #   annotate(geom="text",x=15,y=0,label=sel_sites[i])
  # 
  # #merge two plots
  # evaulation_merge_plot<-plot_grid(gpp_modobs_comp1$gg,
  #                                  gpp_modobs_comp2$gg,gpp_modobs_comp3$gg,
  #                                  widths=15,heights=4,
  #                                  labels = "auto",ncol =3,nrow = 1,label_size = 12,align = "hv")
  # # plot(evaulation_merge_plot)
  # 
  # # put all the plots together:
  # plot_modobs_general[[i]]<-evaulation_merge_plot
}
# names(plot_modobs_general)<-sel_sites

#print the plot
# plot_modobs_general

#(2) For Seasonality
#b. Seasonal course for each sites in different PFTs:
# df_modobs %>%
#   mutate(doy = lubridate::yday(date)) %>%
#   group_by(sitename, doy) %>%
#   summarise(obs = mean(gpp_obs, na.rm = TRUE),
#             mod_old_ori=mean(gpp_mod_old_ori, na.rm = TRUE),
#             mod_recent_ori=mean(gpp_mod_recent_ori, na.rm = TRUE),
#             mod_recent_optim=mean(gpp_mod_recent_optim,na.rm = TRUE)) %>%
#   pivot_longer(c(obs,mod_old_ori,mod_recent_ori,mod_recent_optim), names_to = "Source", values_to = "gpp") %>%
#   ggplot(aes(doy, gpp, color = Source)) +
#   geom_line() +
#   scale_color_manual(values = c("mod_old_ori" = "red","mod_recent_ori"="steelblue2",
#                                 "mod_recent_optim" = "orange", "obs" = "black"),
#                      labels = c("Old P-model","Recent Ori P-model", "Recent Optim P-model","Obs.")) +
#   labs(y = expression( paste("GPP (g C m"^-2, " d"^-1, ")" ) ),
#        x = "Day of year") +
#   facet_wrap(~sitename)
#update using original p-model
season_plot<-df_modobs %>%
  mutate(doy = lubridate::yday(date)) %>%
  group_by(sitename, doy) %>%
  dplyr::summarise(obs = mean(gpp_obs, na.rm = TRUE),
            mod_old_ori=mean(gpp_mod_old_ori, na.rm = TRUE),
            mod_recent_ori=mean(gpp_mod_recent_ori, na.rm = TRUE),
            mod_recent_optim=mean(gpp_mod_recent_optim,na.rm = TRUE)) %>%
  pivot_longer(c(obs,mod_old_ori,mod_recent_ori,mod_recent_optim), names_to = "Source", values_to = "gpp") %>%
  ggplot(aes(doy, gpp, color = Source)) +
  geom_line() +
  scale_color_manual("GPP sources",values = c("mod_old_ori" = "tomato",
                                              "mod_recent_optim" = "green4", "obs" = "gray4"),
                     labels = c("Orig. P-model", "Cali. P-model","Observations")) +
  labs(y = expression( paste("GPP (g C m"^-2, " d"^-1, ")" ) ),
       x = "DoY") +
  facet_wrap(~sitename)+
  theme(
    legend.text = element_text(size=20),
    legend.key.size = unit(2, 'lines'),
    axis.title = element_text(size=24),
    axis.text = element_text(size = 20),
    text = element_text(size=24),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour ="grey",fill="white"),
    # legend.background = element_blank(),
    legend.position = c(0.75,0.05)
  )

####
#save the plot
save.path<-"./manuscript/figures/"
ggsave(paste0(save.path,"FigureS_pmodel_vs_obs_foreachsite_1set_parameter.png"),season_plot,width = 20,height = 20)

#--------
#5b.for Clim-PFTs
#--------
df_modobs<-c()
for(i in 1:length(Clim.PFTs)){
  
  df_modobs_each<-df_final_new %>%
    filter(Clim_PFTs==Clim.PFTs[i]) %>%
    select(sitename,date,Clim_PFTs,gpp_obs_recent,gpp_mod_FULL_ori,gpp_mod_recent_ori,gpp_mod_recent_optim) %>%
    mutate(gpp_obs=gpp_obs_recent,
           gpp_mod_old_ori=gpp_mod_FULL_ori,
           gpp_mod_recent_ori=gpp_mod_recent_ori,
           gpp_mod_recent_optim=gpp_mod_recent_optim) %>%
    mutate(gpp_obs_recent=NULL,
           gpp_mod_FULL_ori=NULL)
  #
  df_modobs<-rbind(df_modobs,df_modobs_each)
  
  #scatter plots to compare the model and observation gpp
  # gpp_modobs_comp1<-df_modobs_each %>%
  #   analyse_modobs2("gpp_mod_old_ori", "gpp_obs", type = "heat")
  # gpp_modobs_comp2<-df_modobs_each %>%
  #   analyse_modobs2("gpp_mod_recent_ori", "gpp_obs", type = "heat")
  # gpp_modobs_comp3<-df_modobs_each %>%
  #   analyse_modobs2("gpp_mod_recent_optim", "gpp_obs", type = "heat")
  # # add the site-name:
  # gpp_modobs_comp1$gg<-gpp_modobs_comp1$gg+
  #   annotate(geom="text",x=15,y=0,label=Clim.PFTs[i])
  # gpp_modobs_comp2$gg<-gpp_modobs_comp2$gg+
  #   annotate(geom="text",x=15,y=0,label=Clim.PFTs[i])
  # gpp_modobs_comp3$gg<-gpp_modobs_comp3$gg+
  #   annotate(geom="text",x=15,y=0,label=Clim.PFTs[i])
  # 
  # #merge two plots
  # evaulation_merge_plot<-plot_grid(gpp_modobs_comp1$gg,
  #                                  gpp_modobs_comp2$gg,gpp_modobs_comp3$gg,
  #                                  widths=15,heights=4,
  #   labels = "auto",ncol =3,nrow = 1,label_size = 12,align = "hv")
  # # plot(evaulation_merge_plot)
  # 
  # # put all the plots together:
  # plot_modobs_general[[i]]<-evaulation_merge_plot
}

#(2) For Seasonality
season_plot<-df_modobs %>%
  mutate(doy = lubridate::yday(date)) %>%
  group_by(Clim_PFTs, doy) %>%
  dplyr::summarise(obs = mean(gpp_obs, na.rm = TRUE),
                   mod_old_ori=mean(gpp_mod_old_ori, na.rm = TRUE),
                   mod_recent_ori=mean(gpp_mod_recent_ori, na.rm = TRUE),
                   mod_recent_optim=mean(gpp_mod_recent_optim,na.rm = TRUE)) %>%
  pivot_longer(c(obs,mod_old_ori,mod_recent_optim), names_to = "Source", values_to = "gpp") %>%
  ggplot(aes(doy, gpp, color = Source)) +
  geom_line() +
  scale_color_manual("GPP sources",values = c("mod_old_ori" = "tomato",
                                              "mod_recent_optim" = "green4", "obs" = "gray4"),
                     labels = c("Orig. P-model", "Cali. P-model","Obseravations")) +
  labs(y = expression( paste("GPP (g C m"^-2, " d"^-1, ")" ) ),
       x = "DoY") +
  # annotate(geom="text",x=200,y=2,label="")+
  facet_wrap(~Clim_PFTs)+
  theme(
    legend.text = element_text(size=20),
    legend.key.size = unit(2, 'lines'),
    axis.title = element_text(size=24),
    axis.text = element_text(size = 20),
    text = element_text(size=24),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour ="grey",fill="white"),
    # legend.background = element_blank(),
    legend.position = c(0.75,0.1)
  )
#save the plot
save.path<-"./manuscript/figures/"
ggsave(paste0(save.path,"Figure5_pmodel_vs_obs_forClimPFTs_1set_parameter.png"),
       season_plot,width = 15,height = 10)

