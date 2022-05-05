#######################################################
##Aim: improve p-model performance
##both for the early spring and peak season
#######################################################
#update in March,07:Normlize the gpp and use two parametes scaling factor( e and f)
#----------
library(dplyr)
library(tidyverse)
library(GenSA)
library(lubridate)
#-------------------------
#(1)load the data and hardening funciton
#-------------------------
#####
#load the data updated p-model data uploaded by Koen
df_recent <- readRDS(paste0("./data-raw/raw_data/P_model_output/model_data.rds")) %>%
  mutate(
    year = format(date, "%Y")
  ) %>%
  na.omit()

sites<-unique(df_recent$sitename)
for (i in 1:length(sites)) {
  df_temp<-df_recent %>%
    filter(sitename==sites[i])
  text.Date<-min(df_temp$date)+c(max(df_temp$date)-min(df_temp$date))*0.1
  #
  df_plot<-df_temp %>%
    ggplot()+
    geom_point(aes(x=date,y=gpp))+
    geom_point(aes(x=date,y=gpp_mod),col="red")+
    annotate(geom = "text",x=text.Date,y=15,label=sites[i])
  #
  k=quantile(df_temp$gpp,probs = c(0.01,0.05,seq(0.1,0.9,0.1)),0.95,0.99)
  df_plot+
    geom_hline(yintercept = k,col="blue")
}

# #filter the observational gpp data:
# df_recent_new<-c()
# for (i in 1:length(sites)) {
#   df_temp<-df_recent %>%
#     filter(sitename==sites[i])
#   #filter the gpp observation data(remove the gpp that below 5 percentile and negative):
#   k=quantile(df_temp$gpp,probs = c(0.01,0.05,seq(0.1,0.9,0.1)),0.95,0.99)
#   df_temp$gpp[df_temp$gpp<as.numeric(k[2])& df_temp$gpp<0]<-NA
#   # df_temp %>%
#   #   ggplot()+
#   #   geom_point(aes(x=date,y=gpp))+
#   #   geom_point(aes(x=date,y=gpp_mod),col="red")+
#   #   annotate(geom = "text",x=text.Date,y=15,label=sites[i])
#   df_recent_new<-rbind(df_recent_new,df_temp)
# }

#do not filter the observation gpp in study:
df_recent_new<-df_recent

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
#(3) optimize for each Clim-PFT
#--------------------------------------------------------------
#first load the PFTs information:
#load the modis data-->tidy from Beni
load(paste0("./data-raw/raw_data/sites_info/","Pre_selected_sites_info.RDA"))
sites.info<-df_sites_sel
#
df_merge<-df_recent_new %>%
left_join(
    sites.info,
    by = "sitename"
  )
#day of the year
library(sirad)  #calculate the day of the year
df_merge$doy<-dayOfYear(df_merge$date)

#main Clim-PFTs
df_merge$Clim_PFTs<-paste0(df_merge$koeppen_code,"-",df_merge$classid)
Clim.PFTs<-sort(unique(df_merge$Clim_PFTs))

#------------------------------------------
#(4)normalized the GPP-->for each site,
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

#---------------------------------
# optimize for each Clim.-PFT
# library(tictoc)#-->record the parameterization time
# tic("start to parameterize")
# par_Clim_PFTs<-c()
# for(i in 1:length(Clim.PFTs)){
#   df_sel<-df_merge.new %>%
#     dplyr::filter(Clim_PFTs==Clim.PFTs[i])
#
#   optim_par <- GenSA::GenSA(
#   par = par,
#   fn = cost,
#   data = df_sel,
#   lower = lower,
#   upper = upper,
#   control = list(max.call=5000))$par
#
#   print(i)
#   par_Clim_PFTs[[i]]<-optim_par
# }
# print("finish parameterization")
# toc()
#
# names(par_Clim_PFTs)<-Clim.PFTs
# print(par_Clim_PFTs)
# # save the optimized data
# save(par_Clim_PFTs,file = paste0(base.path,"data/parameters_MSE_add_baseGDD/test/","optim_par_run5000_beni_Clim_andPFTs_update.rds"))

#--------------------------------------------------------------
#(5) compare the gpp_obs, ori modelled gpp, and gpp modelled using optimated parameters
#--------------------------------------------------------------
#load model parameters
load(paste0("./data/model_parameters/parameters_MSE_add_baseGDD/","optim_par_run5000_beni_Clim_andPFTs_update.rds"))
#check par_Clim_PFTs
print(par_Clim_PFTs)
#a.get the stress factor(calibration factor) for each PFT
df_final<-c()
for (i in 1:length(Clim.PFTs)) {
  df_sel<-df_merge.new %>%
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

# need to remove the sites that do not used in this analysis:
rm.sites<-c("BE-Bra","CA-SF1","CA-SF2","FI-Sod","US-Wi4")
df_final_new<-df_final_new %>%
  filter(sitename!=rm.sites[1] & sitename!=rm.sites[2]&sitename!=rm.sites[3]&sitename!=rm.sites[4]&sitename!=rm.sites[5])

###########test for ts of temp,tmin and tmax############
#Ta
df_final_new %>%
  mutate(doy = lubridate::yday(date)) %>%
  group_by(Clim_PFTs,doy) %>%
  # group_by(sitename,doy) %>%
  dplyr::summarise(gpp_obs=mean(gpp_obs_recent,na.rm=T),
          mean_Ta=mean(temp,na.rm=T),
          mean_Tmin=mean(tmin,na.rm=T),
          mean_Tmax=mean(tmax,na.rm=T),
          VPD=mean(vpd,na.rm=T),
          mean_prec=mean(prec,na.rm=T))%>%
  pivot_longer(c(mean_Ta,mean_Tmin,mean_Tmax),
               names_to = "Ta_source",values_to = "Ta") %>%
  ggplot(aes(doy,Ta,color = Ta_source))+
  geom_line()+
  facet_grid(~Clim_PFTs)
  # facet_grid(~sitename)
#VPD
df_final_new %>%
  mutate(doy = lubridate::yday(date)) %>%
  group_by(Clim_PFTs,doy) %>%
  dplyr::summarise(VPD=mean(vpd,na.rm=T))%>%
  ggplot(aes(doy,VPD))+
  geom_line()+
  facet_grid(~Clim_PFTs)
#prec
df_final_new %>%
  mutate(doy = lubridate::yday(date)) %>%
  group_by(Clim_PFTs,doy) %>%
  dplyr::summarise(mean_prec=mean(prec,na.rm=T))%>%
  ggplot(aes(doy,mean_prec))+
  geom_line()+
  facet_grid(~Clim_PFTs)
#ppfd
df_final_new %>%
  mutate(doy = lubridate::yday(date)) %>%
  group_by(Clim_PFTs,doy) %>%
  dplyr::summarise(mean_ppfd=mean(ppfd,na.rm=T))%>%
  ggplot(aes(doy,mean_ppfd))+
  geom_line()+
  facet_grid(~Clim_PFTs)
#fapar
df_final_new %>%
  mutate(doy = lubridate::yday(date)) %>%
  group_by(Clim_PFTs,doy) %>%
  dplyr::summarise(mean_fapar=mean(fapar_itpl,na.rm=T))%>%
  ggplot(aes(doy,mean_fapar))+
  geom_line()+
  facet_grid(~Clim_PFTs)
#gpp
df_final_new %>%
  mutate(doy = lubridate::yday(date)) %>%
  group_by(Clim_PFTs,doy) %>%
  dplyr::summarise(gpp_obs=mean(gpp_obs_recent,na.rm=T),
            gpp_mod_old=mean(gpp_mod_FULL_ori,na.rm=T),
            gpp_mod_new=mean(gpp_mod_recent_ori,na.rm=T))%>%
  pivot_longer(c(gpp_obs,gpp_mod_old,gpp_mod_new),
               names_to = "gpp_source",values_to = "gpp")%>%
  ggplot(aes(doy,gpp,color=gpp_source))+
  geom_line()+
  facet_grid(~Clim_PFTs)
##through gpp--> the biggest mismatch between obs and mod are in Cfb-DBF(underestimation):
#checking these sites
unique(df_final_new[df_final_new$Clim_PFTs=="Cfb-DBF",]$sitename)

### make evaluation plots
#(1) For General plots
devtools::load_all("D:/Github/rbeni/")
library(rbeni) #-->make the evaluation plot
library(cowplot)
library(grid)

#--------------------------
#modelled and observed gpp:scatter plots:scatter plots-->errors occur, need to find out the reasons later (2022-May,04)
#-------------------------
# plot_modobs_general<-c()
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
# names(plot_modobs_general)<-Clim.PFTs

#print the plot
# plot_modobs_general

#(2) For Seasonality
#a. Seasonal course for different Clim_PFTs:
#plotting:
# 
# season_plot<-df_modobs %>%
#   mutate(doy = lubridate::yday(date)) %>%
#   group_by(Clim_PFTs, doy) %>%
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
#   annotate(geom="text",x=200,y=2,label="")+
#   facet_wrap(~Clim_PFTs)

#update using original p-model
##
nsites<-df_modobs %>%
  group_by(Clim_PFTs)%>%
  dplyr::summarise(nsite=length(unique(sitename)))
nsites$label<-paste0("N = ",nsites$nsite)
sites_num.info<-data.frame(
  doy=rep(20,nrow(nsites)),
  gpp=rep(14,nrow(nsites)),
  nsites
)
#
tag_facet <- function(p, open = "", close = "", tag_pool = letters, x = -Inf, y = Inf, 
                      hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
  
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
  p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE) 
}

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
                     labels = c("Orig. P-model", "Cali. P-model","EC based")) +
  labs(y = expression( paste("GPP (g C m"^-2, " d"^-1, ")" ) ),
       x = "Day of year") +
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
#

#print the plot
library(egg)
season_plot_new<-tag_facet(season_plot,x=sites_num.info$doy,y=sites_num.info$gpp,
          tag_pool = sites_num.info$label,size=5)
# annotate(geom = "text",x=sites_num.info$x,
  #          y=sites_num.info$y,label=sites_num.info$label)
#save the plot
save.path<-"./manuscript/figures/"
ggsave(paste0(save.path,"Figure6_pmodel_vs_obs_forClimPFTs.png"),season_plot_new,width = 15,height = 10)


##########################################################################
#b. Seasonal course for each sites in different PFTs:
# For DBF:
df_modobs %>%
  filter(Clim_PFTs=="Cfa-DBF") %>%
  mutate(doy = lubridate::yday(date)) %>%
  group_by(sitename, doy) %>%
  dplyr::summarise(obs = mean(gpp_obs, na.rm = TRUE),
            mod_old_ori=mean(gpp_mod_old_ori, na.rm = TRUE),
            mod_recent_ori=mean(gpp_mod_recent_ori, na.rm = TRUE),
            mod_recent_optim=mean(gpp_mod_recent_optim,na.rm = TRUE)) %>%
  pivot_longer(c(obs,mod_old_ori,mod_recent_ori,mod_recent_optim), names_to = "Source", values_to = "gpp") %>%
  ggplot(aes(doy, gpp, color = Source)) +
  geom_line() +
  scale_color_manual(values = c("mod_old_ori" = "red","mod_recent_ori"="steelblue2",
                                "mod_recent_optim" = "orange", "obs" = "black"),
                     labels = c("Old P-model","Recent Ori P-model", "Recent Optim P-model","Obs.")) +
  labs(y = expression( paste("GPP (g C m"^-2, " d"^-1, ")" ) ),
       x = "Day of year") +
  facet_wrap(~sitename)

####
library(dplyr)
df_modobs %>%
  filter(Clim_PFTs=="Cfb-DBF")%>%
  group_by(sitename)%>%
  dplyr::summarise(n=n(),
            mean_gpp_obs=mean(gpp_obs),
            mean_gpp_recent_optim=mean(gpp_mod_recent_optim))
