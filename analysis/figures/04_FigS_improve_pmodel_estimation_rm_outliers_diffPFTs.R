#######################################################
##Aim: improve p-model performance
##both for the early spring and peak season
#######################################################
#At the end, adopt the Mekela et al., 2008 to improve the model representation
#----------
library(tidyverse)
library(dplyr)
library(GenSA)
library(lubridate)
# remotes::install_github("computationales/ingestr") #install the package
library(ingestr)
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

sites<-unique(df_recent$sitename)
#
library(ggplot2)
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

#filter the observational gpp data:
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
##------------------
plot_temp<-df_old %>%
  mutate(doy=lubridate::yday(date))%>%
  select(sitename,doy,gpp_mod_FULL,gpp_mod_ORG)%>%
  group_by(sitename,doy)%>%
  summarise(gpp_mod_ORG=mean(gpp_mod_ORG),gpp_mod_FULL=mean(gpp_mod_FULL))%>%
  pivot_longer(c(gpp_mod_FULL,gpp_mod_ORG),values_to = "GPP",names_to = "Source")%>%
  ggplot(aes(x=doy,y=GPP,col=Source))+
  geom_point()+
  facet_wrap(.~sitename)
temp.path<-"./manuscript/test_files/"
ggsave(plot_temp,filename = paste0(temp.path,"gpp_Pmodel_GPP_Full_vs_GPP_ORG.pdf"),width = 15,height = 10)
#####
# source(paste0("./R/functions_in_model/model_hardening_byBeni_addbaseGDD_rev.R"))
source(paste0("./R/functions_in_model/newly_formulated_fun/model_fT_rev.R"))
#--------------------------------------------------------------
#(2) retreive the optimized parameter for the selected sites
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
#(3) optimize for each site
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
##adding day of the year
df_merge$doy<-yday(df_merge$date)

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

# need to remove the sites that do not used in this analysis:
rm.sites<-c("BE-Bra","CA-SF1","CA-SF2","FI-Sod","US-Wi4")
df_merge_new<-df_merge.new %>%
  filter(sitename!=rm.sites[1] & sitename!=rm.sites[2]&sitename!=rm.sites[3]&sitename!=rm.sites[4]&sitename!=rm.sites[5])

#---------------------------------
# optimize for each PFT
#main PFTs
PFTs<-unique(df_merge.new$classid)
# optimize for each PFT
# library(tictoc)#-->record the parameterization time
# tic("start to parameterize")
# par_PFTs<-c()
# for(i in 1:length(PFTs)){
#   df_sel<-df_merge.new %>%
#     dplyr::filter(classid==PFTs[i])
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
#   par_PFTs[[i]]<-optim_par
# }
# print("finish parameterization")
# toc()
# 
# names(par_PFTs)<-PFTs
# print(par_PFTs)
# # save the optimized data
# save(par_PFTs,file = paste0("data/model_parameters/parameters_MAE_newfT/","optim_par_run5000_PFTs.rds"))

#--------------------------------------------------------------
#(5) compare the gpp_obs, ori modelled gpp, and gpp modelled using optimated parameters
#--------------------------------------------------------------
load(paste0("./data/model_parameters/parameters_MAE_newfT/","optim_par_run5000_PFTs.rds"))
#a.get the stress factor(calibration factor) for each PFT
df_final<-c()
for (i in 1:length(PFTs)) {
  df_sel<-df_merge.new %>%
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
  df_final<-rbind(df_final,df_sel_new)
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
df_merge_new<-left_join(df_final_new,df_old,by = c("sitename", "date", "year")) %>%
  mutate(gpp_obs_recent=gpp,
         gpp_obs_old=gpp_obs,
         gpp_mod_FULL_ori=gpp_mod_FULL,
         gpp_mod_recent_ori=gpp_mod,
         gpp_mod_recent_optim=gpp_mod*scaling_factor_optim,
         gpp=NULL,
         gpp_obs=NULL,
         gpp_mod=NULL)
#
###########test for ts of temp,tmin and tmax############
#Ta.....
df_merge_new %>%
  mutate(doy = lubridate::yday(date)) %>%
  group_by(classid,doy) %>%
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
  facet_grid(~classid)
  # facet_grid(~sitename)

### make evaluation plots
#(1) For General plots
# devtools::load_all("D:/Github/rbeni/")
# library(rbeni) #-->make the evaluation plot
library(rbeni)
library(cowplot)
library(grid)

#--------------------------
#modelled and observed gpp:scatter plots-->errors occur, need to find out the reasons later (2022-May,04)
#-------------------------
#source the function:analyse_modobs2
# source(file = "./R/functions_from_beni/analyse_modobs2.R")
plot_modobs_general<-c()
df_modobs<-c()
for(i in 1:length(PFTs)){
  #
  df_modobs_each<-df_merge_new %>%
    filter(classid==PFTs[i]) %>%
    select(sitename,date,classid,gpp_obs_recent,gpp_mod_FULL_ori,gpp_mod_recent_ori,gpp_mod_recent_optim) %>%
    mutate(gpp_obs=gpp_obs_recent,
           gpp_mod_old_ori=gpp_mod_FULL_ori,
           gpp_mod_recent_ori=gpp_mod_recent_ori,
           gpp_mod_recent_optim=gpp_mod_recent_optim) %>%
    mutate(gpp_obs_recent=NULL,
           gpp_mod_FULL_ori=NULL)
  #
  df_modobs<-rbind(df_modobs,df_modobs_each)

  #some setting-->only compare when both gpp_obs and gpp_mod_old_ori are available:
  # df_modobs_each<-df_modobs_each[!is.nan(df_modobs_each$gpp_mod_old_ori) & !is.nan(df_modobs_each$gpp_obs),]
  # #scatter plots to compare the model and observation gpp
  # gpp_modobs_comp1<-df_modobs_each %>%
  #   analyse_modobs2("gpp_mod_old_ori", "gpp_obs", type = "heat")
  # # gpp_modobs_comp2<-df_modobs_each %>%
  # #   analyse_modobs2("gpp_mod_recent_ori", "gpp_obs", type = "heat")
  # gpp_modobs_comp3<-df_modobs_each %>%
  #   analyse_modobs2("gpp_mod_recent_optim", "gpp_obs", type = "heat")
  # # add the site-name:
  # gpp_modobs_comp1$gg<-gpp_modobs_comp1$gg+
  #   annotate(geom="text",x=15,y=0,label=PFTs[i])
  # gpp_modobs_comp2$gg<-gpp_modobs_comp2$gg+
  #   annotate(geom="text",x=15,y=0,label=PFTs[i])
  # gpp_modobs_comp3$gg<-gpp_modobs_comp3$gg+
  #   annotate(geom="text",x=15,y=0,label=PFTs[i])
  # 
  # #merge two plots
  # if(PFTs[i]=="ENF"){
  #   p1<-gpp_modobs_comp1$gg+
  #     ylab(expression("GPP"[EC]*" (g C m"^-2*" d"^-1*")"))+
  #     xlab(expression("GPP"[Orig.P-model]*" (g C m"^-2*" d"^-1*")"))+
  #     annotate(geom="text",x=0,y=20,label="e")
  # 
  #   p3<-gpp_modobs_comp3$gg+
  #     ylab("")+
  #     xlab(expression("GPP"[Cali.P-model]*" (g C m"^-2*" d"^-1*")"))+
  #     annotate(geom="text",x=0,y=20,label="f")
  # }
  # if(PFTs[i]=="DBF"){
  #   p1<-gpp_modobs_comp1$gg+
  #     ylab(expression("GPP"[EC]*" (g C m"^-2*" d"^-1*")"))+
  #     xlab("")+
  #     annotate(geom="text",x=0,y=20,label="a")
  # 
  #   p3<-gpp_modobs_comp3$gg+
  #     ylab("")+
  #     xlab("")+
  #     annotate(geom="text",x=0,y=20,label="b")
  # }
  # if(PFTs[i]=="MF"){
  #   p1<-gpp_modobs_comp1$gg+
  #     ylab(expression("GPP"[EC]*" (g C m"^-2*" d"^-1*")"))+
  #     xlab("")+
  #     annotate(geom="text",x=0,y=20,label="c")
  # 
  #   p3<-gpp_modobs_comp3$gg+
  #     ylab("")+
  #     xlab("")+
  #     annotate(geom="text",x=0,y=20,label="d")
  # }
  # 
  # evaulation_merge_plot<-plot_grid(p1,p3,
  #                                  widths=15,heights=4,
  #                                  ncol =2,nrow = 1,label_size = 12,align = "hv")
  # # plot(evaulation_merge_plot)
  # 
  # #put all the plots together:
  # plot_modobs_general[[i]]<-evaulation_merge_plot
}
# names(plot_modobs_general)<-PFTs

#print the plot
# p_scatterplot_all<-plot_grid(plot_modobs_general$DBF,
#                              plot_modobs_general$MF,
#                              plot_modobs_general$ENF,nrow = 3)
#
# tag_facet <- function(p, open = "", close = "", tag_pool = letters, x = -Inf, y = Inf,
#                       hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
#
#   gb <- ggplot_build(p)
#   lay <- gb$layout$layout
#   tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
#   p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust,
#                 vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE)
# }
# library(egg)
# comp.labels<-data.frame(
#   x=rep(1,6),
#   y=rep(15,6),
#   label=c("a","b","c","d","e","f")
# )
# p_scatterplot_all<-tag_facet(p_scatterplot_all,x=comp.labels$x,y=comp.labels$y,
#                            tag_pool = comp.labels$label,size=5)
# print(p_scatterplot_all)
#save the plot
# save.path<-"./manuscript/figures/"
# ggsave(paste0(save.path,"FigureS_pmodel_vs_obs_scatterplot.png"),p_scatterplot_all,width = 10,height = 15)

#(2) For Seasonality
#a. Seasonal course for different PFTs:
#plotting:
# season_plot<-df_modobs %>%
#   mutate(doy = lubridate::yday(date)) %>%
#   group_by(classid, doy) %>%
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
#   facet_wrap(~classid)

#update using original p-model
##
nsites<-df_modobs %>%
  group_by(classid)%>%
  dplyr::summarise(nsite=length(unique(sitename)))
nsites$label<-paste0("N = ",nsites$nsite)
sites_num.info<-data.frame(
  doy=rep(20,nrow(nsites)),
  gpp=rep(14,nrow(nsites)),
  nsites
)
#
season_plot<-df_modobs %>%
  mutate(doy = lubridate::yday(date)) %>%
  group_by(classid, doy) %>%
  dplyr::summarise(obs = mean(gpp_obs, na.rm = TRUE),
            mod_old_ori=mean(gpp_mod_old_ori, na.rm = TRUE),
            mod_recent_ori=mean(gpp_mod_recent_ori, na.rm = TRUE),
            mod_recent_optim=mean(gpp_mod_recent_optim,na.rm = TRUE)) %>%
  pivot_longer(c(obs,mod_old_ori,mod_recent_optim), names_to = "Source", values_to = "gpp") %>%
  ggplot(aes(doy, gpp, color = Source)) +
  geom_line() +
  scale_color_manual("GPP sources",values = c("mod_old_ori" = "tomato",
                                              "mod_recent_optim" = "steelblue2", "obs" = "gray4"),
                     labels = c("Orig. P-model", "Cali. P-model","Observations")) +
  labs(y = expression( paste("GPP (g C m"^-2, " d"^-1, ")" ) ),
       x = "DoY") +
  # annotate(geom="text",x=200,y=2,label="")+
  facet_wrap(~classid)+
  theme(
    legend.text = element_text(size=20),
    legend.key.size = unit(2, 'lines'),
    axis.title = element_text(size=24),
    axis.text = element_text(size = 20),
    text = element_text(size=24),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour ="grey",fill="white"),
    legend.position = "bottom"
  )
#print the plot
tag_facet <- function(p, open = "", close = "", tag_pool = letters, x = -Inf, y = Inf, 
                      hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
  
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
  p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE) 
}
library(egg)
season_plot_new<-tag_facet(season_plot,x=sites_num.info$doy,y=sites_num.info$gpp,
                           tag_pool = sites_num.info$label,size=5)
# annotate(geom = "text",x=sites_num.info$x,
#          y=sites_num.info$y,label=sites_num.info$label)
#save the plot
save.path<-"./manuscript/figures/"
ggsave(paste0(save.path,"FigureS_pmodel_vs_obs_forPFTs_fT.png"),season_plot,width = 15,height = 10)

#b. Seasonal course for each sites in different PFTs:
# For DBF:
df_modobs %>%
  filter(classid=="DBF") %>%
  mutate(doy = lubridate::yday(date)) %>%
  group_by(sitename, doy) %>%
  dplyr::summarise(obs = mean(gpp_obs, na.rm = TRUE),
            mod_old_ori=mean(gpp_mod_old_ori, na.rm = TRUE),
            mod_recent_ori=mean(gpp_mod_recent_ori, na.rm = TRUE),
            mod_recent_optim=mean(gpp_mod_recent_optim,na.rm = TRUE)) %>%
  pivot_longer(c(obs,mod_old_ori,mod_recent_ori,mod_recent_optim), names_to = "Source", values_to = "gpp") %>%
  ggplot(aes(doy, gpp, color = Source)) +
  geom_line() +
  scale_color_manual(values = c("mod_old_ori" = "red","mod_recent_ori"="dodgerblue",
                                "mod_recent_optim" = "orange", "obs" = "black"),
                     labels = c("Old P-model","Recent Ori P-model", "Recent Optim P-model","Obs.")) +
  labs(y = expression( paste("GPP (g C m"^-2, " d"^-1, ")" ) ),
       x = "Day of year") +
  facet_wrap(~sitename)

# - for MF
df_modobs %>%
  filter(classid=="MF") %>%
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

# - for ENF
df_modobs %>%
  filter(classid=="ENF") %>%
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

#------------------------------------------
#Additional plots(2022-06):compare the modelled GPP with optimilized 3 pars(for diff PFTs)
#with observed GPP-->for different Clim-PFTs
#------------------------------------------
#a.get the stress factor(calibration factor) for each Clim-PFT:using the same parameters 
#for each 
#
test<-df_merge_new %>%
  select(sitename,date,doy,classid,Clim_PFTs,gpp_obs_recent,gpp_mod_FULL_ori,gpp_mod_recent_ori,gpp_mod_recent_optim) %>%
  mutate(gpp_obs=gpp_obs_recent,
         gpp_mod_old_ori=gpp_mod_FULL_ori,
         gpp_mod_recent_ori=gpp_mod_recent_ori,
         gpp_mod_recent_optim=gpp_mod_recent_optim) %>%
  mutate(gpp_obs_recent=NULL,
         gpp_mod_FULL_ori=NULL)
#for Clim-PFTs
season_plot<-test %>%
  group_by(Clim_PFTs, doy) %>%
  dplyr::summarise(obs = mean(gpp_obs, na.rm = TRUE),
                   mod_old_ori=mean(gpp_mod_old_ori, na.rm = TRUE),
                   mod_recent_ori=mean(gpp_mod_recent_ori, na.rm = TRUE),
                   mod_recent_optim=mean(gpp_mod_recent_optim,na.rm = TRUE)) %>%
  pivot_longer(c(obs,mod_old_ori,mod_recent_optim), names_to = "Source", values_to = "gpp") %>%
  ggplot(aes(doy, gpp, color = Source)) +
  geom_line() +
  scale_color_manual("GPP sources",values = c("mod_old_ori" = "tomato",
                                              "mod_recent_optim" = "dodgerblue", "obs" = "gray4"),
                     labels = c("Orig. P-model", "Cali. P-model","Observations")) +
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
    legend.position = "bottom"
  )
##adding the site numbers in each category:
nsites<-test %>%
  group_by(Clim_PFTs)%>%
  dplyr::summarise(nsite=length(unique(sitename)))
nsites$label<-paste0("N = ",nsites$nsite)
sites_num.info<-data.frame(
  doy=rep(20,nrow(nsites)),
  gpp=rep(14,nrow(nsites)),
  nsites
)
season_plot_new<-tag_facet(season_plot,x=sites_num.info$doy,y=sites_num.info$gpp,
                           tag_pool = sites_num.info$label,size=5)

#save the plot
save.path<-"./manuscript/figures/"
ggsave(paste0(save.path,"Figure4_pmodel_vs_obs_forClimPFTs_3set_parameter_fT.png"),
       season_plot_new,width = 15,height = 10)

#for different sites
test %>%
  filter(Clim_PFTs=="Cfb-MF")%>%
  group_by(sitename, doy) %>%
  dplyr::summarise(obs = mean(gpp_obs, na.rm = TRUE),
                   mod_old_ori=mean(gpp_mod_old_ori, na.rm = TRUE),
                   mod_recent_ori=mean(gpp_mod_recent_ori, na.rm = TRUE),
                   mod_recent_optim=mean(gpp_mod_recent_optim,na.rm = TRUE)) %>%
  pivot_longer(c(obs,mod_old_ori,mod_recent_optim), names_to = "Source", values_to = "gpp") %>%
  ggplot(aes(doy, gpp, color = Source)) +
  geom_line() +
  scale_color_manual("GPP sources",values = c("mod_old_ori" = "red",
                                              "mod_recent_optim" = "dodgerblue", "obs" = "gray4"),
                     labels = c("Orig. P-model", "Cali. P-model","Observations")) +
  labs(y = expression( paste("GPP (g C m"^-2, " d"^-1, ")" ) ),
       x = "DoY") +
  # annotate(geom="text",x=200,y=2,label="")+
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
    legend.position = "bottom"
  )

