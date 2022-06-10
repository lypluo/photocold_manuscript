#######################################################
##Aim: improve p-model performance
##both for the early spring and peak season
##update in June, 2022:parameterization for the different PFTs
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
#(3) optimize for each PFT
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
#main PFTs
PFTs<-sort(unique(df_merge$classid))

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
# optimize for each PFT-->parameterized before, hence used the parameters have been parameterized 
# library(tictoc)#-->record the parameterization time
# tic("start to parameterize")
# par_Clim_PFTs<-c()
# for(i in 1:length(PFTs)){
#   df_sel<-df_merge.new %>%
#     dplyr::filter(PFTs==PFTs[i])
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
# save(par_Clim_PFTs,file = paste0(base.path,"data/parameters_MSE_add_baseGDD/test/","optim_par_run5000_beni_PFTs_update.rds"))

#--------------------------------------------------------------
#(5) compare the gpp_obs, ori modelled gpp, and gpp modelled using optimated parameters
#--------------------------------------------------------------
#load model parameters
load(paste0("./data/model_parameters/parameters_MSE_add_baseGDD/","optim_par_run5000_beni_PFTs.rds"))
#check par_Clim_PFTs
print(par_PFTs)
#a.get the stress factor(calibration factor) for each Clim_PFT
df_final<-c()
for (i in 1:length(PFTs)) {
  df_sel<-df_merge.new %>%
    dplyr::filter(classid==PFTs[i])

  scaling_factors <- df_sel %>%
    # group_by(sitename, year) %>%
    do({
      scaling_factor <- model_hardening_2par(.,par_PFTs[[i]])
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

### make evaluation plots
#(1) For General plots
devtools::load_all("D:/Github/rbeni/")
library(rbeni) #-->make the evaluation plot
library(cowplot)
library(grid)

#--------------------------
#modelled and observed gpp:scatter plots:scatter plots-->errors occur, need to find out the reasons later (2022-May,04)
#-------------------------
# only selecting out the Cfb-DBF and Dfb-DBF
#--------
#5a.Cfb-DBF
#--------
season_plot<-df_final_new %>%
  filter(Clim_PFTs=="Cfb-DBF")%>%
  select(sitename,date,classid,gpp_obs_recent,gpp_mod_FULL_ori,gpp_mod_recent_ori,gpp_mod_recent_optim) %>%
  mutate(gpp_obs=gpp_obs_recent,
         gpp_mod_old_ori=gpp_mod_FULL_ori,
         gpp_mod_recent_ori=gpp_mod_recent_ori,
         gpp_mod_recent_optim=gpp_mod_recent_optim) %>%
  mutate(gpp_obs_recent=NULL, ##remove the variables
         gpp_mod_FULL_ori=NULL)%>%
  mutate(doy = lubridate::yday(date)) %>%
  group_by(sitename,doy)%>%
  summarise(obs = mean(gpp_obs, na.rm = TRUE),
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
  annotate(geom="text",x=200,y=2,label="")+
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
    legend.position = c(0.5,0.3)
  )

#save the plot
save.path<-"./manuscript/test_files/Diff_parameterization_approach/"
ggsave(paste0(save.path,"Figure5_pmodel_vs_obs_forPFTs_3sets_parameter_check_Cfb-DBF.png"),
       season_plot,width = 15,height = 10)


#--------
#5b.for Dfb-DBF
#--------
season_plot<-df_final_new %>%
  filter(Clim_PFTs=="Dfb-DBF")%>%
  select(sitename,date,classid,gpp_obs_recent,gpp_mod_FULL_ori,gpp_mod_recent_ori,gpp_mod_recent_optim) %>%
  mutate(gpp_obs=gpp_obs_recent,
         gpp_mod_old_ori=gpp_mod_FULL_ori,
         gpp_mod_recent_ori=gpp_mod_recent_ori,
         gpp_mod_recent_optim=gpp_mod_recent_optim) %>%
  mutate(gpp_obs_recent=NULL, ##remove the variables
         gpp_mod_FULL_ori=NULL)%>%
  mutate(doy = lubridate::yday(date)) %>%
  group_by(sitename,doy)%>%
  summarise(obs = mean(gpp_obs, na.rm = TRUE),
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
  annotate(geom="text",x=200,y=2,label="")+
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
    legend.position = c(0.8,0.3)
  )

#save the plot
save.path<-"./manuscript/test_files/Diff_parameterization_approach/"
ggsave(paste0(save.path,"Figure5_pmodel_vs_obs_forPFTs_3sets_parameter_check_Dfb-DBF.png"),
       season_plot,width = 15,height = 10)
