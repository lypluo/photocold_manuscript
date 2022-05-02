#################################
#Aim: Try to use the constant LUE to estimate the GPP
#constant LUE is estimated by using the data outside the green-up period
#-->to compare if constant LUE also overestimated GPP in the early spring!
#################################
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
##add robust regression:
library(robustbase)
#--------------
#(1)load the merged fluxes data 
#and also used the phenophases(SOS and POS) extracted in "Using_sites_in_Fluxnet2015" analysis 
#to determine the green-up period 
#--------------
#--load the merged fluxes data
load.path<-"D:/data/photocold_project/Merge_Data/Using_sites_in_Fluxnet2015_constant_LUEmodel/"
load(paste0(load.path,"Merged_Flux_and_PhenoCam.RDA"))
#
df_merge$year<-year(df_merge$date)
df_merge$doy<-yday(df_merge$date)

#--load the phenophases data
phenos.path<-"D:/data/photocold_project/event_length/Using_sites_in_Fluxnet2015/"
load(paste0(phenos.path,"df_events_length.RDA"))
#Over_days_length-->days when gpp overestimated 
df_phenos<-df_events_all[,c("sitename","Year","sos","peak","Over_days_length")]
names(df_phenos)<-c("sitename","year","sos","peak","Over_days_length")

#------merge flux data and phenos data------------------------
df_merge<-left_join(df_merge,df_phenos,by=c("sitename","year"))
#only keep the site-year when sos and peak both available
df_final<-df_merge[!is.na(df_merge$sos)&!is.na(df_merge$peak),]
#calculate the PAR =fapar_itpl *ppdf_fluxnet2015
df_final$APAR<-df_final$fapar_itpl*df_final$ppfd_fluxnet2015

#-------------------------------------
#(2)simulate the constant LUE (based on DOY) for each site
#-------------------------------------
###---------
#First part: for all the available sites
###---------
df_final_all<-df_final
sites<-unique(df_final_all$sitename)
LUE_allsites<-c()
for(i in 1:length(sites)){
  df<-df_final_all %>%
    filter(sitename==sites[i])
  #select the data outside of the green-up period to model the LUE
  #df_out_greenup<-df[df$doy<df$sos|df$doy>df$peak,]
  # set the doy>peak in order to capture GPP peak using constant LUE
  df_out_greenup<-df[df$doy>df$peak,]  
  # #using robust regression:-->some error hence change to simple linear regression
  # lm_rob<-lmrob(gpp_obs~APAR - 1,data=df_out_greenup)  #set the intercept =0
  lm_norm<-lm(gpp_obs~APAR - 1,data=df_out_greenup)  #set the intercept =0
  stats_sum<-summary(lm_norm)
  # cbind(coef(lm_rob),confint(lm_rob, level = 0.95))
  LUE_temp<-data.frame("sitename"=sites[i],"LUE"=stats_sum$coefficients[1],"p-value"=stats_sum$coefficients[4])
  LUE_allsites<-rbind(LUE_allsites,LUE_temp)
}
LUE_final<-LUE_allsites[,c("sitename","LUE")]
names(LUE_final)<-c("sitename","conLUE")
#-----merge to df_final dataset----------
df_final_all<-left_join(df_final_all,LUE_final,by="sitename")

#-----------------------------------------------------
#compare the GPP_obs with GPP model through conLUE
#-----------------------------------------------------
df_final_all$gpp_conLUE<-df_final_all$APAR*df_final_all$conLUE
##save the data:
GPP_conLUE<-df_final_all
#
GPP_conLUE<-GPP_conLUE[,c("sitename","date","doy","sos","peak","Over_days_length","APAR","gpp_obs","gpp_mod_FULL","gpp_conLUE")]
# save the constant LUE GPP:
# save.path<-"D:/data/photocold_project/GPP_from_diffSources/site_scale/"
# save(GPP_conLUE,file = paste0(save.path,"GPP_conLUE_FLUX2015_daily.RDA"))
#first have a look at the performance for each site or all sites
#1) for each site:
df_final_all %>%
  group_by(sitename,doy)%>%
  summarise(obs = mean(gpp_obs, na.rm = TRUE),
            obs_sd = sd(gpp_obs,na.rm = T),
            mod_conLUE = mean(gpp_conLUE,na.rm=T),
            mod_conLUE_sd = sd(gpp_conLUE,na.rm = T),
            mod_Pmodel = mean(gpp_mod_FULL,na.rm=T),
            mod_Pmodel_sd = sd(gpp_mod_FULL,na.rm = T)
            ) %>%
    pivot_longer(c(obs,mod_conLUE,mod_Pmodel),
               names_to = "Source", values_to = "gpp") %>%
  ggplot(aes(doy, gpp, color = Source))+
  geom_line()+
  scale_color_manual(values = c("mod_conLUE" = "orange",
   "mod_Pmodel"="steelblue2","obs" = "black"),
   labels = c("constant LUE model","P model","Obs.")) +
  labs(y = expression( paste("GPP (g C m"^-2, " d"^-1, ")" ) ),
       x = "Day of year")+
  facet_wrap(~sitename)

#2) for all the sites:
df_gpp_all<-df_final_all %>%
  group_by(doy)%>%
  summarise(obs = mean(gpp_obs, na.rm = TRUE),
            obs_sd = sd(gpp_obs,na.rm = T),
            mod_conLUE = mean(gpp_conLUE,na.rm=T),
            mod_conLUE_sd = sd(gpp_conLUE,na.rm = T),
            mod_Pmodel = mean(gpp_mod_FULL,na.rm=T),
            mod_Pmodel_sd = sd(gpp_mod_FULL,na.rm = T)) %>%
  mutate(obs_mean=obs,obs_m_sd=obs - obs_sd,obs_a_sd=obs + obs_sd,
         mod_conLUE_mean=mod_conLUE,
         mod_conLUE_m_sd=mod_conLUE - mod_conLUE_sd,
         mod_conLUE_a_sd=mod_conLUE + mod_conLUE_sd,
         mod_Pmodel_mean=mod_Pmodel,
         mod_Pmodel_m_sd=mod_Pmodel - mod_Pmodel_sd,
         mod_Pmodel_a_sd=mod_Pmodel + mod_Pmodel_sd) %>%
  mutate(obs=NULL,mod_conLUE=NULL,mod_Pmodel=NULL) %>%   ##remove the variables do not use
  pivot_longer(c(obs_mean,mod_conLUE_mean,mod_Pmodel_mean),
               names_to = "Source", values_to = "gpp") 
df_gpp_all %>%
  ggplot(aes(doy, gpp, color = Source))+
  # geom_ribbon(aes(ymin = mod_conLUE_m_sd, ymax = mod_conLUE_a_sd), fill="orange",alpha = 0.5,color=NA)+
  # geom_ribbon(aes(ymin = mod_Pmodel_m_sd, ymax = mod_Pmodel_a_sd), fill="skyblue",alpha = 0.5,color=NA)+
  # geom_ribbon(aes(ymin = obs_m_sd, ymax = obs_a_sd),fill="green3", alpha = 0.3,color=NA)+
  geom_line()+
  scale_color_manual(values = c("mod_conLUE_mean" = "orange",
                                "mod_Pmodel_mean"="steelblue2","obs_mean" = "black"),
                     labels = c("constant LUE model","P model","Obs.")) +
  # scale_fill_manual(values = c("mod_conLUE_sd" = "orange",
  #                               "mod_Pmodel_sd"="steelblue2","obs_sd" = "gray"),
  #                    labels = c("constant LUE model","P model","Obs.")) +
  labs(y = expression( paste("GPP (g C m"^-2, " d"^-1, ")" ) ),
       x = "Day of year")+
  theme_classic()


###---------
#Second part: for the sites that identify as the "GPP overestimated"
###---------
df_final_over<-df_final[!is.na(df_final$Over_days_length)&df_final$Over_days_length>20,]
sites<-unique(df_final_over$sitename)
LUE_allsites_new<-c()
for(i in 1:length(sites)){
  df<-df_final_over %>%
    filter(sitename==sites[i])
  #select the data outside of the green-up period to model the LUE
  df_out_greenup<-df[df$doy<df$sos|df$doy>df$peak,]
  Pgpp_perc<-length(df_out_greenup$gpp_obs[!is.na(df_out_greenup$gpp_obs)&df_out_greenup$gpp_obs>1])/length(df_out_greenup$gpp_obs[!is.na(df_out_greenup$gpp_obs)])
  # plot(df_out_greenup$APAR,df_out_greenup$gpp_obs)
  #using normal(simple) regression:
  # lm_rob<-lmrob(gpp_obs~APAR - 1,data=df_out_greenup)  #set the intercept =0
  lm_norm<-lm(gpp_obs~APAR - 1,data=df_out_greenup)  #set the intercept =0
  stats_sum<-summary(lm_norm)
  # cbind(coef(lm_rob),confint(lm_rob, level = 0.95))
  LUE_temp<-data.frame("sitename"=sites[i],"LUE"=stats_sum$coefficients[1],"p-value"=stats_sum$coefficients[4])
  # abline(0,LUE_temp$LUE,col="blue")
  LUE_allsites_new<-rbind(LUE_allsites_new,LUE_temp)
}
LUE_final<-LUE_allsites_new[,c("sitename","LUE")]
names(LUE_final)<-c("sitename","conLUE")
LUE_final<-LUE_final[!is.na(LUE_final$sitename),]
#-----merge to df_final dataset----------
df_final_over<-left_join(df_final_over,LUE_final,by="sitename")

#-----------------------------------------------------
#compare the GPP_obs with GPP model through conLUE
#-----------------------------------------------------
df_final_over$gpp_conLUE<-df_final_over$APAR*df_final_over$conLUE
#first have a look at the performance for each site or all sites
#1) for each site:
df_final_over %>%
  group_by(sitename,doy)%>%
  summarise(obs = mean(gpp_obs, na.rm = TRUE),
            obs_sd = sd(gpp_obs,na.rm = T),
            mod_conLUE = mean(gpp_conLUE,na.rm=T),
            mod_conLUE_sd = sd(gpp_conLUE,na.rm = T),
            mod_Pmodel = mean(gpp_mod_FULL,na.rm=T),
            mod_Pmodel_sd = sd(gpp_mod_FULL,na.rm = T)
  ) %>%
  pivot_longer(c(obs,mod_conLUE,mod_Pmodel),
               names_to = "Source", values_to = "gpp") %>%
  ggplot(aes(doy, gpp, color = Source))+
  geom_line()+
  scale_color_manual(values = c("mod_conLUE" = "orange",
                                "mod_Pmodel"="steelblue2","obs" = "black"),
                     labels = c("constant LUE model","P model","Obs.")) +
  labs(y = expression( paste("GPP (g C m"^-2, " d"^-1, ")" ) ),
       x = "Day of year")+
  facet_wrap(~sitename)

#2) for all the sites:
df_gpp_over<-df_final_over %>%
  group_by(doy)%>%
  summarise(obs = mean(gpp_obs, na.rm = TRUE),
            obs_sd = sd(gpp_obs,na.rm = T),
            mod_conLUE = mean(gpp_conLUE,na.rm=T),
            mod_conLUE_sd = sd(gpp_conLUE,na.rm = T),
            mod_Pmodel = mean(gpp_mod_FULL,na.rm=T),
            mod_Pmodel_sd = sd(gpp_mod_FULL,na.rm = T)) %>%
  mutate(obs_mean=obs,obs_m_sd=obs - obs_sd,obs_a_sd=obs + obs_sd,
         mod_conLUE_mean=mod_conLUE,
         mod_conLUE_m_sd=mod_conLUE - mod_conLUE_sd,
         mod_conLUE_a_sd=mod_conLUE + mod_conLUE_sd,
         mod_Pmodel_mean=mod_Pmodel,
         mod_Pmodel_m_sd=mod_Pmodel - mod_Pmodel_sd,
         mod_Pmodel_a_sd=mod_Pmodel + mod_Pmodel_sd) %>%
  mutate(obs=NULL,mod_conLUE=NULL,mod_Pmodel=NULL) %>%   ##remove the variables do not use
  pivot_longer(c(obs_mean,mod_conLUE_mean,mod_Pmodel_mean),
               names_to = "Source", values_to = "gpp") 
df_gpp_over %>%
  ggplot(aes(doy, gpp, color = Source))+
  # geom_ribbon(aes(ymin = mod_conLUE_m_sd, ymax = mod_conLUE_a_sd), fill="orange",alpha = 0.5,color=NA)+
  # geom_ribbon(aes(ymin = mod_Pmodel_m_sd, ymax = mod_Pmodel_a_sd), fill="skyblue",alpha = 0.5,color=NA)+
  # geom_ribbon(aes(ymin = obs_m_sd, ymax = obs_a_sd),fill="green3", alpha = 0.3,color=NA)+
  geom_line()+
  scale_color_manual(values = c("mod_conLUE_mean" = "orange",
                                "mod_Pmodel_mean"="steelblue2","obs_mean" = "black"),
                     labels = c("constant LUE model","P model","Obs.")) +
  # scale_fill_manual(values = c("mod_conLUE_sd" = "orange",
  #                               "mod_Pmodel_sd"="steelblue2","obs_sd" = "gray"),
  #                    labels = c("constant LUE model","P model","Obs.")) +
  labs(y = expression( paste("GPP (g C m"^-2, " d"^-1, ")" ) ),
       x = "Day of year")+
  theme_classic()

###------------------------
#Third part: check the failed LUE modelled site: US-NR1, US-WCr
##-->because the robust regression failed to fit the right linear regression-->
##hence now change to simple linear regression
###------------------------
# df_final_over %>%
#   filter(sitename=="US-NR1")%>%
#   ggplot()+
#   geom_point(aes(x=date,y=gpp_obs))
