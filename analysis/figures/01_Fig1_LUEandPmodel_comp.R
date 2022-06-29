#############################################
#Aim:comparing the simple LUE, P-model and EC based GPP:
#-->to check if only consider the instanenous environmental vars will result in 
#overestimation of GPP:
##-->need to update the code tomorrrow!
#############################################
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(lme4)
library(tidyverse)
# remotes::install_github("computationales/ingestr") #install the package
library(ingestr)
library(rbeni)

#-----------------
#(1) tidy up the data
#load the merged fluxes data and also used the 
#phenophases(SOS and POS) extracted to determine the green-up period  
#-----------------
#--load the merged fluxes data
load.path<-"./data-raw/raw_data/Merged_data/"
load(paste0(load.path,"Merged_Flux_and_VIs.RDA"))
#
df_merge$year<-year(df_merge$date)
df_merge$doy<-yday(df_merge$date)

#--load the phenophases data
phenos.path<-"./data/event_length/"
load(paste0(phenos.path,"df_events_length.RDA"))
#Over_days_length-->days when gpp overestimated 
df_phenos<-df_events_all[,c("sitename","Year","sos","peak","Over_days_length")]
names(df_phenos)<-c("sitename","year","sos","peak","Over_days_length")

#------merge flux data and phenos data------------------------
df_merge<-left_join(df_merge,df_phenos,by=c("sitename","year"))
#only keep the site-year when sos and peak both available
df_final<-df_merge[!is.na(df_merge$sos)&!is.na(df_merge$peak),]

#--------------------------
#(2)start to fit constant LUE model and simple LUE model using lme()
#--------------------------
ddf <- df_final %>% 
  as_tibble() %>%
  # mutate(doy = lubridate::yday(date)) %>% 
  mutate(greenup = ifelse(doy > sos & doy < peak, TRUE, FALSE))
#check the data avaiablilty:
# visdat::vis_miss(ddf, warn_large_data = FALSE)

## Empirical LUE models
# Determine bias of early season bias of empirical LUE models 
# by fitting models outside greenup period.-->
# two LUE model: constant LUE and LUE with VPD and temp as the drivers
ddf <- ddf %>% 
  #the variables of ppfd_fluxnet2015 has some probelm==> using PPFD_IN_fullday_mean_fluxnet2015
  mutate(lue = gpp_obs / (fapar_itpl * PPFD_IN_fullday_mean_fluxnet2015)) %>% 
  mutate(lue = remove_outliers(lue)) #using the functions in the ingestr

#check the data for each site:
ddf %>%
  group_by(sitename)%>%
  ggplot()+
    # geom_point(aes(x=date,y=lue))+
    ##after the check, we found the PPFD_IN data is not avaiable in CN-Qia-->need to update
    ##update now-->May,3rd
    geom_point(aes(x=date,y=PPFD_IN_fullday_mean_fluxnet2015))+ 
    facet_wrap( ~sitename)

###2a).LUE model
##take mean LUE for constant-LUE model
mod_constlue <- ddf %>% 
  filter(!greenup) %>% 
  pull(lue) %>% 
  mean(., na.rm = TRUE)

## LUE as a linear function of temp and vpd
mod_lue_temp_vpd <- lm(lue ~ temp_day_fluxnet2015 + vpd_day_fluxnet2015, 
                       data = ddf %>% 
                         filter(!greenup))
##2b.) add year and elevation mixed effects model
ddf <- ingestr::siteinfo_fluxnet2015 %>%
  select(sitename, elv) %>%
  right_join(ddf, by = "sitename") %>%
  mutate(year = lubridate::year(date))

##remove the na values and infinite data
tmp <- ddf %>% 
  dplyr::filter(!greenup) %>%
  dplyr::filter(PPFD_IN_fullday_mean_fluxnet2015 > 5) %>% 
  dplyr::select(temp_day_fluxnet2015, vpd_day_fluxnet2015, lue, sitename, year) %>%
  drop_na() %>% 
  dplyr::filter(!is.infinite(lue) & !is.nan(lue)) %>% 
  dplyr::filter(!is.infinite(temp_day_fluxnet2015) & !is.nan(temp_day_fluxnet2015)) %>% 
  dplyr::filter(!is.infinite(vpd_day_fluxnet2015) & !is.nan(vpd_day_fluxnet2015)) %>% 
  filter(vpd_day_fluxnet2015 > 0 & lue > 0) %>% 
  droplevels()

mod_lmer <- lmer(lue ~ temp_day_fluxnet2015+log(vpd_day_fluxnet2015) + (1|sitename),
                   data=tmp)
summary(mod_lmer)
#2c).merge the data
ddf <- ddf %>% 
  # filter(sitename!="CN-Qia")%>%  ## since mod_lmer do not include the CN-Qia site
  #==>glmer prediction has some probelm
  mutate(lue_temp_vpd = predict(mod_lue_temp_vpd,newdata = .))%>%
  mutate(gpp_temp_vpd = lue_temp_vpd * fapar_itpl * PPFD_IN_fullday_mean_fluxnet2015) %>%
  mutate(lue_lmer = predict(mod_lmer, newdata = .)) %>%
  mutate(gpp_lmer = lue_lmer * fapar_itpl * PPFD_IN_fullday_mean_fluxnet2015) %>%
  mutate(gpp_lue_const = mod_constlue * fapar_itpl * PPFD_IN_fullday_mean_fluxnet2015)

#--------------------------
#(3)compare the models
#--------------------------
#change the names of modelled GPP:
ddf<-ddf %>%
  mutate(gpp_pmodel=gpp_mod_FULL,
         gpp_mod_FULL=NULL)
# ## P-model
# ddf %>% 
#   analyse_modobs2("gpp_pmodel", "gpp_obs", type = "hex")
# 
# ## constant LUE model
# ddf %>% 
#   analyse_modobs2("gpp_lue_const", "gpp_obs", type = "hex")
# 
# ## LUE ~ temp + VPD model:lm
ddf %>%
  analyse_modobs2("gpp_temp_vpd", "gpp_obs", type = "hex")
# 
# ## LUE ~ temp + VPD model:glmer
ddf %>%
  filter(!is.nan(gpp_lmer) & !is.infinite(gpp_lmer))%>%
  analyse_modobs2("gpp_lmer", "gpp_obs", type = "hex")

#--------------------------
#(4)Mean seasonal cycle
#--------------------------
df_meandoy <- ddf %>% 
  group_by(sitename, doy) %>% 
  dplyr::summarise(across(starts_with("gpp_"), mean, na.rm = TRUE))

##plot by site:
df_meandoy %>% 
  filter(!is.nan(gpp_lmer) & !is.infinite(gpp_lmer))%>% ##this filter is important
  pivot_longer(c(gpp_obs, gpp_pmodel, gpp_lue_const,gpp_temp_vpd,gpp_lmer), names_to = "model", values_to = "gpp") %>% 
  #fct_relevel: in tidyverse package
  mutate(model = fct_relevel(model, "gpp_obs", "gpp_pmodel", "gpp_lue_const","gpp_temp_vpd","gpp_lmer")) %>% 
  dplyr::filter((model %in% c( "gpp_obs", "gpp_pmodel","gpp_temp_vpd","gpp_lmer"))) %>% 
  # filter(sitename=="DE-Tha")%>%
  ggplot() +
  # geom_ribbon(
  #   aes(x = doy, ymin = obs_min, ymax = obs_max), 
  #   fill = "black", 
  #   alpha = 0.2
  #   ) +
  geom_line(aes(x = doy, y = gpp, color = model), size = 0.4) +
  labs(y = expression( paste("Simulated GPP (g C m"^-2, " d"^-1, ")" ) ), 
       x = "DOY") +
  facet_wrap( ~sitename, ncol = 3 ) +    # , labeller = labeller(climatezone = list_rosetta)
  theme_gray() +
  theme(legend.position = "bottom") +
  scale_color_manual(
    name="Model: ",
    values=c("black", "red", "royalblue", "darkgoldenrod", "springgreen", "orchid4")
  )

# ggsave("./manuscript/test_files/gpp_meandoy.pdf", height = 25, width = 8)

##############################
## Normalise to peak season
##############################
norm_to_peak <- function(df, mod, obs){
  # df<-ddf_t
  # mod<-"gpp_lmer"
  # obs<-"gpp_obs"
  
  q75_obs <- quantile(df[[obs]], probs = 0.75, na.rm = TRUE)
  q75_mod <- quantile(df[[mod]], probs = 0.75, na.rm = TRUE)
  
  ## normalise mod
  #add by YP:first need to change infinite values to NA:
  # df[[mod]][is.nan(df[[mod]])]<-NA
  df[[mod]][is.infinite(df[[mod]])]<-NA
  df[[mod]] <- df[[mod]] * 
    mean(df[[obs]][df[[obs]]>q75_obs], na.rm = TRUE) / 
    mean(df[[mod]][df[[mod]]>q75_mod], na.rm = TRUE)  #YP revised here:some error from Beni's functions
  
  return(df)
}

ddf_norm <- ddf %>% 
  group_by(sitename) %>% 
  nest() %>% 
  mutate(data = purrr::map(data, ~norm_to_peak(., "gpp_pmodel", "gpp_obs"))) %>% 
  # mutate(data = purrr::map(data, ~norm_to_peak(., "gpp_bess", "gpp_obs"))) %>% 
  # mutate(data = purrr::map(data, ~norm_to_peak(., "gpp_rf", "gpp_obs"))) %>% 
  # mutate(data = purrr::map(data, ~norm_to_peak(., "gpp_bf", "gpp_obs"))) %>% 
  mutate(data = purrr::map(data, ~norm_to_peak(., "gpp_lue_const", "gpp_obs"))) %>% 
  mutate(data = purrr::map(data, ~norm_to_peak(., "gpp_temp_vpd", "gpp_obs"))) %>% 
  mutate(data = purrr::map(data, ~norm_to_peak(., "gpp_lmer", "gpp_obs"))) %>% 
  unnest(data)

### Plot normalised by site
df_meandoy_norm <- ddf_norm %>% 
  group_by(sitename, doy) %>% 
  dplyr::summarise(across(starts_with("gpp_"), mean, na.rm = TRUE))

plot_sites<-df_meandoy_norm %>% 
  filter(!is.nan(gpp_lmer) & !is.infinite(gpp_lmer))%>% ##this filter is important
  pivot_longer(c(gpp_obs, gpp_pmodel, gpp_lue_const, gpp_temp_vpd,gpp_lmer), names_to = "model", values_to = "gpp") %>%
  mutate(model = fct_relevel(model, "gpp_obs", "gpp_pmodel", "gpp_lue_const", "gpp_temp_vpd","gpp_lmer")) %>%
  dplyr::filter((model %in% c( "gpp_obs", "gpp_pmodel","gpp_lmer"))) %>%  ##only select one model
  # pivot_longer(c(gpp_obs, gpp_pmodel, gpp_lue_const, gpp_temp_vpd), names_to = "model", values_to = "gpp") %>% 
  # mutate(model = fct_relevel(model, "gpp_obs", "gpp_pmodel", "gpp_lue_const", "gpp_temp_vpd")) %>% 
  # dplyr::filter((model %in% c( "gpp_obs", "gpp_pmodel", "gpp_temp_vpd"))) %>% 
  # filter(sitename=="US-UMB")%>%
  ggplot() +
  # geom_ribbon(
  #   aes(x = doy, ymin = obs_min, ymax = obs_max), 
  #   fill = "black", 
  #   alpha = 0.2
  #   ) +
  geom_line(aes(x = doy, y = gpp, color = model)) +
  labs(y = expression( paste("GPP (g C m"^-2, " d"^-1, ")" ) ),
       x = "DoY") +
  facet_wrap( ~sitename) +
  # theme_gray() +
  scale_color_manual("GPP sources",values = c("gpp_obs" = "black",
    "gpp_pmodel" = "red","gpp_lmer"="dodgerblue"),
    labels = c("Obervations","P-model","LUE-lme"))+
  # scale_color_manual(
  #   name="Model: ",
  #   values=c("black", "red", "royalblue", "darkgoldenrod", "springgreen", "orchid4")
  # )
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

# ggsave("./manuscript/test_files/gpp_meandoy_norm.pdf", height = 25, width = 8)
ggsave("./manuscript/figures/FigS_eachsite_gpp_meandoy_norm.png",width = 20,height = 20)

####################
#(5)prepare the official plots:
#using norm_GPP results
####################
#----------------Figure :plotting for 
#load the site infos:
load(paste0("./data-raw/raw_data/sites_info/","Pre_selected_sites_info.RDA"))
sites.info<-df_sites_sel;rm(df_sites_sel)
#
ddf_norm<-left_join(ddf_norm,sites.info)
ddf_norm<-ddf_norm %>%
  mutate(Clim_PFTs=paste0(koeppen_code,"-",classid))

#check sites in different Clim.-PFTs
sites.info%>%
  mutate(Clim_PFTs=paste0(koeppen_code,"-",classid)) %>%
  select(sitename,Clim_PFTs)%>%
  # filter(Clim_PFTs=="Cfa-ENF")
  filter(Clim_PFTs=="Dfb-ENF")

### I. Plot normalised by Clim_PFTs
#test why the norm results is strange for sites:US-NR1, US-PFa, and JP-MBF
#-->the probelm is because of the normalization has probelm:
#original data:
# ddf_t<-ddf%>%
#   filter(sitename=="US-NR1")%>%
#   select(date,gpp_obs,gpp_pmodel,gpp_lue_const,gpp_temp_vpd,gpp_lmer)
#   # pivot_longer(starts_with("gpp_"),names_to = "model",values_to = "gpp")%>%
#   # ggplot()+
#   # geom_line(aes(x=date,y=gpp,col=model))
# ddf_t%>%
#   filter(!is.nan(gpp_lmer) & !is.infinite(gpp_lmer))%>% #important syntax
#   ggplot()+
#   geom_point(aes(x=date,y=gpp_lmer,col="lmer(temp,vpd)"))+
#   geom_point(aes(x=date,y=gpp_temp_vpd,col="lm(temp,vpd)"))
# #norm data
# ddf_norm_t<-ddf_norm%>%
#   filter(sitename=="US-NR1")%>%
#   select(date,gpp_obs,gpp_pmodel,gpp_lue_const,gpp_temp_vpd,gpp_lmer)
# ddf_norm_t%>%
#   filter(!is.nan(gpp_lmer) & !is.infinite(gpp_lmer))%>% #important syntax
#   ggplot()+
#   geom_point(aes(x=date,y=gpp_lmer,col="lmer(temp,vpd)"))+
#   geom_point(aes(x=date,y=gpp_temp_vpd,col="lm(temp,vpd)"))

df_meandoy_norm_Clim_PFTs <- ddf_norm %>%
  filter(sitename!="JP-MBF" & sitename!="US-NR1" & sitename!="US-PFa")%>% ##remove the sites do not performed when lmer
  group_by(Clim_PFTs, doy) %>% 
  dplyr::summarise(across(starts_with("gpp_"), mean, na.rm = TRUE))

#Figure for Clim.-PFTs
plot_final<-df_meandoy_norm_Clim_PFTs %>% 
  filter(!is.nan(gpp_lmer) & !is.infinite(gpp_lmer))%>% ##this filter is important
  pivot_longer(c(gpp_obs, gpp_pmodel, gpp_lue_const, gpp_temp_vpd,gpp_lmer), names_to = "model", values_to = "gpp") %>%
  mutate(model = fct_relevel(model, "gpp_obs", "gpp_pmodel", "gpp_lue_const", "gpp_temp_vpd","gpp_lmer")) %>%
  dplyr::filter((model %in% c( "gpp_obs", "gpp_pmodel","gpp_lmer"))) %>% ##select only one lm
  ggplot() +
  geom_line(aes(x = doy, y = gpp, color = model)) +
  labs(y = expression( paste("GPP (g C m"^-2, " d"^-1, ")" ) ),
       x = "DoY") +
  # facet_wrap( ~Clim_PFTs, ncol = 3, scales = "free_y" ) +
  facet_wrap( ~Clim_PFTs) +
  # theme_gray() +
  # theme(legend.position = "bottom") +
  # scale_color_manual(
  #   name="Model: ",
  #   values=c("black", "red", "royalblue", "darkgoldenrod", "springgreen", "orchid4"))+
  scale_color_manual("GPP sources",values = c("gpp_obs" = "black",
     "gpp_pmodel" = "red", "gpp_lmer" = "dodgerblue"),
     labels = c("Observations","P-model","LUE-lme")) +
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
ggsave("./manuscript/figures/Figure1_gpp_meandoy_norm_forClimPFTs.png",plot_final,width = 15,height = 10)
