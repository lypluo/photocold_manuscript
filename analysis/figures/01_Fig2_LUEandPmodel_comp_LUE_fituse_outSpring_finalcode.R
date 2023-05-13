#############################################
#Aim:comparing the simple LUE, P-model and EC based GPP:
#-->to check if only consider the instanenous environmental vars will result in 
#overestimation of GPP:
#############################################
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(lme4)
library(tidyverse)
# remotes::install_github("computationales/ingestr") #install the package
library(ingestr)
devtools::load_all("D:/Github/rbeni/")
# library(rbeni)
library(plotrix) #calculate the standard error 

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
##add the variable month
df_final$month<-month(df_final$date)
#--------------------------
#(2)start to fit constant LUE model and simple LUE model using lme()
#--------------------------
ddf <- df_final %>% 
  as_tibble() %>%
  # mutate(doy = lubridate::yday(date)) %>% 
  mutate(greenup = ifelse(doy > sos & doy < peak, TRUE, FALSE))%>%
  mutate(spring = ifelse(month >=3 & month<=5,TRUE,FALSE))
#check the data avaiablilty:
# visdat::vis_miss(ddf, warn_large_data = FALSE)

## Empirical LUE models
# Determine bias of early season bias of empirical LUE models 
# by fitting models outside spring (March-May).-->
# two LUE model: constant LUE and LUE with VPD and temp as the drivers
ddf <- ddf %>% 
  #the variables of ppfd_fluxnet2015 has some probelm==> using PPFD_IN_fullday_mean_fluxnet2015
  mutate(lue = gpp_obs / (fapar_itpl * PPFD_IN_fullday_mean_fluxnet2015)) %>% 
  mutate(lue = remove_outliers(lue)) #using the functions in the ingestr
#save the ddf-->using for further analysis in the future
# save(ddf,file = paste0("D:/Github/photocold_manuscript/test/test_dataset/","GPP_andLUE.RDA"))

#check the data for each site:
# ddf %>%
#   group_by(sitename)%>%
#   ggplot()+
#     # geom_point(aes(x=date,y=lue))+
#     ##after the check, we found the PPFD_IN data is not avaiable in CN-Qia-->need to update
#     ##update now-->May,3rd
#     geom_point(aes(x=date,y=PPFD_IN_fullday_mean_fluxnet2015))+ 
#     facet_wrap( ~sitename)

###2a).LUE model-->using the data outside of the green-up period
##take mean LUE for constant-LUE model
mod_constlue <- ddf %>% 
  filter(!spring) %>% 
  pull(lue) %>% 
  mean(., na.rm = TRUE)

## LUE as a linear function of temp and vpd
mod_lue_temp_vpd <- lm(lue ~ temp_day_fluxnet2015 + vpd_day_fluxnet2015, 
                       data = ddf %>% 
                         filter(!spring))
##2b.) add year and elevation mixed effects model
ddf <- ingestr::siteinfo_fluxnet2015 %>%
  select(sitename, elv) %>%
  right_join(ddf, by = "sitename") %>%
  mutate(year = lubridate::year(date))

##remove the na values and infinite data
tmp <- ddf %>% 
  dplyr::filter(!spring) %>%
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
# # ## P-model
# ddf %>%
#   analyse_modobs2("gpp_pmodel", "gpp_obs", type = "hex")
# # 
# # ## constant LUE model
# ddf %>%
#   analyse_modobs2("gpp_lue_const", "gpp_obs", type = "hex")
# # 
# # ## LUE ~ temp + VPD model:lm
# ddf %>%
#   analyse_modobs2("gpp_temp_vpd", "gpp_obs", type = "hex")
# # 
# # ## LUE ~ temp + VPD model:glmer
# ddf %>%
#   filter(!is.nan(gpp_lmer) & !is.infinite(gpp_lmer))%>%
#   analyse_modobs2("gpp_lmer", "gpp_obs", type = "hex")

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

# Residual analysis for reply to reviewer 1 ----------------------------
tmp <- ddf |> 
  mutate(res_pmodel = gpp_pmodel - gpp_obs,
         res_lmm = gpp_lmer - gpp_obs) |> 
  mutate(doy = lubridate::yday(date)) |> 
  mutate(spring = ifelse(doy < 365/2, TRUE, FALSE)) |> 
  mutate(fapar_bin = cut(fapar_itpl, breaks = seq(0, 1, by = 0.1)))

tmp_agg <- tmp |> 
  group_by(sitename, fapar_bin, spring) |> 
  summarise(res_pmodel = mean(res_pmodel, na.rm = TRUE)) |> 
  pivot_wider(names_from = c("spring"), values_from = "res_pmodel") |> 
  mutate(diff_spring = `TRUE` - `FALSE`) |> 
  ungroup() |> 
  group_by(sitename) |> 
  summarise(diff_spring = mean(diff_spring, na.rm = TRUE))

plot_1 <- ggplot(
  data = tmp_agg,
  aes(x = sitename, y = diff_spring)) +
  geom_bar(stat = "identity") +
  labs(y = expression( paste("Mean difference in spring vs rest bias within fAPAR bin (g C m"^-2, " d"^-1, ")" ) ),
       x = "") +
  theme_classic() +
  coord_flip()


## Residual vs fAPAR bin by site -----------------
### P-model bias ----------------
tmp |> 
  filter(sitename %in% unique(tmp$sitename)[1:18]) |> 
  ggplot(aes(x = fapar_bin, y = res_pmodel, fill = spring)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_fill_manual(
    name="Spring",
    values=c("#777055ff", "#29a274ff")
  ) +
  facet_wrap( ~sitename, ncol = 3 ) +
  labs() +
  labs(y = expression( paste("P-model GPP residual (g C m"^-2, " d"^-1, ")" ) ),
       x = "fAPAR bin") +
  theme_classic()
ggsave("./manuscript/figures/FigADD_residual_fapar_pmodel_1.png", width = 12, height = 18)

#### Example ----------
gg1 <- tmp |> 
  filter(sitename %in% c("US-UMd", "IT-Ren")) |> 
  ggplot(aes(x = fapar_bin, y = res_pmodel, fill = spring)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_fill_manual(
    name="Spring",
    values=c("#777055ff", "#29a274ff")
  ) +
  facet_wrap( ~sitename, ncol = 3 ) +
  labs() +
  labs(y = expression( paste("P-model GPP bias (g C m"^-2, " d"^-1, ")" ) ),
       x = "fAPAR bin") +
  theme_classic()

gg2 <- tmp |> 
  filter(sitename %in% c("US-UMd", "IT-Ren")) |> 
  ggplot(aes(x = fapar_bin, y = res_lmm, fill = spring)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_fill_manual(
    name="Spring",
    values=c("#777055ff", "#29a274ff")
  ) +
  facet_wrap( ~sitename, ncol = 3 ) +
  labs() +
  labs(y = expression( paste("LMM GPP bias (g C m"^-2, " d"^-1, ")" ) ),
       x = "fAPAR bin") +
  theme_classic()

cowplot::plot_grid(gg1, gg2, ncol = 1, labels = c("a", "b"))
ggsave("./manuscript/figures/FigADD_residual_fapar_EXAMPLE.pdf", width = 12, height = 10)
# ggsave("./manuscript/figures/FigADD_residual_fapar_EXAMPLE.png", width = 12, height = 10)


tmp |> 
  filter(sitename %in% unique(tmp$sitename)[19:length(unique(tmp$sitename))]) |> 
  ggplot(aes(x = fapar_bin, y = res_pmodel, fill = spring)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_fill_manual(
    name="Spring",
    values=c("#777055ff", "#29a274ff")
  ) +
  facet_wrap( ~sitename, ncol = 3 ) +
  labs() +
  labs(y = expression( paste("P-model GPP residual (g C m"^-2, " d"^-1, ")" ) ),
       x = "fAPAR bin") +
  theme_classic()
ggsave("./manuscript/figures/FigADD_residual_fapar_pmodel_2.png", width = 12, height = 18)

### LMM model bias ------------------
tmp |> 
  filter(sitename %in% unique(tmp$sitename)[1:18]) |> 
  ggplot(aes(x = fapar_bin, y = res_lmm, fill = spring)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_fill_manual(
    name="Spring",
    values=c("#777055ff", "#29a274ff")
  ) +
  facet_wrap( ~sitename, ncol = 3 ) +
  labs() +
  labs(y = expression( paste("P-model GPP residual (g C m"^-2, " d"^-1, ")" ) ),
       x = "fAPAR bin") +
  theme_classic()
ggsave("./manuscript/figures/FigADD_residual_fapar_lmm_1.png", width = 12, height = 18)

tmp |> 
  filter(sitename %in% unique(tmp$sitename)[19:length(unique(tmp$sitename))]) |> 
  ggplot(aes(x = fapar_bin, y = res_lmm, fill = spring)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_fill_manual(
    name="Spring",
    values=c("#777055ff", "#29a274ff")
  ) +
  facet_wrap( ~sitename, ncol = 3 ) +
  labs() +
  labs(y = expression( paste("P-model GPP residual (g C m"^-2, " d"^-1, ")" ) ),
       x = "fAPAR bin") +
  theme_classic()
ggsave("./manuscript/figures/FigADD_residual_fapar_lmm_2.png", width = 12, height = 18)


## Difference between spring and autumn bias across sites
tmp |> 
  ggplot(aes(x = sitename, y = res_pmodel, fill = spring)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(
    name="Spring",
    values=c("#777055ff", "#29a274ff")
  ) +
  coord_flip()
  

# Normalise to peak season ----------------------------
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
  mutate(lue_pmodel=gpp_pmodel/c(PPFD_IN_fullday_mean_fluxnet2015*fapar_itpl))%>%
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

ddf_lue<-ddf %>% 
  mutate(lue_pmodel=gpp_pmodel/c(PPFD_IN_fullday_mean_fluxnet2015*fapar_itpl))%>%
  #update in May, 2023-->remove the inf value in LUE
  filter(is.finite(lue_pmodel) & is.finite(lue_lmer))%>%
  filter(lue_pmodel>-0.1 & lue_pmodel<0.08)%>%
  filter(lue_lmer>-0.1 & lue_lmer<0.08)

ddf_norm_lue <- ddf_lue%>%
  # filter(lue_pmodel>0 & lue_lmer>0)%>%
  group_by(sitename) %>% 
  nest() %>% 
  #add in May, 2023:
  mutate(data = purrr::map(data, ~norm_to_peak(., "lue_pmodel", "lue"))) %>% 
  mutate(data = purrr::map(data, ~norm_to_peak(., "lue_temp_vpd", "lue"))) %>% 
  mutate(data = purrr::map(data, ~norm_to_peak(., "lue_lmer", "lue"))) %>% 
  unnest(data)

### Plot normalised by site
df_meandoy<-ddf %>%
  group_by(sitename, doy) %>% 
  dplyr::summarise(across(starts_with("gpp_"), mean, na.rm = TRUE))
df_meandoy_norm <- ddf_norm %>% 
  group_by(sitename, doy) %>% 
  dplyr::summarise(across(starts_with("gpp_"), mean, na.rm = TRUE))

plot_sites_norm <- df_meandoy_norm %>% 
  filter(!is.nan(gpp_lmer) & !is.infinite(gpp_lmer))%>% ##this filter is important
  pivot_longer(c(gpp_obs, gpp_pmodel, gpp_lue_const, gpp_temp_vpd, gpp_lmer), names_to = "model", values_to = "gpp") %>%
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
  geom_line(aes(x = doy, y = gpp, color = model),size=0.8) +
  labs(y = expression( paste("GPP (g C m"^-2, " d"^-1, ")" ) ),
       x = "DoY") +
  facet_wrap( ~sitename) +
  # theme_gray() +
  scale_color_manual("GPP sources",values = c("gpp_obs" = "black",
    "gpp_pmodel" = "orange","gpp_lmer"=adjustcolor("brown2",0.8)),
    labels = c(expression(GPP[obs]),expression(GPP[Pmodel]),expression(GPP[LME])))+
    # labels = c("Obervations","P-model","LME"))+
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
  )+
theme(legend.text.align = 0)  #align the legend (all the letter start at the same positoin)
# ggsave("./manuscript/test_files/gpp_meandoy_norm.pdf", height = 25, width = 8)
ggsave("./manuscript/figures/FigS_eachsite_gpp_meandoy_norm.png",width = 20,height = 20)

####################
# Publication figure ------------------------------
# (5) prepare the official plots:
# using norm_GPP results
####################
#----------------Figure :plotting for 
#load the site infos:
load(paste0("./data-raw/raw_data/sites_info/","Pre_selected_sites_info.RDA"))
sites.info<-df_sites_sel;rm(df_sites_sel)
#normalization:
ddf_norm<-left_join(ddf_norm,sites.info)
ddf_norm<-ddf_norm %>%
  mutate(Clim_PFTs=paste0(koeppen_code,"-",classid))
#add in May, 2023:
ddf_lue<-left_join(ddf_lue,sites.info)
ddf_lue<-ddf_lue %>%
  mutate(Clim_PFTs=paste0(koeppen_code,"-",classid))
ddf_norm_lue<-left_join(ddf_norm_lue,sites.info)
ddf_norm_lue<-ddf_norm_lue %>%
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
  dplyr::summarise(
    gpp_obs_min=min(gpp_obs,na.rm = T),
    gpp_obs_max=max(gpp_obs,na.rm = T),
    gpp_obs_sd=sd(gpp_obs,na.rm = T),
    gpp_obs_mean=mean(gpp_obs,na.rm=T),
    lue_obs_sd=sd(lue,na.rm = T),
    lue_obs_mean=mean(lue,na.rm=T),
    #add sos_mean and eos_mean
    sos_mean=round(mean(sos,na.rm=T),0),
    peak_min=round(min(peak,na.rm=T),0),
    fapar = mean(fapar_itpl, na.rm = TRUE),
    across(starts_with("gpp_"), mean, na.rm = TRUE),
    across(starts_with("lue"),mean,na.rm=TRUE)
    ) %>%
  #add in May,2023
  mutate(avg_period = ifelse(doy >= sos_mean & doy <= peak_min, TRUE, FALSE))

#---->for LUE:
df_meandoy_Clim_PFTs_lue <- ddf_lue %>%
  filter(sitename!="JP-MBF" & sitename!="US-NR1" & sitename!="US-PFa")%>% ##remove the sites do not performed when lmer
  #update in May, 2023:adding the LUE
  group_by(Clim_PFTs, doy) %>% 
  dplyr::summarise(
    lue_obs_sd=sd(lue,na.rm = T),
    lue_obs_mean=mean(lue,na.rm=T),
    #add sos_mean and eos_mean
    sos_mean=round(mean(sos,na.rm=T),0),
    peak_min=round(min(peak,na.rm=T),0),
    fapar = mean(fapar_itpl, na.rm = TRUE),
    across(starts_with("lue"),mean,na.rm=TRUE)
  ) %>%
  #add in May,2023
  mutate(avg_period = ifelse(doy >= sos_mean & doy <= peak_min, TRUE, FALSE))


df_meandoy_norm_Clim_PFTs_lue <- ddf_norm_lue %>%
  filter(sitename!="JP-MBF" & sitename!="US-NR1" & sitename!="US-PFa")%>% ##remove the sites do not performed when lmer
  #update in May, 2023:adding the LUE
  group_by(Clim_PFTs, doy) %>% 
  dplyr::summarise(
    lue_obs_sd=sd(lue,na.rm = T),
    lue_obs_mean=mean(lue,na.rm=T),
    #add sos_mean and eos_mean
    sos_mean=round(mean(sos,na.rm=T),0),
    peak_min=round(min(peak,na.rm=T),0),
    fapar = mean(fapar_itpl, na.rm = TRUE),
    across(starts_with("lue"),mean,na.rm=TRUE)
  ) %>%
  #add in May,2023
  mutate(avg_period = ifelse(doy >= sos_mean & doy <= peak_min, TRUE, FALSE))
#------------------
#update in May,2023-->add LUE plots:
#------------------
plot_final_lue_1 <- df_meandoy_Clim_PFTs_lue %>% 
  # filter(!is.nan(gpp_lmer) & !is.infinite(gpp_lmer))%>% ##this filter is important
  mutate(lue_obs=lue,lue=NULL)%>%
  pivot_longer(c(lue_obs, lue_pmodel,lue_temp_vpd,lue_lmer), names_to = "model", values_to = "lue") %>%
  mutate(model = fct_relevel(model, "lue_obs", "lue_pmodel","lue_temp_vpd","lue_lmer")) %>%
  # dplyr::filter((model %in% c( "lue_obs", "lue_pmodel","lue_lmer"))) %>% ##select the lmer
  dplyr::filter((model %in% c( "lue_obs", "lue_pmodel"))) %>% 
  ggplot() +
  geom_line(aes(x = doy, y = lue, color = model),size=0.8) +
  geom_ribbon(aes(x=doy,ymin=lue_obs_mean - lue_obs_sd,ymax=lue_obs_mean + lue_obs_sd),fill="black",alpha=0.3)+
  labs(y = expression( paste("LUE (g C ppfd"^-1*" m"^-2, " d"^-1, ")" ) ),
       x = "DoY") +
  # facet_wrap( ~Clim_PFTs, ncol = 3, scales = "free_y" ) +
  facet_wrap( ~Clim_PFTs) +
  ylim(0,0.05)+
  # theme_gray() +
  # theme(legend.position = "bottom") +
  # scale_color_manual(
  #   name="Model: ",
  #   values=c("black", "red", "royalblue", "darkgoldenrod", "springgreen", "orchid4"))+
  scale_color_manual("LUE sources",values = c("lue_obs" = "black",
                                              "lue_pmodel" = "orange"),
                                              # "lue_lmer" = "brown2"),
                     labels = c(expression(LUE[obs]),expression(LUE[Pmodel]),expression(LUE[LME])))+
  # labels = c("Observations","P-model","LME")) +
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
  )+
  theme(legend.text.align = 0)
ggsave("./manuscript/figures/FigureS_gpp_meandoy_forClimPFTs_lue.png",
       plot_final_lue_1,width = 15,height = 10)


#
plot_final_lue_2 <- df_meandoy_norm_Clim_PFTs_lue %>% 
  # filter(!is.nan(gpp_lmer) & !is.infinite(gpp_lmer))%>% ##this filter is important
  mutate(lue_obs=lue,lue=NULL)%>%
  pivot_longer(c(lue_obs, lue_pmodel,lue_temp_vpd,lue_lmer), names_to = "model", values_to = "lue") %>%
  mutate(model = fct_relevel(model, "lue_obs", "lue_pmodel","lue_temp_vpd","lue_lmer")) %>%
  # dplyr::filter((model %in% c( "lue_obs", "lue_pmodel","lue_lmer"))) %>% ##select the lmer
  dplyr::filter((model %in% c( "lue_obs", "lue_pmodel"))) %>% ##select the lmer
  ggplot() +
  geom_line(aes(x = doy, y = lue, color = model),size=0.8) +
  geom_ribbon(aes(x=doy,ymin=lue_obs_mean - lue_obs_sd,ymax=lue_obs_mean + lue_obs_sd),fill="black",alpha=0.3)+
  labs(y = expression( paste("LUE (g C ppfd"^-1*" m"^-2, " d"^-1, ")" ) ),
       x = "DoY") +
  # facet_wrap( ~Clim_PFTs, ncol = 3, scales = "free_y" ) +
  facet_wrap( ~Clim_PFTs) +
  ylim(0,0.05)+
  # theme_gray() +
  # theme(legend.position = "bottom") +
  # scale_color_manual(
  #   name="Model: ",
  #   values=c("black", "red", "royalblue", "darkgoldenrod", "springgreen", "orchid4"))+
  scale_color_manual("LUE sources",values = c("lue_obs" = "black",
                                              "lue_pmodel" = "orange"), 
                                            # "lue_lmer" = "brown2"),
                     labels = c(expression(LUE[obs]),expression(LUE[Pmodel]),expression(LUE[LME])))+
  # labels = c("Observations","P-model","LME")) +
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
  )+
  theme(legend.text.align = 0)  #align the legend (all the letter start at the same positoin)
ggsave("./manuscript/figures/FigureS_gpp_meandoy_norm_forClimPFTs_lue.png",
       plot_final_lue_2,width = 15,height = 10)


##following the Wolfhart's suggestion, adding the mean bias betweeen model and obs -> add May,2023
gpp_mean_bias <- df_meandoy_norm_Clim_PFTs%>%
  filter(avg_period==TRUE)%>%
  group_by(Clim_PFTs)%>%
  dplyr::summarise(gpp_bias_pmodel_greenup=mean(gpp_bias_pmodel,na.rm=T),
                   gpp_bias_lmer_greenup=mean(gpp_bias_lmer,na.rm=T))

#merge gpp_mean_bias to datasets:
df_meandoy_norm_Clim_PFTs<-left_join(df_meandoy_norm_Clim_PFTs,
          gpp_mean_bias)

#Figure for Clim.-PFTs
plot_final <- df_meandoy_norm_Clim_PFTs %>% 
  filter(!is.nan(gpp_lmer) & !is.infinite(gpp_lmer))%>% ##this filter is important
  pivot_longer(c(gpp_obs, gpp_pmodel, gpp_lue_const, gpp_temp_vpd,gpp_lmer), names_to = "model", values_to = "gpp") %>%
  mutate(model = fct_relevel(model, "gpp_obs", "gpp_pmodel", "gpp_lue_const", "gpp_temp_vpd","gpp_lmer")) %>%
  dplyr::filter((model %in% c( "gpp_obs", "gpp_pmodel","gpp_lmer"))) %>% ##select only one lm
  ggplot() +
  geom_line(aes(x = doy, y = gpp, color = model),size=0.8) +
  geom_ribbon(aes(x=doy,ymin=gpp_obs_mean - gpp_obs_sd,ymax=gpp_obs_mean + gpp_obs_sd),fill="black",alpha=0.3)+
  #adding the gpp bias-->May,2023:
  geom_segment(aes(x=120,y=0,xend=120+gpp_bias_pmodel_greenup*50,yend=0),col="orange",size=1.1)+
  geom_segment(aes(x=120,y=-1,xend=120+gpp_bias_lmer_greenup*50,yend=-1),col="brown2",size=1.1)+
  geom_segment(aes(x=120,y=-1.8,xend=120,yend=0.8),col=adjustcolor("black",0.8),size=1.02)+
  #adding the bias value:
  geom_text(aes(x=160,y=1,label=round(gpp_bias_pmodel_greenup,2)))+
  geom_text(aes(x=160,y=-2,label=round(gpp_bias_lmer_greenup,2)))+
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
     "gpp_pmodel" = "orange", "gpp_lmer" = "brown2"),
     labels = c(expression(GPP[obs]),expression(GPP[Pmodel]),expression(GPP[LME])))+
     # labels = c("Observations","P-model","LME")) +
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
  )+
theme(legend.text.align = 0)  #align the legend (all the letter start at the same positoin)
##adding the site number in each Clim-PFTs panel:
#update using original p-model
##
nsites<-ddf_norm %>%
  group_by(Clim_PFTs)%>%
  dplyr::summarise(nsite=length(unique(sitename)))
nsites$label<-paste0("N = ",nsites$nsite)
sites_num.info<-data.frame(
  doy=rep(20,nrow(nsites)),
  gpp=rep(14,nrow(nsites)),
  nsites
)
#function to add the site number info:
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
plot_final_addN<-tag_facet(plot_final,x=sites_num.info$doy,y=sites_num.info$gpp,
                           tag_pool = sites_num.info$label,size=5)
#
ggsave("./manuscript/figures/Figure2_gpp_meandoy_norm_forClimPFTs.png",plot_final_addN,width = 15,height = 10)


# Residual analysis for reply to reviewer 1 ----------------------------
tmp <- df_meandoy_norm_Clim_PFTs |> 
  mutate(res_pmodel = gpp_pmodel - gpp_obs,
         res_lmm = gpp_lmer - gpp_obs) |> 
  mutate(spring = ifelse(doy < 365/2, TRUE, FALSE)) |> 
  mutate(fapar_bin = cut(fapar, breaks = seq(0, 1, by = 0.1)))

## Residual vs fAPAR bin Clim-PFT -----------------
### P-model bias ----------------
tmp |> 
  ggplot(aes(x = fapar_bin, y = res_pmodel, fill = spring)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_fill_manual(
    name="Spring",
    values=c("#777055ff", "#29a274ff")
  ) +
  facet_wrap( ~Clim_PFTs, ncol = 2 ) +
  labs() +
  labs(y = expression( paste("P-model GPP residual (g C m"^-2, " d"^-1, ")" ) ),
       x = "fAPAR bin") +
  theme_classic()

ggsave("./manuscript/figures/FigADD_residual_fapar_pmodel_ClimPFT.png", width = 12, height = 12)

