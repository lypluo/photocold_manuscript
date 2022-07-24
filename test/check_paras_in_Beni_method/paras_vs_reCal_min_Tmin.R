##---------------------------------------
#Aim: check how the parameter(e.g. a1,a2) variation along the temperature gradient
##using the parameter calibrated from Beni's method-->
##In this script-->to compare the paras vs Tmin, recalculate the Timin in the early season:
#before the green-up(b_GP), green-up(GP), early season (ES;b_GP+GP period)
##---------------------------------------
library(dplyr)
# devtools::load_all("D:/Github/rbeni/")
# library(rbeni) #-->make the evaluation plot
library(tidyverse)
library(cowplot)
library(grid)
library(ggpubr)
library(lubridate)
#---------------------------
#(1)load the calibrated parameters for each site
#---------------------------
load(paste0("./data/model_parameters/parameters_MSE_add_baseGDD/","optim_par_run5000_beni_eachsite_updated.rds"))

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
rownames(pars_final)<-c("a1","b1","a2","b2","e","f","k")

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
source(paste0("./R/functions_in_model/model_hardening_byBeni_addbaseGDD_rev.R"))
#a.get the stress factor(calibration factor) for each site
df_final<-c()
df_recent$doy<-yday(df_recent$date)
for (i in 1:length(sel_sites)) {
  df_sel<-df_recent %>%
    dplyr::filter(sitename==sel_sites[i])
  
  scaling_factors <- df_sel %>%
    # group_by(sitename, year) %>%
    do({
      scaling_factor <- model_hardening_2par(.,par_mutisites[[i]])
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

#----
#load the phenology data
#----
pheno.path<-"./data/event_length/"
load(file=paste0(pheno.path,"df_events_length.RDA"))
df.pheno<-df_events_all %>%
      select(sitename:peak)%>%
      mutate(year=Year,Year=NULL)
#----merge the data------
df_merge_new<-left_join(df_merge_new,df.pheno,by=c("sitename","year"))
df_avg<-df_merge_new %>%
  group_by(sitename)%>%
  summarise(sos_mean=mean(sos,na.rm=T),
            sos_median=median(sos,na.rm = T))
hist(x=unlist(df_avg[,"sos_mean"]))
hist(x=unlist(df_avg[,"sos_median"]))
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
#adding in July,22
df_sum_yearly_3<-df_merge_new%>%
  group_by(sitename,year)%>%
  filter(doy>=1 & doy<=60)%>%
  dplyr::summarise(tmin_b_GP=min(tmin,na.rm=T))
df_sum_yearly_4<-df_merge_new%>%
  group_by(sitename,year)%>%
  filter(doy>=61 & doy<=182)%>%
  dplyr::summarise(tmin_GP=min(tmin,na.rm=T))
df_sum_yearly_5<-df_merge_new%>%
  group_by(sitename,year)%>%
  filter(doy>=1 & doy<=182)%>%
  dplyr::summarise(tmin_ES=min(tmin,na.rm=T))
#
df_sum_add<-left_join(df_sum_yearly_3,df_sum_yearly_4,by=c("sitename","year"))
df_sum_add<-left_join(df_sum_add,df_sum_yearly_5,by=c("sitename","year"))  
#------------merging------------
df_sum_yearly<-left_join(df_sum_yearly,df_sum_add,by=c("sitename","year"))

#---
#summary site-years for site
#---
df_sum_1<-df_sum_yearly %>%
  group_by(sitename) %>%
  summarise_at(vars(temp:fapar_spl,tmin_b_GP:tmin_ES),mean,na.rm=T)
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
library(ggpmisc)
plot_data<-df_sum_new %>%
  # dplyr::select(sitename,temp,prec,tmin,fT,classid,Clim.PFTs,a1:k)%>%
  mutate(PFT=classid,classid=NULL)

#--------------------
#plots:environmental drivers vs parameters
#--------------------
#1)for all the sites:
plot_para_each<-function(df,Tvar,para,do_legend){
  # df<-plot_data
  # Tvar<-"tmin_b_GP"
  # para<-"a2"
  # do_legend=FALSE
  
  df<-df %>%
    select(sitename,Tvar,PFT,para)
  names(df)<-c("sitename","tmin","PFT","para")
  
  df_plot<-df %>%
    ggplot(aes(x=tmin,y=para,col=PFT))+
    geom_point(size=5)+
    # geom_smooth(method = lm)+
    scale_color_manual(values = c("DBF"="orange","MF"="cyan","ENF"="magenta"))+
    xlab(paste0(Tvar,"(°C)"))+
    ylab(paste0(para," (°C)"))+
    theme(
      legend.position = c(0.15,0.8),
      legend.text = element_text(size=22),
      legend.key.size = unit(2, 'lines'),
      axis.title = element_text(size=26),
      axis.text = element_text(size = 22),
      text = element_text(size=24),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(colour ="grey",fill="white"))
  if(do_legend==FALSE){
    df_plot<-df_plot+
      theme(legend.position = "none")
  }
  return(df_plot)
}
################only test the parameter a2#############
p_tmin_a2<-plot_para_each(plot_data,"tmin","a2",TRUE)
p_tmin_b_GP_a2<-plot_para_each(plot_data,"tmin_b_GP","a2",TRUE)
p_tmin_GP_a2<-plot_para_each(plot_data,"tmin_GP","a2",TRUE)
p_tmin_ES_a2<-plot_para_each(plot_data,"tmin_ES","a2",TRUE)

#merge the plots:
paras_general_range<-cowplot::plot_grid(p_tmin_a2,p_tmin_b_GP_a2,
      p_tmin_GP_a2,p_tmin_ES_a2,
      nrow=2, ncol = 2,labels = "auto",label_size = 20,align = "hv")
#save the plot
save.path<-"./manuscript/test_files/check_paras_in_Beni_method/"
ggsave(paste0(save.path,"Check_paras_general_beni_method_min_Tmin_vs_a2.png"),
       paras_general_range,width = 18,height = 15)

#2)for sites and groups(PFT):
library(ggforce)
library(ggrepel)
plot_paras<-function(df,Env_var,para,do_legend){
  # df<-plot_data
  # Env_var<-"tmin_b_GP"
  # para<-"a2"
  # do_legend=FALSE
  # for example: tmin_b_GP vs a2
  #I.site-level
  df_site_level<-df %>%
    dplyr::select(sitename,PFT,koeppen_code,Clim.PFTs,para,Env_var,temp)
  names(df_site_level)<-c("sitename","PFT","Clim.","Clim.PFTs","para","Env_var","tmean")
  #
  t_pos<-match(Env_var,names(df_site_level))
  df_site_level_new<-df_site_level
  names(df_site_level_new)[t_pos]<-"Env_var"
  
  #III.PFT level---
  df_PFT_level<-df_site_level%>%
    mutate(PFT=factor(PFT,levels = c("DBF","MF","ENF")))%>%
    group_by(PFT)%>%
    dplyr::summarise(Env_var=mean(Env_var,na.rm=T),
                     tmean=mean(tmean,na.rm=T),
                     para=mean(para,na.rm = T))
  
  ##----plotting----##
  library(ggpmisc)  
  # library(ggpubr)
  #linear regression:
  #----DBF------
  lm_DBF<-lm(data=df_site_level_new[df_site_level_new$PFT=='DBF',],
             para~Env_var)
  stat_lm_DBF<-summary(lm_DBF)
  stat_DBF_label<-data.frame(r.squared=round(stat_lm_DBF$r.squared,2),
                             p.value=round(coef(stat_lm_DBF)[2,4],4))
  #----Dfc-ENF-----
  lm_Dfc_ENF<-lm(data=df_site_level_new[df_site_level_new$Clim.PFTs=='Dfc-ENF',],
                 para~Env_var)
  stat_lm_Dfc_ENF<-summary(lm_Dfc_ENF)
  stat_Dfc_ENF_label<-data.frame(r.squared=round(stat_lm_Dfc_ENF$r.squared,2),
                                 p.value=round(coef(stat_lm_Dfc_ENF)[2,4],4))
  
  pars_final<-ggplot()+
    geom_point(data=df_site_level_new,aes(x=Env_var,y=para,col=PFT),size=3)+
    # scale_color_discrete_sequential(palette = "Viridis")+
    geom_text_repel(data=df_site_level_new,aes(x=Env_var,y=para,label=sitename),size=5)+
    stat_poly_line(data=df_site_level_new[df_site_level_new$PFT=='DBF',],
                   aes(x=Env_var,y=para,col=PFT),
                   fill=adjustcolor("goldenrod1"),method = "lm",formula = y ~ x,lty=2)+
    # stat_poly_eq(data=df_site_level_new[df_site_level_new$PFT=='DBF',],
    #                aes(x=Env_var,y=para,col=PFT,
    #                    label = paste(
    #                                  # after_stat(grp.label), "*\"：\"*",
    #                                  # after_stat(eq.label), "*\", \"*",
    #                                  after_stat(rr.label), 
    #                                  after_stat(p.value.label),
    #                                  sep = "*\", \"*"),
    #                    label.x=0.5,label.y="bottom"))+
    # ggforce::geom_mark_ellipse(data=df_site_level_new[df_site_level_new$PFT=='DBF',],
    #                            aes(x=Env_var,y=para,label=PFT,group=PFT,col=PFT),label.fill = "goldenrod1",
    #                            con.border = "one",con.cap = 0,con.size = 1.1,con.colour = "goldenrod1",
    #                            con.arrow = grid::arrow(angle=30,ends = "last",length = unit(0.1,"inches")))+  ##DBF
    stat_poly_line(data=df_site_level_new[df_site_level_new$Clim.PFTs=='Dfc-ENF',],
                   aes(x=Env_var,y=para,col=PFT),fill=adjustcolor("magenta1"),
                   method = "lm",formula = y ~ x,lty=2)+
    # stat_poly_eq(data=df_site_level_new[df_site_level_new$Clim.PFTs=='Dfc-ENF',],
    #                  aes(x=Env_var,y=para,col=PFT,
    #                     label = paste(
    #                       after_stat(rr.label), 
    #                       after_stat(p.value.label),
    #                       sep = "*\", \"*"),
    #                 label.x=0.5,label.y="bottom"))+
    # ggforce::geom_mark_ellipse(data=df_site_level_new[df_site_level_new$Clim.PFTs=="Dfc-ENF",],
    #                            aes(x=Env_var,y=para,label=Clim.PFTs,group=Clim.PFTs,col=PFT),label.fill = "magenta1",
    #                            con.border = "one",con.cap = 0,con.size = 1.1,con.colour = "magenta1",
    #                            con.arrow = grid::arrow(angle=30,ends = "last",length = unit(0.1,"inches")))+  ##Dfc-ENF
    scale_color_manual(values = c("DBF"="orange","MF"="cyan","ENF"="magenta"))+
    xlab(paste0(Env_var," (°C)"))+
    ylab(paste0(para," (°C)"))+
    theme(
      legend.text = element_text(size=22),
      legend.position = c(0.15,0.8),
      legend.key.size = unit(2, 'lines'),
      axis.title = element_text(size=26),
      axis.text = element_text(size = 22),
      text = element_text(size=24),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(colour ="grey",fill="white")
    )
  if(para=="a1"){
    pars_final<-pars_final+
      annotate(geom = "text",x=-4.1,y=33,label = paste0("italic(R) ^ 2 == ",
                                                         stat_DBF_label$r.squared),parse=TRUE,col="orange",size=5)+
      annotate(geom = "text",x=1,y=33,label = paste0("italic(p) ==",
                                                      round(stat_DBF_label$p.value,2)),parse=TRUE,col="orange",size=5)+
      annotate(geom = "text",x=-4.1,y=29,label = paste0("italic(R) ^ 2 == ",
                                                         stat_Dfc_ENF_label$r.squared),parse=TRUE,col="magenta",size=5)+
      annotate(geom = "text",x=1,y=29,label = paste0("italic(p) == ",
                                                      round(stat_Dfc_ENF_label$p.value,2)),parse=TRUE,col="magenta",size=5)
  }
  ##tmin vs a2
  if(para=="a2" & Env_var=="tmin"){
    pars_final<-pars_final+
      annotate(geom = "text",x=-4.1,y=210,label = paste0("italic(R) ^ 2 == ",
                                                        stat_DBF_label$r.squared),parse=TRUE,col="orange",size=5)+
      annotate(geom = "text",x=1,y=210,label = paste0("italic(p) ==",
                                                     round(stat_DBF_label$p.value,2)),parse=TRUE,col="orange",size=5)+
      annotate(geom = "text",x=-4.1,y=192,label = paste0("italic(R) ^ 2 == ",
                                                        stat_Dfc_ENF_label$r.squared),parse=TRUE,col="magenta",size=5)+
      annotate(geom = "text",x=1,y=192,label = paste0("italic(p) == ",
                                                     round(stat_Dfc_ENF_label$p.value,2)),parse=TRUE,col="magenta",size=5)
  }
  if(para=="a2" & Env_var=="tmin_b_GP"){
    pars_final<-pars_final+
      annotate(geom = "text",x=-25.1,y=210,label = paste0("italic(R) ^ 2 == ",
                                                         stat_DBF_label$r.squared),parse=TRUE,col="orange",size=5)+
      annotate(geom = "text",x=-20,y=210,label = paste0("italic(p) ==",
                                                      round(stat_DBF_label$p.value,2)),parse=TRUE,col="orange",size=5)+
      annotate(geom = "text",x=-25.1,y=192,label = paste0("italic(R) ^ 2 == ",
                                                         stat_Dfc_ENF_label$r.squared),parse=TRUE,col="magenta",size=5)+
      annotate(geom = "text",x=-20,y=192,label = paste0("italic(p) == ",
                                                      round(stat_Dfc_ENF_label$p.value,2)),parse=TRUE,col="magenta",size=5)
  }
  if(para=="a2" & Env_var=="tmin_GP"){
    pars_final<-pars_final+
      annotate(geom = "text",x=-4.1,y=210,label = paste0("italic(R) ^ 2 == ",
                                                         stat_DBF_label$r.squared),parse=TRUE,col="orange",size=5)+
      annotate(geom = "text",x=1,y=210,label = paste0("italic(p) ==",
                                                      round(stat_DBF_label$p.value,2)),parse=TRUE,col="orange",size=5)+
      annotate(geom = "text",x=-4.1,y=192,label = paste0("italic(R) ^ 2 == ",
                                                         stat_Dfc_ENF_label$r.squared),parse=TRUE,col="magenta",size=5)+
      annotate(geom = "text",x=1,y=192,label = paste0("italic(p) == ",
                                                      round(stat_Dfc_ENF_label$p.value,2)),parse=TRUE,col="magenta",size=5)
  }
  if(para=="a2" & Env_var=="tmin_ES"){
    pars_final<-pars_final+
      annotate(geom = "text",x=-15.1,y=210,label = paste0("italic(R) ^ 2 == ",
                                                          stat_DBF_label$r.squared),parse=TRUE,col="orange",size=5)+
      annotate(geom = "text",x=-10,y=210,label = paste0("italic(p) ==",
                                                        round(stat_DBF_label$p.value,2)),parse=TRUE,col="orange",size=5)+
      annotate(geom = "text",x=-15.1,y=192,label = paste0("italic(R) ^ 2 == ",
                                                          stat_Dfc_ENF_label$r.squared),parse=TRUE,col="magenta",size=5)+
      annotate(geom = "text",x=-10,y=192,label = paste0("italic(p) == ",
                                                        round(stat_Dfc_ENF_label$p.value,2)),parse=TRUE,col="magenta",size=5)
  }
################
  
  if(para=="b1"){
    pars_final<-pars_final+
      annotate(geom = "text",x=-4.1,y=22,label = paste0("italic(R) ^ 2 == ",
                                                        stat_DBF_label$r.squared),parse=TRUE,col="orange",size=5)+
      annotate(geom = "text",x=1,y=22,label = paste0("italic(p) ==",
                                                     round(stat_DBF_label$p.value,2)),parse=TRUE,col="orange",size=5)+
      annotate(geom = "text",x=-4.1,y=20,label = paste0("italic(R) ^ 2 == ",
                                                        stat_Dfc_ENF_label$r.squared),parse=TRUE,col="magenta1",size=5)+
      annotate(geom = "text",x=1,y=20,label = paste0("italic(p) == ",
                                                     round(stat_Dfc_ENF_label$p.value,2)),parse=TRUE,col="magenta1",size=5)
  }
  if(para=="b2"){
    pars_final<-pars_final+
      annotate(geom = "text",x=-4.1,y=20,label = paste0("italic(R) ^ 2 == ",
                                                        stat_DBF_label$r.squared),parse=TRUE,col="orange",size=5)+
      annotate(geom = "text",x=1,y=20,label = paste0("italic(p) ==",
                                                     round(stat_DBF_label$p.value,2)),parse=TRUE,col="orange",size=5)+
      annotate(geom = "text",x=-4.1,y=18,label = paste0("italic(R) ^ 2 == ",
                                                        stat_Dfc_ENF_label$r.squared),parse=TRUE,col="magenta1",size=5)+
      annotate(geom = "text",x=1,y=18,label = paste0("italic(p) == ",
                                                     round(stat_Dfc_ENF_label$p.value,2)),parse=TRUE,col="magenta1",size=5)
  }
  if(para=="k"){
    pars_final<-pars_final+
      annotate(geom = "text",x=-10.1,y=15,label = paste0("italic(R) ^ 2 == ",
                                                      stat_DBF_label$r.squared),parse=TRUE,col="orange",size=5)+
      annotate(geom = "text",x=-5,y=15,label = paste0("italic(p) ==",
                                                     round(stat_DBF_label$p.value,2)),parse=TRUE,col="orange",size=5)+
      annotate(geom = "text",x=-10.1,y=13,label = paste0("italic(R) ^ 2 == ",
                                                      stat_Dfc_ENF_label$r.squared),parse=TRUE,col="magenta1",size=5)+
      annotate(geom = "text",x=-5,y=13,label = paste0("italic(p) == ",
                                                     round(stat_Dfc_ENF_label$p.value,2)),parse=TRUE,col="magenta1",size=5)
  }

  if(do_legend==FALSE){
    pars_final<-pars_final+
      theme(legend.position = "none")
  }
  #
  return(pars_final)
}
############
# df<-plot_data
# Env_var<-"tmin"
# para<-"a1"
# do_legend=FALSE

p_tmin_a2<-plot_paras(df = plot_data,Env_var = "tmin",
                      para = "a2",FALSE)  
p_tmin_b_GP_a2<-plot_paras(df = plot_data,Env_var = "tmin_b_GP",
                      para = "a2",FALSE)  
p_tmin_GP_a2<-plot_paras(df = plot_data,Env_var = "tmin_GP",
                           para = "a2",FALSE)  
p_tmin_ES_a2<-plot_paras(df = plot_data,Env_var = "tmin_ES",
                         para = "a2",FALSE)  



#change the x labels:

#merge the plots:
paras_range<-cowplot::plot_grid(p_tmin_a2,p_tmin_b_GP_a2,
                                p_tmin_GP_a2,p_tmin_ES_a2,
      nrow=2, ncol = 2,labels = "auto",label_size = 20,align = "hv")


#save the plot
save.path<-"./manuscript/test_files/check_paras_in_Beni_method/"
ggsave(paste0(save.path,"Check_paras_beni_method_min_Tmin_vs_a2.png"),
       paras_range,width = 18,height = 15)


#--------------------
#plots:fT vs mean
#--------------------
plot_fT<-ggplot()+
  geom_point(data = plot_data,aes(x=temp,y=fT,col=PFT,size=prec))+
  geom_text_repel(data = plot_data,aes(x=temp,y=fT,col=PFT,label=sitename),size=4)+
  scale_color_manual(values = c("DBF"="orange","MF"="cyan","ENF"="magenta"))+
  geom_smooth(data=plot_data,aes(x=temp,y=fT),col="blue",
              method = "lm",formula = y ~ x,lty=2,fill=adjustcolor("steelblue2",0.2))+
  stat_poly_eq(data=plot_data,
                 aes(x=temp,y=fT,
                     label = paste(
                                   after_stat(rr.label),
                                   after_stat(p.value.label),
                                   sep = "*\", \"*"),
                     ),col="blue")+
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

