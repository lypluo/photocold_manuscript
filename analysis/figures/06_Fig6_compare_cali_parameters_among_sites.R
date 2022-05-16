##---------------------------------------
#Aim: To compare the parameters among different groups(e.g.PFTs) after calibrating
#parameters for each site
##---------------------------------------
library(dplyr)
devtools::load_all("D:/Github/rbeni/")
library(rbeni) #-->make the evaluation plot
library(tidyverse)
library(cowplot)
library(grid)
library(ggpubr)
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
#---
#check the variables distribution and boxplots
#---
vars.names<-c("a1","b1","a2","b2","e","f","k")
for (i in 1:length(vars.names)) {
  hist(as.numeric(unlist(df_final_new[,vars.names[i]])),xlab = vars.names[i])
}

##----------boxplot---------------------
#a. first for site-level parameters
data_sel_sites<-df_final_new %>%
  select(sitename,classid,a1:k)%>%
  pivot_longer(c(a1:k),names_to = "parameter",values_to = "parameter_value")
#only focus on the a,b,c,d,k
data_sel_sites<-data_sel_sites %>%
  filter(parameter %in% c("a1","b1","a2","b2","k"))%>%
  mutate(PFT=classid,
         classid=NULL)
data_sel_sites$flag=rep("site",nrow(data_sel_sites))
#
data_sel_sites$parameter<-factor(data_sel_sites$parameter,
                                 levels = c("a1","b1","a2","b2","k"))
#-----------
#b.also load the parameters for diff PFTs:
load(paste0("./data/model_parameters/parameters_MSE_add_baseGDD/","optim_par_run5000_beni_PFTs.rds"))
paras_PFTs<-data.frame(DBF=par_PFTs$DBF,
                       MF=par_PFTs$MF,
                       EN=par_PFTs$ENF)
paras_PFTs<-as.data.frame(t(paras_PFTs))
paras_PFTs$PFT<-c("DBF","MF","ENF")
#also change the parameters names:
names(paras_PFTs)<-c("a1","b1","a2","b2","e","f","k","PFT")
#
data_sel_PFTs<-paras_PFTs %>%
  select(a1:b2,k,PFT)%>%
  pivot_longer(c(a1,b1,a2,b2,k),names_to = "parameter",values_to = "parameter_value")
data_sel_PFTs$flag=rep("PFT",nrow(data_sel_PFTs))
#
data_sel_PFTs$parameter<-factor(data_sel_PFTs$parameter,
                                levels = c("a1","b1","a2","b2","k"))
#c.load the parameters for diff Clim-PFTs:
load(paste0("./data/model_parameters/parameters_MSE_add_baseGDD/","optim_par_run5000_beni_Clim_andPFTs_update.rds"))
#
paras_Clim_PFTs<-c()
N<-length(names(par_Clim_PFTs))
for(i in 1:N){
  temp<-t(as.data.frame(par_Clim_PFTs[i]))
  paras_Clim_PFTs<-rbind(paras_Clim_PFTs,temp)
}
paras_Clim_PFTs<-as.data.frame(paras_Clim_PFTs)
names(paras_Clim_PFTs)<-c("a1","b1","a2","b2","e","f","k")
#
paras_Clim_PFTs$Clim_PFTs<-names(par_Clim_PFTs)
#
data_sel_Clim_PFTs<-paras_Clim_PFTs %>%
  select(a1:b2,k,Clim_PFTs)%>%
  pivot_longer(c(a1,b1,a2,b2,k),names_to = "parameter",values_to = "parameter_value")
data_sel_Clim_PFTs$flag=rep("Clim-PFT",nrow(data_sel_Clim_PFTs))

#########################
#merge site and PFTs data
data_sel_final<-bind_rows(data_sel_sites,data_sel_PFTs)
data_sel_final<-bind_rows(data_sel_final,data_sel_Clim_PFTs)
  
#box plot with gitter 
data_sel_final$PFT<-factor(data_sel_final$PFT,levels = c("DBF","MF","ENF"))
##parameters distributions of different sites
#
# devtools::install_github("zeehio/facetscales")
# library(facetscales)
# scales_y<-list(
#   "a"=scale_y_continuous(limits = c(-40,25)),
#   "b"=scale_y_continuous(limits = c(-5,25)),
#   "c"=scale_y_continuous(limits = c(-0,200)),
#   "d"=scale_y_continuous(limits = c(-5,20)),
#   "k"=scale_y_continuous(limits = c(-10,15))
# )

para_sites<-ggplot(data=data_sel_final[data_sel_final$flag=="site",],aes(x=parameter,y=parameter_value,fill=PFT,col=PFT))+
  geom_point(position = position_jitterdodge())+
  geom_boxplot(alpha=0.6)+
  xlab("")+
  facet_wrap(~parameter,scales = "free",nrow = 3)+
  xlab("Parameters")+
  ylab("")+
  theme_bw()+
  theme(legend.position = c(0.75,0.18),
        legend.background = element_blank(),
        legend.title = element_text(size=24),
        legend.text = element_text(size=22),
        axis.title = element_text(size=24),
        axis.text.y = element_text(size=22),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 22)) ##change the facet label size
#---------------------------------------------
tag_facet <- function(p, open = "", close = "", tag_pool = letters, x = -Inf, y = Inf, 
                      hjust = "", vjust = "", fontface = 2, family = "", ...) {
  
  # p<-para_sites
  # x=paras_PFTs_new$x
  # y=paras_PFTs_new$parameter_value
  # open=""
  # close=""
  # tag_pool=paras_PFTs_new$label
  
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
  p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE) 
}
library(egg)

##add the PFTs parameters onto each panel:
paras_PFTs_new<-data_sel_final[data_sel_final$flag=="PFT",-1]
paras_PFTs_new$x<-rep(NA,nrow(paras_PFTs_new))
paras_PFTs_new[paras_PFTs_new$PFT=="MF",]$x<-1
paras_PFTs_new[paras_PFTs_new$PFT=="DBF",]$x<-0.75
paras_PFTs_new[paras_PFTs_new$PFT=="ENF",]$x<-1.25
##
paras_PFTs_new$label<-rep("*",nrow(paras_PFTs_new))
paras_PFTs_new$col<-c(rep("red",5),rep("forestgreen",5),rep("blue",5))
#
paras_PFTs_new$PFT<-factor(paras_PFTs_new$PFT,levels = c("DBF","MF","ENF"))
paras_PFTs_new$parameter<-factor(paras_PFTs_new$parameter,levels = c("a1","b1","a2","b2","k"))
paras_boxplot<-tag_facet(para_sites,x=paras_PFTs_new$x,y=paras_PFTs_new$parameter_value,
                         #here I add 1 for y axis since there parameters for PFTs did not display properly
                           tag_pool = paras_PFTs_new$label,size=10,col=paras_PFTs_new$col)

#save the plot
save.path<-"./manuscript/figures/"
ggsave(paste0(save.path,"Figure6_parameters.png"),paras_boxplot)

#----------------scatter plot------------------
#environmental drivers vs parameters
library(ggforce)
plot_paras<-function(df_meteo,df_paras,Env_var,para,do_legend){
  # df_meteo<-df_final_new
  # df_paras<-data_sel_final
  # Env_var<-"tmean"
  # para<-"a1"
  # do_legend=FALSE
  # for example: Tmin vs a1
  #I.site-level
  df_site_level<-df_meteo %>%
    select(sitename,classid,koeppen_code,Clim.PFTs,para,tmin,temp)
  names(df_site_level)<-c("sitename","PFT","Clim.","Clim.PFTs","para","tmin","tmean")
  #
  t_pos<-match(Env_var,names(df_site_level))
  df_site_level_new<-df_site_level
  names(df_site_level_new)[t_pos]<-"Env_var"
  #II.Clim-PFT level---
  df_Clim_PFT_level<-df_site_level %>%
    group_by(Clim.PFTs)%>%
    summarise(tmin=mean(tmin,na.rm=T),
              tmean=mean(tmean,na.rm=T))
  #
  par_Clim_PFT_level<-df_paras[df_paras$flag=="Clim-PFT",]
  par_Clim_PFT_level<-par_Clim_PFT_level %>%
    filter(parameter==para)%>%
    mutate(sitename=NULL,PFT=NULL,flag=NULL,
           Clim.PFTs=Clim_PFTs,Clim_PFTs=NULL,parameter=NULL)
  names(par_Clim_PFT_level)<-c("para","Clim.PFTs")
  #
  df_Clim_PFT_level_new<-left_join(df_Clim_PFT_level,par_Clim_PFT_level)
  df_Clim_PFT_level_new$Clim<-substr(df_Clim_PFT_level_new$Clim.PFTs,1,3)
  df_Clim_PFT_level_new$Clim<-factor(df_Clim_PFT_level_new$Clim,
                                     levels = c("Cfa","Cfb","Dfb","Dfc"))
  df_Clim_PFT_level_new$PFT<-substr(df_Clim_PFT_level_new$Clim.PFTs,5,7)
  df_Clim_PFT_level_new$PFT<-factor(df_Clim_PFT_level_new$PFT,levels = c("DBF","MF","ENF"))
  #
  t_pos<-match(Env_var,names(df_Clim_PFT_level_new))
  names(df_Clim_PFT_level_new)[t_pos]<-"Env_var"
  
  #III.PFT level---
  df_PFT_level<-df_site_level %>%
    group_by(PFT)%>%
    summarise(tmin=mean(tmin,na.rm=T),
              tmean=mean(tmean,na.rm=T))
  par_PFT_level<-df_paras[df_paras$flag=="PFT",] %>%
    mutate(sitename=NULL,Clim_PFTs=NULL,flag=NULL)%>%
    filter(parameter==para)
  par_PFT_level<-par_PFT_level[,-1]
  names(par_PFT_level)<-c("para","PFT")
  #
  df_PFT_level_new<-left_join(df_PFT_level,par_PFT_level)
  #
  t_pos<-match(Env_var,names(df_PFT_level_new))
  names(df_PFT_level_new)[t_pos]<-"Env_var"
  
  ##----plotting----##
  pars_final<-ggplot()+
    geom_point(data=df_site_level_new,aes(x=Env_var,y=para,shape=PFT,size="site"))+
    ggforce::geom_mark_ellipse(data=df_site_level_new,
    aes(x=Env_var,y=para,label=PFT,shape=PFT),label.fill = "grey",
    con.border = "one",con.cap = 0,con.size = 1.1,
    con.arrow = grid::arrow(angle=30,ends = "last",length = unit(0.1,"inches")))+  ##site-level
    geom_point(data=df_Clim_PFT_level_new,
               aes(x=Env_var,y=para,shape=PFT,col=Clim,size="Clim-PFT"),alpha=0.8)+ ##Clim-PFT level
    geom_point(data=df_PFT_level_new,
               aes(x=Env_var,y=para,shape=PFT,size="PFT"),col="blue",alpha=0.4)+##PFT level
    scale_size_manual("",values = c("site"=5,"Clim-PFT"=8,"PFT"=15))+
    xlab(paste0(Env_var," (°C)"))+
    ylab(paste0(para," (°C)"))+
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
  if(do_legend==FALSE){
    pars_final<-pars_final+
      theme(legend.position = "none")
  }
  #
  return(pars_final)
}
############
# df_meteo<-df_final_new
# df_paras<-data_sel_final
# Env_var<-"tmin"
# para<-"a1"
# do_legend=FALSE

p_tmin_a1<-plot_paras(df_meteo = df_final_new,df_paras = data_sel_final,Env_var = "tmin",
           para = "a1",FALSE)
p_tmean_b1<-plot_paras(df_meteo = df_final_new,df_paras = data_sel_final,Env_var = "tmean",
                      para = "b1",FALSE)  
p_tmin_a2<-plot_paras(df_meteo = df_final_new,df_paras = data_sel_final,Env_var = "tmin",
                      para = "a2",FALSE)
p_tmean_b2<-plot_paras(df_meteo = df_final_new,df_paras = data_sel_final,Env_var = "tmean",
                       para = "b2",FALSE)  
p_tmean_k<-plot_paras(df_meteo = df_final_new,df_paras = data_sel_final,Env_var = "tmean",
                      para = "k",TRUE)  
#change the x labels:
p_tmin_a1<-p_tmin_a1+xlab(expression("T"[min]*" (°C)"))
p_tmean_b1<-p_tmean_b1+xlab(expression("T"[mean]*" (°C)"))
p_tmin_a2<-p_tmin_a2+xlab(expression("T"[min]*" (°C)"))
p_tmean_b2<-p_tmean_b2+xlab(expression("T"[mean]*" (°C)"))
p_tmean_k<-p_tmean_k+xlab(expression("T"[mean]*" (°C)"))+
  theme(legend.position = c(1.7,0.4),
        legend.spacing.y = unit(0.4, 'cm'),
        legend.key.size = unit(2,"cm"))+
  guides(size=guide_legend(ncol=3,byrow=T))+
  guides(col=guide_legend(ncol = 4,byrow = T,override.aes = list(size = 6)))+
  guides(shape=guide_legend(ncol = 3,byrow = T,override.aes = list(size = 6)))
  
#merge the plots:
paras_range<-cowplot::plot_grid(p_tmin_a1,p_tmean_b1,p_tmin_a2,p_tmean_b2,p_tmean_k,
          nrow = 3,ncol=2,labels = "auto",label_size = 20,align = "hv")
######save the plot###########
save.path<-"./manuscript/figures/"
ggsave(paste0(save.path,"Figure7_parameters_ranges.png"),paras_range,height = 18,width = 20)

#############################additional code ###########################
#----
#check the paraters difference among different groups
#----
library(ggpubr)
library(cowplot)
check_groups<-function(df,par_name){
  # df<-df_final_new
  # par_name<-"a1"

  df_t<-df %>%
    select(sitename,classid,koeppen_code,Clim.PFTs,par_name)
  names(df_t)<-c("sitename","classid","koeppen_code","Clim.PFTs","par")
  #for different PFTs
  p_PFTs<-ggplot(data=df_t,aes(x=par,color=classid,fill=classid))+
    geom_histogram(aes(y=..density..,),position = "identity",binwidth = 1,alpha=0.5)+
    geom_density(alpha=.2)+
    xlab(par_name)
  #for different Clim.
  p_Clim<-ggplot(data=df_t,aes(x=par,color=koeppen_code,fill=koeppen_code))+
    geom_histogram(aes(y=..density..,),position = "identity",binwidth = 1,alpha=0.5)+
    geom_density(alpha=.2)+
    xlab(par_name)
  #for different Clim.-PFTs
  p_Clim.PFTs<-ggplot(data=df_t,aes(x=par,color=Clim.PFTs,fill=Clim.PFTs))+
    geom_histogram(aes(y=..density..,),position = "identity",binwidth = 1,alpha=0.5)+
    xlab(par_name)
  # geom_density(alpha=.2)
  #
  p_merge<-plot_grid(p_PFTs,p_Clim,p_Clim.PFTs)
  return(p_merge)
}
##
check_groups(df_final_new,"a1")
check_groups(df_final_new,"b1")
check_groups(df_final_new,"a2")
check_groups(df_final_new,"b2")
check_groups(df_final_new,"e")
check_groups(df_final_new,"f")
check_groups(df_final_new,"k")

#----
#check environmental drivers relationship between parameters
#----
check_relation<-function(df,par_name){
  # df<-df_final_new
  # par_name<-"a1"

  df_t<-df %>%
    select(sitename:fapar_spl,classid,koeppen_code,Clim.PFTs,par_name)
  names(df_t)<-c("sitename","temp","prec","vpd","ppdf","elv",
                 "tmin","tmax","fapar_itpl","fapar_spl",
                 "classid","koeppen_code","Clim.PFTs","par")
  #
  p_ta<-ggplot(data=df_t,aes(x=temp,y=par,color=classid))+
    geom_point()+
    xlab("ta")+
    ylab(par_name)
  #
  p_prec<-ggplot(data=df_t,aes(x=prec,y=par,color=classid))+
    geom_point()+
    xlab("prec")+
    ylab(par_name)
  p_vpd<-ggplot(data=df_t,aes(x=vpd,y=par,color=classid))+
    geom_point()+
    xlab("vpd")+
    ylab(par_name)
  p_ppfd<-ggplot(data=df_t,aes(x=ppdf,y=par,color=classid))+
    geom_point()+
    xlab("ppfd")+
    ylab(par_name)
  p_tmin<-ggplot(data=df_t,aes(x=tmin,y=par,color=classid))+
    geom_point()+
    xlab("tmin")+
    ylab(par_name)
  p_tmax<-ggplot(data=df_t,aes(x=tmax,y=par,color=classid))+
    geom_point()+
    xlab("tmax")+
    ylab(par_name)
  p_fapar<-ggplot(data=df_t,aes(x=fapar_itpl,y=par,color=classid))+
    geom_point()+
    xlab("fapar_itpl")+
    ylab(par_name)
    #
  p_merge<-plot_grid(p_ta,p_prec,p_vpd,
                     p_ppfd,p_tmin,p_fapar,nrow = 2,align = "hv")
  return(p_merge)
}

#
check_relation(df_final_new,"a1")
check_relation(df_final_new,"b1")
check_relation(df_final_new,"a2")
check_relation(df_final_new,"b2")
check_relation(df_final_new,"e")
check_relation(df_final_new,"f")
check_relation(df_final_new,"k")

