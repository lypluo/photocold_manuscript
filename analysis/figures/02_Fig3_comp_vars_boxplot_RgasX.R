#############################################################################
#Aim:compare the variables in "event" and "non-event" site years after aligning the data-->specifically for Temp
#-------------------------------------------------------------------------
#(1)load the data that includes the "is_event" information
#-------------------------------------------------------------------------
#---------------------
#A.load the event_length data
#---------------------
load.path<-"./data/event_length/"
load(paste0(load.path,"df_events_length.RDA"))
# a function to separate out the site-year that event_length higher than some thresholds(e.g. 30 days)
sep_siteyears<-function(df,sep_var,sep_thrs){
  # df<-df_events_all
  # sep_var<-"Over_days_length"
  # sep_thrs<-30
  #
  df_sep<-df[df[,sep_var]>sep_thrs & !is.na(df[,sep_var]),]
  pos_sep<-as.vector(which(df[,sep_var]>sep_thrs))
  #as some site the sep_var==NA, hence using this way to separate the data
  pos_left<-setdiff(as.vector(1:length(df[,sep_var])),pos_sep)
  df_left<-df[pos_left,]
  df_new<-list(df_sep,df_left)
  names(df_new)<-c("event_siteyears","noevent_siteyears")
  return(df_new)
}
##separate the site-year when the over_days_length > 20 days
df.sep20<-sep_siteyears(df_events_all,"Over_days_length",20)

#---------------------
#B.load the data 
#---------------------
load.path<-"./data/data_used/"
#from new method:
load(paste0(load.path,"ddf_labeled_norm_trs_newmethod_all_overestimation_Fluxnet2015_sites.RDA"))
df_all_sites<-ddf_labeled;rm(ddf_labeled) 
df_norm_all<-df_all_sites
  
#-------------------------------------------------------------------------
#(2)start to align the data according to Beni's functions of "align_events" and "get_consecutive"
#-------------------------------------------------------------------------
#---------------------------------
#source the functions to detect the consecutive events
#---------------------------------
fun.path<-"./R/functions_from_beni/"
#source get_consecutive.R
source(paste0(fun.path,"get_consecutive.R"))
#source get_consecutive_greenup.R
source(paste0(fun.path,"get_consecutive_greenup.R"))

#source align_event.R -->for event site years
source(paste0(fun.path,"align_events_df.R"))
#source align_nonevent.R -->for non-event site years
source(paste0(fun.path,"align_nonevents_df.R"))
#-------------------------------------------------
#add additional code for PFTs: add PFTs information for each sites-->2021-12-13
#------------------------------------------------
#load the modis data for classifying the vegetation types:
library(tidyverse) #-->function read_rds
PFT.path<-"./data-raw/raw_data/sites_info/"
load(paste0(PFT.path,"Available_sites_info.RDA"))
#
df_norm_all<-df_norm_all%>%
  left_join(
    df_sites_avai,
    by = "sitename"
  )

sep_siteyears_data<-function(df.data,dovars,df.sep,leng_threshold,before,after,nbins,do_norm){
  # df.data<-df_norm_all
  # dovars<-c("gpp_obs")
  # df.sep<-df.sep20
  # leng_threshold<-5
  # before=30
  # after=0
  # nbins=10
  # do_norm=TRUE
  
  #-------------------------------------------------------------------------
  #A.separating the datasets:"event" and "nonevent" site years
  #-------------------------------------------------------------------------
  N<-nrow(df.sep$event_siteyears)
  df_events.years<-c()    #target
  df_nonevents.years<-c() #target
  #for is.event years
  pos.isevent<-c()
  for(i in 1:N){
    df.temp<-subset(df.data,sitename==df.sep$event_siteyears$sitename[i] & Year==df.sep$event_siteyears$Year[i])
    pos.temp<-c(which(df.data$sitename==df.temp$sitename & df.data$Year==df.temp$Year))
    
    df_events.years<-rbind(df_events.years,df.temp)
    pos.isevent<-c(pos.isevent,pos.temp)
  }
  #for no is.event years
  pos.non_isevent<-setdiff(c(1:nrow(df.data)),pos.isevent)
  df_nonevents.years<-df.data[pos.non_isevent,] #target
  
  #---------------------------------
  #B.start to align different events(for event site years) and green-up(for non-event site years) in each site-year
  #---------------------------------
  # df<-df_norm_all
  # leng_threshold<-20
  # before=30
  # after=30
  # nbins=10
  # do_norm=FALSE
  
  #function format
  # align_events <- function( df, leng_threshold, before, after, nbins, do_norm=FALSE )
  #at this stage-->bins do not used(now set default value as 10)
  #-------------
  #set different length_threshold-->5 days(consecutive 5 days overestimation will be an event)
  #-------------
  #and select the 30 days before the events 
  df_len_events<-align_events(df_events.years,dovars,leng_threshold = leng_threshold,before = before,after = after,
                              nbins = nbins,do_norm = do_norm)
  print("ok")
  df_len_nonevents<-align_nonevents(df_nonevents.years,dovars,leng_threshold = leng_threshold,before = before,after = after,
                                    nbins = nbins,do_norm = do_norm)
  #---------------------------------
  #C.separate the df_event and df_nonevent;
  #---------------------------------
  df_all<-c()
  #for event site years
  df_dday<-df_len_events$df_dday
  #for non_event site years, take the doy belongs to green-up period:
  df_noevent_dday<-df_len_nonevents$df_dday
  #
  df_all<-list(df_dday=df_dday,df_noevent_dday=df_noevent_dday)
  return(df_all)
}
#!!important step:leng_threshold=5-->merge the events(consecutive days are over 5 days)
#do_vars-->the variables that are going to be processed(normalized)-->actually do not processed in this process 
# names(df_norm_all)
# do_vars<-c("gpp_obs","fapar_itpl","fapar_spl",paste0(c("ppfd","temp_day","temp_min","temp_max",
#                                                        "vpd_day","prec","patm","SW_IN","ws",paste0("TS_",1:7),paste0("SWC_",1:5)),"_fluxnet2015"),
#            "gcc_90","rcc_90")
#set the before events days from 30 days to 60 days
df_len5_nonnorm<-sep_siteyears_data(df_norm_all,do_vars,df.sep20,5,60,0,10,FALSE)

#-------------------------------------------------------------------------
#(3)calculating some important variables like LUE...
#-------------------------------------------------------------------------
#--------------------------
#calculating more variables
#--------------------------
#1)LUE=GPP/fAPAR*ppfd
#unit-->GPP: umol m-2 s-1; ppdf-->umol m-2 s-1; fAPRA: unitless
#2)GRVI=(gcc-rcc)/c(gcc+rcc)
#3)albedo:alpha_SW<-SW_OUT/SW_IN; alpha_ppdf<-PPFD_IN/PPFD_OUT
#4)approximated fAPARchl=EVI*factor(factor = 1)
for(i in 1:length(df_len5_nonnorm)){
  #
  df_proc<-df_len5_nonnorm[[i]]
  df_proc$LUE<-df_proc$gpp_obs/c(df_proc$fapar_itpl*df_proc$PPFD_IN_fullday_mean_fluxnet2015)
  # df_proc$GRVI<-c(df_proc$gcc_90-df_proc$rcc_90)/c(df_proc$gcc_90+df_proc$rcc_90)
  df_proc$alpha_SW<-df_proc$SW_OUT_fullday_mean_fluxnet2015/df_proc$SW_IN_fullday_mean_fluxnet2015
  df_proc$alpha_PPFD<-df_proc$PPFD_OUT_fullday_mean_fluxnet2015/df_proc$PPFD_IN_fullday_mean_fluxnet2015
  df_proc$fAPAR_chl<-df_proc$evi*1
  #assign value back:
  df_len5_nonnorm[[i]]<-df_proc
}

#-------------------------------------------------------------------------
#(4)going to compare the "event" and "non-event" site with boxplot
#-------------------------------------------------------------------------
#--------------
#first to classify the "df_dday" and "df_noevent_dday" into different N classes
#--------------
library(plyr)  
library(ggplot2)
library(cowplot)
library(grid)
library(ggpubr)
Classify_Nbins_andPlot<-function(df,class_dday_range,class_var,N,do_manual_class,PFT_name){
  # df<-df_len5_nonnorm
  # class_dday_range<-c(-60,70)
  # class_var<-"SW_IN_midday_mean_fluxnet2015"
  # N<-10
  # do_manual_class<-FALSE
  # PFT_name<-"all PFTs"
  #----------------
  #do classification:
  #----------------
  df.event<-df$df_dday
  df.noevent<-df$df_noevent_dday
  #merge df.event and df.noevent
  df.event$year_flag<-rep("GPP overestimated sites",nrow(df.event))
  df.noevent$year_flag<-rep("GPP non-overestimated sites",nrow(df.noevent))
  #
  df.all<-rbind(df.event,df.noevent)
  df.all$year_flag<-factor(df.all$year_flag,levels = c("GPP overestimated sites","GPP non-overestimated sites"))
  #only select the data dday between -60 to 70:
  df.all<-df.all[df.all$dday>=class_dday_range[1] & df.all$dday<=class_dday_range[2],]
  #only targeting on one specific PFT
  if(PFT_name!="All PFTs"){
    df.all<-df.all %>% filter(classid==PFT_name)
  }
  
  # Bins for different variables
  xmin<-round(min(df.all[,class_var],na.rm = T),10)
  xmax<-round(max(df.all[,class_var],na.rm = T),10)         
  bins  <- seq( from=xmin, 
                to=xmax, by=(xmax-xmin)/N )
  
  #this is range setting for Tmin not for Rg
  if(do_manual_class==TRUE){
    xmin<-c(-40)
    xmax<-30
    bins  <- seq( from=xmin, 
                  to=xmax, by=(xmax-xmin)/N )
  }
  
  #
  df.all<-df.all %>% mutate(inbin=cut( as.numeric(df.all[,class_var]), breaks = bins ))
  #
  return(df.all)
}

##classify the df with the Rg
# df_bins<-Classify_Nbins_andPlot(df_len5_nonnorm,c(-60,70),"SW_IN_midday_mean_fluxnet2015",10,FALSE,"all")

#-----------------
#start plotting
#----------------
comp_boxplot<-function(df,comp_yvar,do_legend,end_xylab,PFT_name){
  # df<-df_bins
  # comp_yvar<-c("temp_min_fluxnet2015")
  # do_legend<-FALSE
  # end_xylab<-c("SW_IN midday mean (W m-2)","")
  # PFT_name<-"All PFTs"
  #------------------------------------
  #plotting 
  #------------------------------------
  ##
  df.tidy<-df[,c("sitename","date",comp_yvar,"inbin","year_flag")]
  names(df.tidy)<-c("sitename","date","comp_vary","inbin","year_flag")
  #remove the rows when inbin ==NA
  df.tidy<-df.tidy[!is.na(df.tidy$inbin),]
  ##
  p_plot_main<-ggplot(df.tidy,aes(x = inbin, y = comp_vary,fill=year_flag))+
    geom_boxplot()+
    # annotate("rect",xmin=0,xmax=70,ymin = -Inf,ymax = Inf,alpha=0.2)+  #
    scale_fill_manual("",values = c("GPP overestimated sites"=adjustcolor("tomato",1),
                    "GPP non-overestimated sites"=adjustcolor("green4",1)))+
    theme_classic()+
    theme(legend.position = c(0.2,0.95),legend.background = element_blank(),
          legend.text = element_text(size=20),
          axis.title = element_text(size = 30),
          axis.text = element_text(size = 24),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    rotate_x_text(30)+
    ylim(-45,25)+
    xlab(end_xylab[1])+
    ylab(end_xylab[2])
    
  if(comp_yvar=="alpha_SW"|comp_yvar=="alpha_PPFD"){
    p_plot_main<-p_plot_main+
      ylim(0,1)
  }
  #legend
  if(do_legend==FALSE){
    p_plot_main<-p_plot_main+
      theme(legend.position = "none")
  }
  #add the PFT name:
  N_x<-length(unique(df.tidy$inbin))
  p_plot_main<-p_plot_main+
    annotate(geom = "text",x=N_x,y=-30,label=PFT_name,col="blue",size=10)
  
   #
  return(p_plot_main)
}
##
# df<-df_bins
# comp_yvar<-c("SW_IN_fullday_mean_fluxnet2015")
# do_legend<-FALSE
# end_xylab<-c("Minium Ta (degreeC)","SW_IN full-dday mean (W m-2)")
# comp_boxplot(df_bins,c("SW_IN_fullday_mean_fluxnet2015"),TRUE,c("Minium Ta (degreeC)","SW_IN full-dday mean (W m-2)"))
#-------------------------------------------------------------------------
#(5) officially to compare the vars in two different group ("event" and "non_event")
#-------------------------------------------------------------------------
##mainly compare between the SW and Ta
#--------------
#classify based on the tmin(x),and plot others
#-------------
########
#A1.For all PFTs
########

#A1a).classification:
df_bins_SW_IN_all<-Classify_Nbins_andPlot(df_len5_nonnorm,c(-60,70),"SW_IN_midday_mean_fluxnet2015",8,FALSE,"All PFTs")
#change "ppfd_fluxnet2015" to "PPFD_IN_fullday_mean_fluxnet2015"/PPFD_IN_midday_mean_fluxnet2015
df_bins_ppfd_all<-Classify_Nbins_andPlot(df_len5_nonnorm,c(-60,80),"PPFD_IN_fullday_mean_fluxnet2015",8,FALSE,"All PFTs")

#A1b).plotting--->update in March,2022
#----------------
##For SW
#SW_midday_mean and tmin
p_SW_midday_mean_tmin_all<-comp_boxplot(df_bins_SW_IN_all,c("temp_min_fluxnet2015"),FALSE,
                                   c("SW_IN midday mean (W m-2)","Minimum Ta (°C)"),"All PFTs")
#SW_midday_mean and max
p_SW_midday_mean_tmax_all<-comp_boxplot(df_bins_SW_IN_all,c("temp_max_fluxnet2015"),FALSE,
                                    c("SW_IN middday mean (W m-2)","Maximum Ta (°C)"),"All PFTs")
#SW_midday_mean and tmean
p_SW_midday_mean_tmean_all<-comp_boxplot(df_bins_SW_IN_all,c("temp_day_fluxnet2015"),FALSE,
                                   c("SW_IN middday mean (W m-2)","mean Ta (°C)"),"All PFTs")
#tmin and alpha_SW
p_SW_midday_mean_alpha_SW_all<-comp_boxplot(df_bins_SW_IN_all,c("alpha_SW"),FALSE,c("SW_IN middday mean (W m-2)","alpha_SW"),"All PFTs")
# tmin and alpha_PPFD
# p_SW_midday_mean_alpha_PPFD<-comp_boxplot(df_bins_SW_IN,c("alpha_PPFD"),FALSE,c("SW_IN middday mean (W m-2)","alpha_PPFD"))
#-------------
##For ppfd
p_ppfd_midday_mean_tmin_all<-comp_boxplot(df_bins_ppfd_all,c("temp_min_fluxnet2015"),FALSE,
                                    c("PAR midday mean (umol m-2 s-1)","Minimum Ta (°C)"),"All PFTs")
p_ppfd_midday_alpha_PPFD_all<-comp_boxplot(df_bins_ppfd_all,c("alpha_PPFD"),FALSE,
                                       c("PAR middday mean (W m-2)","alpha_PPFD"),"All PFTs")
#------------------------------
#A1c).change the labels 
p_SW_midday_mean_tmin_all<-p_SW_midday_mean_tmin_all+
  xlab("")+
  ylab(expression("T"[min]*" (°C)"))
  # xlab(expression("SW"[IN]*" midday mean (W m"^-2*")"))
  
p_ppfd_midday_mean_tmin_all<-p_ppfd_midday_mean_tmin_all+
  xlab("")+
  ylab(expression("T"[min]*" (°C)"))
  # xlab(expression("PAR (umol m"^-2*"s"^-1*")"))

########
#A2.For different PFTs:DBF, MF, and ENF
########
df_bins_SW_IN_DBF<-Classify_Nbins_andPlot(df_len5_nonnorm,c(-60,70),"SW_IN_midday_mean_fluxnet2015",8,FALSE,"DBF")
df_bins_SW_IN_MF<-Classify_Nbins_andPlot(df_len5_nonnorm,c(-60,70),"SW_IN_midday_mean_fluxnet2015",8,FALSE,"MF")
df_bins_SW_IN_ENF<-Classify_Nbins_andPlot(df_len5_nonnorm,c(-60,70),"SW_IN_midday_mean_fluxnet2015",8,FALSE,"ENF")
#
df_bins_ppfd_DBF<-Classify_Nbins_andPlot(df_len5_nonnorm,c(-60,80),"PPFD_IN_fullday_mean_fluxnet2015",8,FALSE,"DBF")
df_bins_ppfd_MF<-Classify_Nbins_andPlot(df_len5_nonnorm,c(-60,80),"PPFD_IN_fullday_mean_fluxnet2015",8,FALSE,"MF")
df_bins_ppfd_ENF<-Classify_Nbins_andPlot(df_len5_nonnorm,c(-60,80),"PPFD_IN_fullday_mean_fluxnet2015",8,FALSE,"ENF")

##For SW
p_SW_midday_mean_tmin_DBF<-comp_boxplot(df_bins_SW_IN_DBF,c("temp_min_fluxnet2015"),FALSE,
                                        c("SW_IN midday mean (W m-2)","Minimum Ta (°C)"),"DBF")
p_SW_midday_mean_tmin_MF<-comp_boxplot(df_bins_SW_IN_MF,c("temp_min_fluxnet2015"),FALSE,
                                        c("SW_IN midday mean (W m-2)","Minimum Ta (°C)"),"MF")
p_SW_midday_mean_tmin_ENF<-comp_boxplot(df_bins_SW_IN_ENF,c("temp_min_fluxnet2015"),FALSE,
                                        c("SW_IN midday mean (W m-2)","Minimum Ta (°C)"),"ENF")
##For ppfd
p_ppfd_midday_mean_tmin_DBF<-comp_boxplot(df_bins_ppfd_DBF,c("temp_min_fluxnet2015"),FALSE,
                                          c("PAR midday mean (umol m-2 s-1)","Minimum Ta (°C)"),"DBF")
p_ppfd_midday_mean_tmin_MF<-comp_boxplot(df_bins_ppfd_MF,c("temp_min_fluxnet2015"),FALSE,
                                          c("PAR midday mean (umol m-2 s-1)","Minimum Ta (°C)"),"MF")
p_ppfd_midday_mean_tmin_ENF<-comp_boxplot(df_bins_ppfd_ENF,c("temp_min_fluxnet2015"),FALSE,
                                          c("PAR midday mean (umol m-2 s-1)","Minimum Ta (°C)"),"ENF")
#for albedo:
p_SW_midday_mean_alpha_SW_DBF<-comp_boxplot(df_bins_SW_IN_all,c("alpha_SW"),FALSE,c("SW_IN middday mean (W m-2)","alpha_SW"),"DBF")
p_ppfd_midday_alpha_PPFD_DBF<-comp_boxplot(df_bins_ppfd_all,c("alpha_PPFD"),FALSE,
                                           c("PAR middday mean (W m-2)","alpha_PPFD"),"DBF")
p_SW_midday_mean_alpha_SW_MF<-comp_boxplot(df_bins_SW_IN_all,c("alpha_SW"),FALSE,c("SW_IN middday mean (W m-2)","alpha_SW"),"MF")
p_ppfd_midday_alpha_PPFD_MF<-comp_boxplot(df_bins_ppfd_all,c("alpha_PPFD"),FALSE,
                                           c("PAR middday mean (W m-2)","alpha_PPFD"),"MF")
p_SW_midday_mean_alpha_SW_ENF<-comp_boxplot(df_bins_SW_IN_all,c("alpha_SW"),FALSE,c("SW_IN middday mean (W m-2)","alpha_SW"),"ENF")
p_ppfd_midday_alpha_PPFD_ENF<-comp_boxplot(df_bins_ppfd_all,c("alpha_PPFD"),FALSE,
                                           c("PAR middday mean (W m-2)","alpha_PPFD"),"ENF")

#change the labels 
#
p_SW_midday_mean_tmin_DBF<-p_SW_midday_mean_tmin_DBF+
  # xlab(expression("SW"[IN]*" midday mean (W m"^-2*")"))+
  xlab("")+
  ylab("")
p_SW_midday_mean_tmin_MF<-p_SW_midday_mean_tmin_MF+
  xlab(expression("SW"[IN]*" midday mean (W m"^-2*")"))+
  ylab(expression("T"[min]*" (°C)"))
p_SW_midday_mean_tmin_ENF<-p_SW_midday_mean_tmin_ENF+
  ylab("")+
  xlab(expression("SW"[IN]*" midday mean (W m"^-2*")"))
#
p_ppfd_midday_mean_tmin_DBF<-p_ppfd_midday_mean_tmin_DBF+
  xlab("")+
  ylab("")
  # xlab(expression("PAR (umol m"^-2*"s"^-1*")"))
p_ppfd_midday_mean_tmin_MF<-p_ppfd_midday_mean_tmin_MF+
  xlab(expression("PAR midday mean (umol m"^-2*"s"^-1*")"))+
  ylab(expression("T"[min]*" (°C)"))
p_ppfd_midday_mean_tmin_ENF<-p_ppfd_midday_mean_tmin_ENF+
  ylab("")+
  xlab(expression("PAR midday mean (umol m"^-2*"s"^-1*")"))
#

#B. merge the plots and save:
save.path<-"./manuscript/figures/"
p_SW_IN_Tmin<-plot_grid(p_SW_midday_mean_tmin_all,p_SW_midday_mean_tmin_DBF,
                        p_SW_midday_mean_tmin_MF,p_SW_midday_mean_tmin_ENF,
                        nrow = 2,ncol=2,labels = "auto",label_size = 20)
ggsave(paste0(save.path,"Figure3_boxplot_SW_Tmin.png"),p_SW_IN_Tmin,width = 25,height = 15)
#
p_ppfd_Tmin<-plot_grid(p_ppfd_midday_mean_tmin_all,p_ppfd_midday_mean_tmin_DBF,
                        p_ppfd_midday_mean_tmin_MF,p_ppfd_midday_mean_tmin_ENF,
                        nrow = 2,ncol=2,labels = "auto",label_size = 20)
ggsave(paste0(save.path,"FigureS_boxplot_ppfd_Tmin.png"),p_ppfd_Tmin,width = 25,height = 15)
#test albedo
p_albedo1<-plot_grid(p_SW_midday_mean_alpha_SW_all,p_SW_midday_mean_alpha_SW_DBF,
                    p_SW_midday_mean_alpha_SW_MF,p_SW_midday_mean_alpha_SW_ENF,
                       nrow = 2,ncol=2)
p_albedo2<-plot_grid(p_ppfd_midday_alpha_PPFD_all,p_ppfd_midday_alpha_PPFD_DBF,
                     p_ppfd_midday_alpha_PPFD_MF,p_ppfd_midday_alpha_PPFD_ENF,
                    nrow = 2,ncol=2)

