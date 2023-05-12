#############################################################################
#Aim:compare the variables in "event" and "non-event" site years after aligning the data-->specifically for Temp
library(dplyr)
library(lubridate)
library(tidyverse)
#-------------------------------------------------------------------------
#(1)load the data that includes the "is_event" information
#-------------------------------------------------------------------------
#---------------------
#A.load the event_length data
#---------------------
load.path<-"./data/event_length/"
load(paste0(load.path,"df_events_length.RDA"))
#
tt<-df_events_all %>%
  group_by(sitename) %>%
  dplyr::summarise(Years=range(Year))
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
df.sep30<-sep_siteyears(df_events_all,"Over_days_length",30)
#---------------------
#B.load the data 
#---------------------
##load the data
load.path<-"./data/data_used/"
#from new method:
load(paste0(load.path,"ddf_labeled_norm_trs_newmethod_all_overestimation_Fluxnet2015_sites.RDA"))
df_all_sites<-ddf_labeled;rm(ddf_labeled) 
df_norm_all<-df_all_sites

#load the site infos:(including the PFTs and Cimate info)
load(paste0("./data-raw/raw_data/sites_info/","Pre_selected_sites_info.RDA"))
sites.info<-df_sites_sel;rm(df_sites_sel)
#
df_norm_all<-left_join(df_norm_all,sites.info)
df_norm_all<-df_norm_all %>%
  mutate(Clim_PFTs=paste0(koeppen_code,"-",classid))
#
##-->update in May,2023-->also add the snow data from modis
snow.path<-"./data-raw/raw_data/snow_data/"
load(paste0(snow.path,"snow_MOD.RDA"))
snow_MOD_sel$date<-as.POSIXct(snow_MOD_sel$date)
#only select the sites are used for the analysis:
sel_sites<-unique(df_events_all$sitename)
snow_MOD_sel<-snow_MOD_sel %>%
  filter(sitename %in% sel_sites)
#---merge df_norm_all with snow data:
#!! the merge does not really mrege the snowval-->do not know exact reason!!
df_norm_all<-left_join(df_norm_all,snow_MOD_sel[,c("sitename","date","snowval","snowval_QA")])
#test:
# df_norm_all%>%
# group_by(sitename)%>%
# ggplot(aes(x=date,y=snowval))+
#   geom_point()+
#   xlim(0,100)+
#   facet_wrap(~sitename)

#test!
df.test<-left_join(df.sep30$noevent_siteyears,sites.info)
df.test1<-df.test %>%
  filter(classid=="DBF")
df.test<-left_join(df.sep30$event_siteyears,sites.info)
df.test2<-df.test %>%
  filter(classid=="DBF")

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

sep_siteyears_data<-function(df.data,dovars,df.sep,leng_threshold,before,after,nbins,do_norm){
  # df.data<-df_norm_all
  # dovars<-c("gpp_obs","fapar_itpl","gcc_90")
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
#do_vars-->the variables that are going to be processed:
names(df_norm_all)
do_vars<-c("gpp_obs","fapar_itpl","fapar_spl",paste0(c("ppfd","PPFD_IN_fullday_mean","temp_day","temp_min","temp_max",
            "vpd_day","prec","patm","SW_IN","ws",paste0("TS_",1:7),paste0("SWC_",1:5)),"_fluxnet2015"),
           "gcc_90","rcc_90","snowval")
#set the before events days from 30 days to 60 days
df_len5_nonnorm<-sep_siteyears_data(df_norm_all,do_vars,df.sep30,5,60,0,10,FALSE)
#calculate the site-years for each group (event and non-event sites):
#
event_siteyears<-df_len5_nonnorm$df_dday %>%
  select(sitename,Year)%>%
  mutate(site_years=paste0(sitename,Year))
length(unique(event_siteyears$site_years))
#
noevent_siteyears<-df_len5_nonnorm$df_noevent_dday %>%
  select(sitename,Year)%>%
  mutate(site_years=paste0(sitename,Year))
length(unique(noevent_siteyears$site_years))

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
  df_proc$GRVI<-c(df_proc$gcc_90-df_proc$rcc_90)/c(df_proc$gcc_90+df_proc$rcc_90)
  df_proc$alpha_SW<-df_proc$SW_OUT_fullday_mean_fluxnet2015/df_proc$SW_IN_fullday_mean_fluxnet2015
  df_proc$alpha_PPFD<-df_proc$PPFD_OUT_fullday_mean_fluxnet2015/df_proc$PPFD_IN_fullday_mean_fluxnet2015
  df_proc$fAPAR_chl<-df_proc$evi*1
  ##adding the mean values of soil mositure of first 3 layers:
  df_proc$TS_top3_fluxnet2015<-rowMeans(df_proc[,c("TS_1_fluxnet2015",
          "TS_2_fluxnet2015","TS_3_fluxnet2015")])
  #assign value back:
  df_len5_nonnorm[[i]]<-df_proc
}
#check the data aviablity for each variable:
library(visdat)
pos_TS<-grep("TS_",names(df_len5_nonnorm$df_dday))
pos_SWC<-grep("SWC_",names(df_len5_nonnorm$df_dday))
# visdat::vis_miss(df_len5_nonnorm$df_dday[,pos_TS], warn_large_data = FALSE)
# visdat::vis_miss(df_len5_nonnorm$df_dday[,pos_SWC], warn_large_data = FALSE)
#----calculate the mean T between "overestimated site" and "non-overestimated site"--
#only using the data bebtween Jan and June:
##the code is in "photocold_manuscript/test/"

#-------------------------------------------------------------------------
#(4)going to compare the "event" and "non-event" site
#-------------------------------------------------------------------------
  #-------------------------------------------------------------------------
  #start plotting
  #-------------------------------------------------------------------------
  library(plyr)  
  library(ggplot2)
  library(cowplot)
  library(grid)
  
#function used to compare the "event" and "non_event" sites using geom_ribbon 
#source the comparsion test written by me:
fun.path<-"D:/Github/photocold_lyp/R/Step2_identify_events/Functions/functions_from_YP/"
source(paste0(fun.path,"Difference_test_for_2Classes.R"))

#
# plot(df_len5_nonnorm$df_dday$date,df_len5_nonnorm$df_dday$TS_1_fluxnet2015,
#      xlab = "Date",ylab="SY_PSB Tsoil(degC)")
plot_2groups<-function(df,comp_var,var_unit,do_norm,do_legend){
    # df<-df_len5_nonnorm
    # comp_var<-"gpp_obs"
    # var_unit<-"degC"
    # do_norm<-FALSE
    # do_legend<-TRUE

    #
    df.event<-df$df_dday
    df.noevent<-df$df_noevent_dday
    #------------------------------------
    #remove the sites that have the point numbers smaller than 10
    #------------------------------------
    N_value<-ddply(df.noevent,.(sitename),summarise,N=length(date))
    rm_sites<-subset(N_value,N<10)
    if(nrow(rm_sites)>=1){
      for(i in nrow(rm_sites)){
        rm_site_name<-rm_sites[i]
        df.noevent<-df.noevent[df.noevent$sitename!=as.character(rm_site_name),]
      }  
    }
    df.noevent<-df.noevent
    #--------------------------------
    #selected the most relevant vars
    #for the comparison,select the original variable if "do_norm"==FALSE,otherwise "do_norm"==TRUE
    #--------------------------------
    if(do_norm==FALSE){
      #update in May, 2023-->add classid variable
      df.event_sel<-df.event[,c("sitename","Year","date","dday","classid",comp_var)]
      names(df.event_sel)<-c("sitename","Year","date","dday","PFT","comp_var")
      df.nonevent_sel<-df.noevent[,c("sitename","Year","date","dday","classid",comp_var)]
      names(df.nonevent_sel)<-c("sitename","Year","date","dday","PFT","comp_var")
    }
    if(do_norm==TRUE){
      comp_var_norm<-paste0("ds",comp_var)
      df.event_sel<-df.event[,c("sitename","Year","date","dday","classid",comp_var_norm)]
      names(df.event_sel)<-c("sitename","Year","date","dday","PFT","comp_var")
      #
      df.nonevent_sel<-df.noevent[,c("sitename","Year","date","dday","classid",comp_var_norm)]
      names(df.nonevent_sel)<-c("sitename","Year","date","dday","PFT","comp_var")
      #also make var_unit()==""
      var_unit<-c()
    }
    #--------------------------------------
    #plotting
    #--------------------------------------
    # Create a text
    grob_event <- grobTree(textGrob("Event site-year", x=0.6,  y=0.9, hjust=0,
                                    gp=gpar(col="red", fontsize=18, fontface="italic")))
    # grob_event_name <- grobTree(textGrob(df_name, x=0.8,  y=0.1, hjust=0,
    #                                      gp=gpar(col="red", fontsize=18, fontface="italic")))
    grob_nonevent <- grobTree(textGrob("Non-Event site-year", x=0.6,  y=0.9, hjust=0,
                                       gp=gpar(col="green4", fontsize=18, fontface="italic")))
    # grob_nonevent_name<-grobTree(textGrob(df_name, x=0.8,  y=0.1, hjust=0,
    #                                       gp=gpar(col="red", fontsize=18, fontface="italic")))
    #y axis range:
    ymin<-min(range(df.event_sel$comp_var,na.rm = T),range(df.nonevent_sel$comp_var,na.rm = T))
    ymax<-max(range(df.event_sel$comp_var,na.rm = T),range(df.nonevent_sel$comp_var,na.rm = T))
    #x axis range
    # x_range_event<-range(df.event_sel$doy)
    # x_range_nonevent<-range(df.nonevent_sel$doy)
    ###start to make the quantiile plot:
    ##2022,Aug-->change the "GPP overestimated sites" to "SY_PSB", the other to "SY_ASB"
    # df.event_sel$flag<-rep("GPP overestimated sites",nrow(df.event_sel))
    # df.nonevent_sel$flag<-rep("GPP non-overestimated sites",nrow(df.nonevent_sel))
    df.event_sel$flag<-rep("SY_PSB",nrow(df.event_sel))
    df.nonevent_sel$flag<-rep("SY_ASB",nrow(df.nonevent_sel))
    #!one thing need to pay attention-->use q90 to limit the dday range in sites as the sites diff a lot
    #to ensure the results do not impact by one specific site
    dday_range_event<-ddply(df.event_sel,.(sitename),summarize,min_dday=min(dday),max_dday=max(dday))
    dday_range_nonevent<-ddply(df.nonevent_sel,.(sitename),summarize,min_dday=min(dday),max_dday=max(dday))
    sel_dday_event<-floor(quantile(dday_range_event$max_dday,0.75))
    sel_dday_nonevent<-floor(quantile(dday_range_nonevent$max_dday,0.75))
    #
    df.event_sel<-df.event_sel[df.event_sel$dday<=sel_dday_event,]
    df.nonevent_sel<-df.nonevent_sel[df.nonevent_sel$dday<=sel_dday_nonevent,]
    ##merge and classify 
    df.all<-rbind(df.event_sel,df.nonevent_sel)
    df.all$flag<-factor(df.all$flag,levels = c("SY_PSB","SY_ASB"))
    ##
    #for min temperature-->first find the mininum temperature for each site then calculate the quantile
    if(comp_var=="temp_min_fluxnet2015"){
      df.all<-df.all %>%
        group_by(sitename,dday,flag) %>%
        dplyr::summarise(comp_var=min(comp_var,na.rm = T))
    }
    
    #merge the different sites in "event" and "non-event" sites->calculate the quantiles at the same time
    df.all_q<-ddply(df.all,.(flag,dday),summarize,q10=quantile(comp_var,0.1,na.rm = T),q25=quantile(comp_var,0.25,na.rm = T),
          q50=quantile(comp_var,0.5,na.rm = T),q75=quantile(comp_var,0.75,na.rm = T),q90=quantile(comp_var,0.9,na.rm = T))
    df.all_q_PFTs<-ddply(df.all,.(flag,PFT,dday),summarize,q10=quantile(comp_var,0.1,na.rm = T),q25=quantile(comp_var,0.25,na.rm = T),
                         q50=quantile(comp_var,0.5,na.rm = T),q75=quantile(comp_var,0.75,na.rm = T),q90=quantile(comp_var,0.9,na.rm = T),
                         mean=mean(comp_var,na.rm=T),sd=sd(comp_var,na.rm = T))
    #-----------------------
    #add a additional test-->to compare if the vars are significant different(q50) in two categories
    #----------------------
    # compare_diff<-function(df_comp,sel_perc,flag){
    #   # df_comp<-df.all
    #   # sel_perc<-1             #which period (percentage of the minimum days) used to compare between vars in two category
    #   # flag<-"all"
    # 
    #   df_comp_event<-df_comp[df_comp$flag=="SY_PSB",]
    #   df_comp_nonevent<-df_comp[df_comp$flag=="SY_ASB",]
    #   #------
    #   #select the data
    #   #------
    #   range_dday_event<-range(df_comp_event$dday)
    #   range_dday_nonevent<-range(df_comp_nonevent$dday)
    #   range_sel<-c(min(range_dday_event[1],range_dday_nonevent[1]),
    #                floor(sel_perc*min(range_dday_event[2],range_dday_nonevent[2])))
    #   #
    #   df_comp_event_sel<-subset(df_comp_event,dday>=range_sel[1]&dday<=range_sel[2])
    #   df_comp_nonevent_sel<-subset(df_comp_nonevent,dday>=range_sel[1]&dday<=range_sel[2])
    #   #additional-->for before events(b0):
    #   df_comp_event_b0<-subset(df_comp_event,dday>=range_sel[1]&dday<=0)
    #   df_comp_nonevent_b0<-subset(df_comp_nonevent,dday>=range_sel[1]&dday<=0)
    #   #make the statistical comparision:
    #   if(flag=="q50"){
    #     stat.test_1<-difftest_function(df_comp_event_sel$q50,df_comp_nonevent_sel$q50)
    #     stat.test_2<-difftest_function(df_comp_event_b0$q50,df_comp_nonevent_b0$q50)
    #     stat.test<-cbind(stat.test_1,stat.test_2)
    #     names(stat.test)<-c("All selecte period","Before events")
    #   }
    #   if(flag=="all"){
    #     stat.test_1<-difftest_function(df_comp_event_sel$comp_var,df_comp_nonevent_sel$comp_var)
    #     stat.test_2<-difftest_function(df_comp_event_b0$comp_var,df_comp_nonevent_b0$comp_var)
    #     stat.test<-cbind(stat.test_1,stat.test_2)
    #     names(stat.test)<-c("All selecte period","Before events")
    #   }
    #   return(stat.test)
    # }
    # #
    # stat.alldata_for2<-compare_diff(df.all,1,"all")
    # stat.q50_for2<-compare_diff(df.all_q,1,"q50")
    #making the quantile plot using ribbon function:
    p_plot<-ggplot(df.all_q)+
      # annotate("rect",xmin=0,xmax=max(df.all_q$dday),ymin = -Inf,ymax = Inf,alpha=0.2)+
      #some changes here
      annotate("rect",xmin=0,xmax=70,ymin = -Inf,ymax = Inf,alpha=0.2)+
      geom_line(aes(x=dday,y=q50,col=flag),size=1.05)+
      scale_color_manual("",values = c("SY_PSB"="red","SY_ASB"="blue"),
                         labels=c(expression(SY[DSPR]),expression(SY[0])))+ ##add subscript in the legend
      # geom_ribbon(aes(x=dday,ymin=q10,ymax=q90,fill=flag),alpha=0.15)+
      geom_ribbon(aes(x=dday,ymin=q25,ymax=q75,fill=flag),alpha=0.4)+
      scale_fill_manual("",values = c("SY_PSB"="red","SY_ASB"="dodgerblue"),
                        labels=c(expression(SY[DSPR]),expression(SY[0])))+
      ylab(paste0(comp_var," ",var_unit))+
      xlab("rday")+
      theme_classic()+
      theme(legend.position = c(0.3,0.9),legend.background = element_blank(),
            legend.key.size = unit(2, 'lines'),
            legend.text = element_text(size=20),
            axis.title = element_text(size=24),
            axis.text = element_text(size = 20))+
      theme(legend.text.align = 0)+  #align the legend (all the letter start at the same positoin)
      xlim(-60,70)  #add x range in 2021-09-25
      
    #legend
    if(do_legend==FALSE){
    p_plot<-p_plot+
      theme(legend.position = "none")
    }
    # print(p_plot)
    #returun object
    out<-c()
    out$data_PFTs<-df.all_q_PFTs #add in May, 2023-->to separate the data in different PFTs
    # out$stats_q50<-stat.q50_for2
    # out$stats_alldata<-stat.alldata_for2
    out$plot<-p_plot
    
    return(out)
  }

#-------------------------------------------------------------------------
#(5) officially to compare the vars in two different group ("event" and "non_event")
#-------------------------------------------------------------------------
  ###potentially variables need to be checked:
  #fAPAR-itpl;fAPRA_spl
  #ppfd
  #temp_day
  #vpd_day
  #prec
  #patm
  #gpp_obs,gpp_mod_full,gpp_obs_norm,gpp_mod_norm,gpp_res

#2021-09-25:check the temp_min,temp_max,temp_day, and LUE
#update in 2021-11-16:update the gpp_obs, gpp_biaes
#--------------
#I.GPP and LUE, PAR, and fAPAR
#-------------
#gpp_obs
p_gpp_obs_len5_b60<-plot_2groups(df_len5_nonnorm,"gpp_obs","(umol m-2 s-1)",do_norm = FALSE,do_legend = TRUE)
#gpp biaes
p_gpp_biaes_len5_b60<-plot_2groups(df_len5_nonnorm,"gpp_res","(umol m-2 s-1)",do_norm = FALSE,FALSE)
#LUE
p_LUE_len5_b60<-plot_2groups(df_len5_nonnorm,"LUE","",do_norm = FALSE,TRUE)
#for ppfd: change the "ppfd_fluxnet2015" to "PPFD_IN_fullday_mean_fluxnet2015"
# p_ppfd_len5_b60<-plot_2groups(df_len5_nonnorm,"ppfd_fluxnet2015","(u mol m-2 s-1)",do_norm = FALSE,do_legend = TRUE)
p_ppfd_len5_b60<-plot_2groups(df_len5_nonnorm,"PPFD_IN_fullday_mean_fluxnet2015","(u mol m-2 s-1)",do_norm = FALSE,do_legend = FALSE)
#fapar_spl and fapar_itpl
p_fapar_itpl_len5_b60<-plot_2groups(df_len5_nonnorm,"fapar_itpl","",do_norm = FALSE,do_legend = FALSE)
#some modifying in the plot:
p_gpp_obs_len5_b60$plot<-p_gpp_obs_len5_b60$plot+
  xlab("")+
  ylab(expression("GPP"[obs]*" (g "*"m"^-2*" d"^-1*")"))
p_gpp_biaes_len5_b60$plot<-p_gpp_biaes_len5_b60$plot+
  xlab("")+
  ylab(expression("GPP bias"*" (g "*"m"^-2*" d"^-1*")"))
p_ppfd_len5_b60$plot<-p_ppfd_len5_b60$plot+
  # xlab("dday")+
  ylab(expression("PAR"*" (u mol "*"m"^-2*" s"^-1*")"))
p_fapar_itpl_len5_b60$plot<-p_fapar_itpl_len5_b60$plot+
  xlab("")+
  ylim(0.1,1)+
  ylab("fAPAR")
#------------------
#II. fAPAR in different PFTs
#------------------
p_fapar_itpl_len5_b60<-plot_2groups(df_len5_nonnorm,"fapar_itpl","",do_norm = FALSE,do_legend = FALSE)
#data for different PFTs--->update in May, 2023:
df.all_q_PFTs<-p_fapar_itpl_len5_b60$data_PFTs
#making the quantile plot using ribbon function:
p_fAPR_PFT<-ggplot(df.all_q_PFTs)+
  # annotate("rect",xmin=0,xmax=max(df.all_q$dday),ymin = -Inf,ymax = Inf,alpha=0.2)+
  #some changes here
  annotate("rect",xmin=0,xmax=70,ymin = -Inf,ymax = Inf,alpha=0.2)+
  geom_line(aes(x=dday,y=q50,col=flag),size=1.05)+
  scale_color_manual("",values = c("SY_PSB"="red","SY_ASB"="blue"),
                     labels=c(expression(SY[DSPR]),expression(SY[0])))+ ##add subscript in the legend
  # geom_ribbon(aes(x=dday,ymin=q10,ymax=q90,fill=flag),alpha=0.15)+
  geom_ribbon(aes(x=dday,ymin=q25,ymax=q75,fill=flag),alpha=0.4)+
  scale_fill_manual("",values = c("SY_PSB"="red","SY_ASB"="dodgerblue"),
                    labels=c(expression(SY[DSPR]),expression(SY[0])))+
  ylab("fAPAR")+
  xlab("rday")+
  facet_wrap( ~PFT) +
  theme_classic()+
  theme(legend.position = c(0.1,0.9),legend.background = element_blank(),
        legend.key.size = unit(2, 'lines'),
        legend.text = element_text(size=18),
        axis.title = element_text(size=20),
        axis.text = element_text(size = 18),
        #increase the facet label font
        strip.text.x = element_text(size=20)
        )+
  theme(legend.text.align = 0)+  #align the legend (all the letter start at the same positoin)
  xlim(-60,70)  #add x range in 2021-09-25

##add the d,e,f indicate the panels:
dat_text<-data.frame(label=c("a","b","c"),classid=c("DBF","ENF","MF"),
                     x=rep(-58,3),y=rep(0.9,3))
# dat_text$parameter<-factor(dat_text$parameter,levels = c("tau","X0","Smax"))
p_fAPR_PFT<-p_fAPR_PFT+geom_text(
  data=dat_text,
  size=6,col="black",
  mapping = aes(x=x,y=y,label=c(rep("a",3),rep("b",3),rep("c",3))) #interesting-->need to add 3
)+
theme(axis.title.x = element_blank())
#------------------
#III. snow in different PFTs
#------------------
#contrast for snow infomation for "event(GPP overestimation)" site and "nonevent" site
snow_MOD_sel$snowval[!is.na(snow_MOD_sel$snowval)]
##
event_siteyear<-df.sep30$event_siteyears%>%
  mutate(siteyear=paste0(sitename,"-",Year))
nonevent_siteyear<-df.sep30$noevent_siteyears%>%
  mutate(siteyear=paste0(sitename,"-",Year))

#
snow_event_siteyear<-snow_MOD_sel%>%
  mutate(year=year(date),siteyear=paste0(sitename,"-",year))%>%
  filter(siteyear %in% unique(event_siteyear$siteyear))

snow_nonevent_siteyear<-snow_MOD_sel%>%
  mutate(year=year(date),siteyear=paste0(sitename,"-",year))%>%
  filter(siteyear %in% unique(nonevent_siteyear$siteyear))

#adding the site PFT:
snow_event_siteyear<-left_join(snow_event_siteyear,sites.info[,c("sitename","classid")])
snow_event_siteyear$flag<-rep("SY_PSB",nrow(snow_event_siteyear))
#
snow_nonevent_siteyear<-left_join(snow_nonevent_siteyear,sites.info[,c("sitename","classid")])
snow_nonevent_siteyear$flag<-rep("SY_ASB",nrow(snow_nonevent_siteyear))
#------------roughly adding the mean sos values:
df_events_all_sel<-df_events_all[,c("sitename","Year",'sos',"peak")]
names(df_events_all_sel)<-c("sitename","year",'sos',"peak")
snow_event_siteyear<-left_join(snow_event_siteyear,df_events_all_sel)
snow_nonevent_siteyear<-left_join(snow_nonevent_siteyear,df_events_all_sel)

#merge
df.all_snow<-rbind(snow_event_siteyear,snow_nonevent_siteyear)

#---plotting----
df.all_snow_mean_PFT<-df.all_snow %>%
  mutate(doy=yday(date))%>%
  group_by(classid,flag,doy)%>%
  dplyr::summarise(snow_frac=mean(snowval,na.rm=T),
                   snow_frac_sd=sd(snowval,na.rm=T),
                   sos_mean=round(mean(sos,na.rm=T),0),
                   peak_mean=round(mean(peak,na.rm=T),0))%>%
  mutate(rday=doy-sos_mean)
  

#data for different PFTs--->update in May, 2023:
#making the quantile plot using ribbon function:
p_snow_PFT<-df.all_snow_mean_PFT%>%
  group_by(classid)%>%
  ggplot()+
  annotate("rect",xmin=0,xmax=70,ymin = -Inf,ymax = Inf,alpha=0.2)+
  geom_point(aes(x=rday,y=snow_frac,col=flag),size=1.05)+
  scale_color_manual("",values = c("SY_PSB"="red","SY_ASB"="blue"),
                     labels=c(expression(SY[DSPR]),expression(SY[0])))+ ##add subscript in the legend
  # geom_ribbon(aes(x=dday,ymin=q10,ymax=q90,fill=flag),alpha=0.15)+
  geom_ribbon(aes(x=rday,ymin=snow_frac-snow_frac_sd,ymax=snow_frac+snow_frac_sd,fill=flag),alpha=0.4)+
  scale_fill_manual("",values = c("SY_PSB"="red","SY_ASB"="dodgerblue"),
                    labels=c(expression(SY[DSPR]),expression(SY[0])))+
  ylab("Snow fraction (%)")+
  xlab("rday")+
  facet_wrap( ~classid) +
  theme_classic()+
  theme(legend.position = c(0.1,0.9),legend.background = element_blank(),
        legend.key.size = unit(2, 'lines'),
        legend.text = element_text(size=18),
        axis.title = element_text(size=20),
        axis.text = element_text(size = 18),
        strip.text.x = element_text(size=20))+
  theme(legend.text.align = 0)+  #align the legend (all the letter start at the same positoin)
  xlim(-60,70)  #add x range 
#
##add the d,e,f indicate the panels:
dat_text<-data.frame(label=c("d","e","f"),classid=c("DBF","ENF","MF"),
                     x=rep(-58,3),y=rep(90,3))
# dat_text$parameter<-factor(dat_text$parameter,levels = c("tau","X0","Smax"))
p_snow_PFT<-p_snow_PFT+geom_text(
  data=dat_text,
  size=6,col="black",
  mapping = aes(x=x,y=y,label=label)
)+
theme(legend.position =  "none")

##merge fAPAR and snow fraction:
plot_merge<-cowplot::plot_grid(p_fAPR_PFT,p_snow_PFT,
          ncol = 1,labels = "auto",label_size = 20,align = "hv")
#save the plots:
save.path<-"./manuscript/figures/"
ggsave(paste0(save.path,"FigureS_fAPAR_snow.png"),plot_merge,
       width = 15,height = 10)


