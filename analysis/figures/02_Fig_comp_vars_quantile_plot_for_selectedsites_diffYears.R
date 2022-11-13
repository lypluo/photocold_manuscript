##---------------------------------------
#Aim: To each site, align different years with start of season and compare different
#years' environmental variables-->to see if if lower Tmin with more severe overestimation
##---------------------------------------

#-------------------------------------------------------------------------
##(1)selecting which several sites that are representative for the analysis:
#-------------------------------------------------------------------------
#dding the information of how many percentage of year with spring bias%
event.data.path<-"./data/event_length/"
load(file=paste0(event.data.path,"df_events_length.RDA"))
#first add a flag to indicate if the sites have overestimation or not;
df_events_all<-df_events_all %>%
  mutate(event_flag=rep(NA,nrow(df_events_all)),
         event_flag=ifelse(Over_days_length<20|is.na(Over_days_length),"no","yes"))
t1<-df_events_all %>%
  group_by(sitename)%>%
  dplyr::summarise(nyear=length(unique(Year)))
t2<-df_events_all%>%
  filter(event_flag=="yes")%>%
  group_by(sitename)%>%
  dplyr::summarise(nyear_event=length(unique(Year)),
                   event_day_mean=mean(Over_days_length),
                   event_day_mean=round(event_day_mean,2))
t_merge<-left_join(t1,t2)
event_merge<-t_merge%>%
  mutate(event_day_mean=ifelse(is.na(event_day_mean),0,event_day_mean))%>%
  mutate(event_perc=nyear_event/nyear)%>%
  mutate(event_perc=ifelse(is.na(event_perc),0,event_perc),
         event_perc=round(event_perc,2))
#-->select several sites for further analysis:
#selection crition:1)sites with high and intermediate gpp overestimation;
#2)different PFT;(3)different continent:North America; Eurasia
#3)as many years data as possible
#DBF(US-UMB;RU-Fyo);ENF(US-NR1,FI-Hyy);MF(US-PFa,US-Syv)

#-------------------------------------------------------------------------
##(2).load the data
#-------------------------------------------------------------------------
load(paste0("./data/data_used/","ddf_align_data_from_sos.RDA"))
#put all the years data from each site together:
df_all<-rbind(df_len5_nonnorm$df_dday,df_len5_nonnorm$df_noevent_dday)

#-------------------------------------------------------------------------
##(3).making the plots
#-------------------------------------------------------------------------
#prepare the plot function:
plot_2groups<-function(df,site,comp_var,var_unit,do_legend){
  # df<-df_all
  # site<-"US-UMB"
  # comp_var<-"temp_min_fluxnet2015"
  # var_unit<-"(degreeC)"
  # do_legend<-TRUE
  
  #------------------------------------
  #data tidy and select the relevant variables
  #------------------------------------
  df_sel<-df%>%
    filter(sitename==site)
  ###
  #selected the most relevant vars
  #for the comparison,select the original variable if "do_norm"==FALSE,otherwise "do_norm"==TRUE
  df.sel_temp<-df_sel[,c("sitename","Year","date","dday","gpp_res",comp_var)] #gpp_res:gpp_mod-gpp_obs
  names(df.sel_temp)<-c("sitename","Year","date","dday","bias","comp_var")
  ##do the summary analysis to determine the orders of years that are overestimated:
  summary.bias<-df.sel_temp %>%
    filter(dday>=0)%>%
    group_by(Year)%>%
    dplyr::summarise(spring_bias=round(mean(bias,na.rm=T),3))
  df.all<-left_join(df.sel_temp,summary.bias)
  df.all$spring_bias<-as.factor(df.all$spring_bias)
  #--------------------------------------
  #plotting
  #--------------------------------------
  #y axis range:
  ymin<-min(range(df.sel_temp$comp_var,na.rm = T))
  ymax<-max(range(df.sel_temp$comp_var,na.rm = T))
  #x axis range
  # x_range_event<-range(df.event_sel$doy)
  # x_range_nonevent<-range(df.nonevent_sel$doy)
  #

  ####start to make the plots:
  p_plot<-df.all%>%
    ggplot(aes(x=dday,y=comp_var,col=spring_bias))+
    geom_point(size=1)+
    geom_smooth(se=FALSE,method = "lm")+
    scale_color_viridis_d()+
    ylab(paste0(comp_var," ",var_unit))+
    xlab("gday")+
    annotate("rect",xmin=0,xmax=70,ymin = -Inf,ymax = Inf,alpha=0.2)+
    theme_classic()+
    theme(legend.background = element_blank(),
          legend.key.size = unit(2, 'lines'),
          legend.title = element_text(size=24),
          legend.text = element_text(size=18),
          axis.title = element_text(size=20),
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
  
  return(p_plot)
}

##
# df<-df_all
# site<-"US-UMB"
# comp_var<-"temp_min_fluxnet2015"
# 
# do_legend<-TRUE
#DBF(US-UMB;RU-Fyo);ENF(US-NR1,FI-Hyy);MF(US-PFa,US-Syv)
plot_2groups(df_all,"US-Syv","temp_min_fluxnet2015","(degreeC)",TRUE)
