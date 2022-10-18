##---------------------------------------
#Aim: To add barplots/different colors on the map to indicate different springtime bias 
#in different site years
##---------------------------------------
library(ggplot2)
library(ggmap)
library(rgdal)
library(rgeos)
library(maptools)
library(dplyr)
library(tidyr)
library(tmap)
library(mapdata)
library(rworldmap)
library(rworldxtra)
library(colorRamps)
library(graphics)
library(jpeg)
#I.load the map 
data(coastsCoarse)
#prepration for map
# newmap <- getMap(resolution = "high")[getMap()$ADMIN!='Antarctica',]
newmap <- getMap(resolution = "high")

#(1)load the EC sites used from Fluxnet 2015
library("readxl")
sites.path<-"./data-raw/raw_data/sites_info/"
load(paste0(sites.path,"Pre_selected_sites_info.RDA"))
ori_coord_sites<-df_sites_sel
#at the end, seveals sites do not used for analysis because of the data availabilty 
#-->no complete spring 
#load the final sites name
load.path<-"./data/event_length/"
load(paste0(load.path,"df_events_length.RDA"))
final.sites<-unique(df_events_all$sitename)
match.pos<-match(final.sites,ori_coord_sites$sitename)
final_coord_sites<-ori_coord_sites[match.pos,]

#(2)add the flag and set the factors
ori_coord_sites$flag<-rep("Ori-sites",nrow(ori_coord_sites))
final_coord_sites$flag<-rep("Final_sites",nrow(final_coord_sites))

####set the factors
#for sites from ori sites
ori_coord_sites$koeppen_code<-factor(ori_coord_sites$koeppen_code,levels = c("Cfa","Cfb","Dfb","Dfc"))
ori_coord_sites$PFT<-factor(ori_coord_sites$classid,levels = c("DBF","MF","ENF"))
# #for final sites
final_coord_sites$koeppen_code<-factor(final_coord_sites$koeppen_code,levels = c("Cfa","Cfb","Dfb","Dfc"))
final_coord_sites$PFT<-factor(final_coord_sites$classid,levels = c("DBF","MF","ENF"))

#(3)also adding the information of how many percentage of year with spring bias%
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
#merge the datasets
final_coord_sites<-left_join(final_coord_sites,event_merge)

#(4)load the GPP data(has been characterized as
#"overestimated" or "non-overestimated")
load.path<-"./data/data_used/"
#from new method:
load(paste0(load.path,"ddf_labeled_norm_trs_newmethod_all_overestimation_Fluxnet2015_sites.RDA"))
df_all_sites<-ddf_labeled;rm(ddf_labeled) 
#select the needed variables 
df_GPP<-df_all_sites%>%
  select(sitename,date,doy,greenup,gpp_obs,gpp_mod_FULL,gpp_res)%>%
  filter(greenup=="yes")%>% ##selecting the green-up period data
  group_by(sitename)%>%
  dplyr::summarise(gpp_obs_GP=mean(gpp_obs,na.rm=T),
         gpp_mod_GP=mean(gpp_mod_FULL,na.rm=T),
         gpp_bias=mean(c(gpp_mod_GP - gpp_obs_GP),na.rm=T)
         )   #GP:for green-up period
##merge the data:
final_coord_sites<-left_join(final_coord_sites,df_GPP)

#---------------------------------------------
#II.plotting
#---------------------------------------------
library(RColorBrewer)
library(grDevices)
############
# map theme
############
#can refer:http://www.sthda.com/english/wiki/ggplot2-themes-and-background-colors-the-3-elements
theme_map <- 
  # theme_dark() +    # theme_minimal()
  theme(
    #add by YP:
    # panel.background = element_rect(fill = "gray60",
    #                                 colour = "gray60",
    #                                 size = 0.5, linetype = "solid"),
    # #add by YP:
    # plot.background = element_rect(fill="gray60"),
    #
    plot.title = element_text(hjust = 0, face="bold", size = 18),
    
    # legend.position = "right", # c(0.07, 0.35), #"left"
    # legend.key.size = unit(c(5, 1), "mm"),
    legend.title=element_text(size=12),
    legend.text=element_text(size=10),
    
    # axis.line = element_blank(),
    # axis.text = element_blank(),
    # axis.title = element_blank(),
    
    # panel.grid.major = element_line(colour="black",size = 0.5,linetype = "solid"),
    panel.grid.minor = element_blank(),
    # plot.margin = unit( c(0, 0, 0, 5) , "mm")
  )

# define labels
lat.labels <- seq(30, 90, 30)
lat.short  <- seq(30, 90, 10)
lon.labels <- seq(-180, 180, 60)
lon.short  <- seq(-180, 180, 10)

a <- sapply( lat.labels, function(x) if (x>0) {parse(text = paste0(x, "*degree ~ N"))} else if (x==0) {parse(text = paste0(x, "*degree"))} else {parse(text = paste0(-x, "*degree ~ S"))} )
b <- sapply( lon.labels, function(x) if (x>0) {parse(text = paste0(x, "*degree ~ E"))} else if (x==0) {parse(text = paste0(x, "*degree"))} else {parse(text = paste0(-x, "*degree ~ W"))} )
#update in March, 2022
a<-expression("30°" ~ N, "60°" ~ N, "90°" ~ N)
b<-expression("30°"  ~ W, "120°" ~ W, "60°" ~ W,"0°",
              "60°" ~ E, "120°" ~ E, "180°" ~ E)
#---------------------------------------------
# 1. Create ggplot object
#---------------------------------------------
lonmin=-180
lonmax=180
latmin=30
latmax=90
#group=group-->results in the wrong map background:ask Beni's advices
#something need to be paid attention-->to make sure the plot looks right
#-->should set latmin=-90; latmax=90
#-->and also leave some place for latitude and longtitude-->set the limits in scale_x/y_continous adding or minus some numbers
gg <- ggplot() +
  theme_map+
  # background countries
  # geom_polygon(data=newmap, aes(long, lat, group=group), color=NA, fill='grey75') +
  # Coastline
  geom_path(data=coastsCoarse, aes(long, lat, group=group), color='black',size=1.02) +
  # 
  scale_x_continuous(expand = c(0,0), limits = c(-1+lonmin,lonmax+1), breaks = lon.labels, labels = b) +
  scale_y_continuous(expand = c(0,0), limits = c(-1+latmin,latmax+1), breaks = lat.labels, labels = a) +
  labs( x = "longtitude", y = "latitude")
##update in 2022,Oct-->using Beni' code to create a empty global map
source("./R/functions_from_beni/plot_map_simpl.R")
gg<-plot_map_simpl(-180,180,30,90)
#---------------------------------------------
# 2. add sites information
#---------------------------------------------
library(ggrepel)  #add the site labels
barwidth = 1.5
barheight = 2

#figure format 1:barplot-->
p_final<-gg+
  # geom_point(data=final_coord_sites,aes(x=lon,y=lat,shape=PFT,col=event_perc),size=4,pch=16)+
    geom_rect(data = final_coord_sites,
            aes(xmin = lon - barwidth,
                xmax = lon + barwidth,
                ymin = lat,
                ymax = lat+barheight*event_perc,fill=PFT)) +
  scale_color_gradientn(
    colours = c("red", "white", "blue"),
    values = c(0, 0.5, 1))+
  scale_fill_manual(values = c("DBF"="orange","MF"="cyan","ENF"="magenta"))+
  geom_point(data=final_coord_sites%>% filter(event_perc==0),
             aes(x=lon,y=lat),size=2,pch=16,col="black")+
  geom_text(data = final_coord_sites %>% filter(sitename=="DK-Sor"),
            aes(x = lon+8,
                y = lat + 0.1*event_perc*barheight,
                label = paste0(event_perc*100," %")
                ),
            size = 3,col="orange")+ #only indicate the sites identified with non-overestimation(0%)
  geom_label_repel(data=final_coord_sites,
                   aes(x=lon,y=lat,label = sitename),col="blue",
                   label_size=NA,alpha=0.8,label.padding = .1,
                   max.overlaps = 50,label.size = 0.1,
                   arrow = arrow(ends = "first",length = unit(0.05,"inch")),
                   size = 2.8)
#save the plots
save.path<-"./manuscript/figures/"
ggsave(file=paste0(save.path,"FigureS_sites_distribution_with_barplots.png"),
       p_final,dev="png",width = 12,height=7)

#figure format 2:color for different bias-->
p_final<-gg+
  geom_point(data=final_coord_sites,aes(x=lon,y=lat,col=gpp_bias),size=4,pch=16)+
    scale_color_gradientn("GPP bias",
    colours = c("red", "white", "blue"),
    values = c(0, 0.5, 1))
  # scale_fill_manual(values = c("DBF"="orange","MF"="cyan","ENF"="magenta"))+
  # geom_point(data=final_coord_sites%>% filter(event_perc==0),
  #            aes(x=lon,y=lat),size=2,pch=16,col="black")+
  # geom_text(data = final_coord_sites %>% filter(sitename=="DK-Sor"),
  #           aes(x = lon+8,
  #               y = lat + 0.1*event_perc*barheight,
  #               label = paste0(event_perc*100," %")
  #           ),
  #           size = 3,col="orange")+ #only indicate the sites identified with non-overestimation(0%)
  # geom_label_repel(data=final_coord_sites,
  #                  aes(x=lon,y=lat,label = sitename),col="blue",
  #                  label_size=NA,alpha=0.8,label.padding = .1,
  #                  max.overlaps = 50,label.size = 0.1,
  #                  arrow = arrow(ends = "first",length = unit(0.05,"inch")),
  #                  size = 2.8)
#save the plots
save.path<-"./manuscript/figures/"
ggsave(file=paste0(save.path,"FigureS_sites_distribution_with_colored_bias.png"),
       p_final,dev="png",width = 12,height=7)
