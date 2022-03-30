#############################################################################
#Aim:tidy the stats from the "is_event":the dates when GPP is overestimation.
#-------------------------------------------------------------------------
#(1)load the data
#-------------------------------------------------------------------------
load.path<-"./data/data_used/"
#from new method:
load(paste0(load.path,"ddf_labeled_norm_trs_newmethod_all_overestimation_Fluxnet2015_sites.RDA"))
df_norm_trs_newM_allover<-ddf_labeled;rm(ddf_labeled)

#(2)select the growing season period and stats such as "is_event"
df_norm_newM_allover_sel<-df_norm_trs_newM_allover[,c("sitename","date","doy","greenup","GS","over_estim","is_event","is_event_less10")]

#
library(reshape2)
df_all<-df_norm_newM_allover_sel
df.event<-melt(df_all,id.vars = c("sitename","date","doy"))
names(df.event)<-c("sitename","date","doy","category","is_event")
#------------------------------------------------------------------------
#(2)plotting
#---------------------------------
##histogram for is_event
#---------------------------------
library(ggplot2)
library(plyr)
#calculate the mean of doy
df.plot<-df.event[df.event$is_event=="yes",]
df.plot$category<-factor(df.plot$category,levels = c("GS","over_estim","greenup","is_event","is_event_less10"))
mu <- ddply(df.plot, "category", summarise, grp.mean=mean(doy))
p_hist<-ggplot(df.plot, aes(x=doy, col=category))+
  geom_histogram(binwidth = 1,fill="white",
                 position="identity")+
  # geom_histogram(aes(y=..density..), alpha=0.5,binwidth = 5,
  #                position="identity")+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=category),
             linetype="dashed",size=1.2)+
  scale_color_manual(values=c("GS"="green3","over_estim"="grey","greenup"="skyblue","is_event"="tomato","is_event_less10"="orange"))+
  theme_minimal()

p_density<-ggplot(df.plot, aes(x=doy, col=category,fill=category))+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=category),
             linetype="dashed",size=1.2)+
  # geom_histogram(aes(y=..density..), alpha=0.5, 
  #                position="identity")+
  geom_density(alpha=.4,size=1.2)+
  scale_color_manual(values=c("GS"="green3","over_estim"="grey","greenup"="skyblue","is_event"="tomato","is_event_less10"="orange"))+
  scale_fill_manual(values=c("GS"="green3","over_estim"="grey","greenup"="skyblue","is_event"="tomato","is_event_less10"="orange"))+
  theme_minimal()
#
#library(cowplot)
library(ggpubr)
library(gridExtra)
p_merge<-ggarrange(p_hist,p_density,ncol=2,align = "hv",common.legend=TRUE)
print(p_merge)
# p_merge<-plot_grid(p_hist,p_density,
#                    labels = "auto",ncol=2,label_size = 12,align = "hv")
#save the plot
# save.path<-"C:/Users/yluo/Documents/GitHub/photocold/plot/"
# pdf(paste0(save.path,"is_event(doy).pdf"),width = 8,height=5)
# print(p_merge)
# dev.off()
# plot(p_merge)
# ggsave(file=paste0(save.path,"distribution_of_GS_overestim_greenup_is_event(doy).png"),p_merge,dev="png",width = 10,height=5)

