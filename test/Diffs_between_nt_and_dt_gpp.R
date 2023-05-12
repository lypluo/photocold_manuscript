######################################################################
#Aim: to test if the GPP_nt(nighttime partioned) and GPP_dt(daytime partioned)
# differs largely after the filtering their extremes values (removing the upper 
#and lower 2.5% quantiles of the difference between GPP values estimated based on these two methods)
######################################################################
library(ggpubr)
library(ggplot2)
library(ggpmisc)

#--load the merged fluxes data
load.path<-"./data-raw/raw_data/Merged_data/"
load(paste0(load.path,"Merged_Flux_and_VIs.RDA"))

#-------------------
#compare between GPP_nt and GPP_dt after filtering:
#-------------------
#linear regression between gpp_dt and gpp_nt
lm_stat<-lm(df_merge$gpp_dt~df_merge$gpp_nt,data=df_merge)
stat.sum<-summary(lm_stat)
p_comp<-ggplot(data=df_merge,aes(x=gpp_nt,y=gpp_dt,col=sitename))+
  geom_point()+
  xlim(-5,20)+
  geom_abline(slope = 1,intercept = 0,lty=2,col="blue")+
  geom_abline(slope = coef(lm_stat)[2],intercept = coef(lm_stat)[1],size=1.1)+
  labs(y = expression( paste("GPP"[DT]* " (g C m"^-2, " d"^-1, ")" ) ),
       x = expression( paste("GPP"[NT]* " (g C m"^-2, " d"^-1, ")" ) ))+
#adding the R2 and linear regression equations:
  annotate(geom="text",x=0,y=19,label = paste0("Slope = ", round(coef(lm_stat)[2],2)))+
  annotate(geom = "text",x=0,y=17,label = paste0("italic(R) ^ 2 == ",round(stat.sum$r.squared,2)),
           parse=TRUE)+
  theme(
    axis.title = element_text(size=22),
    axis.text = element_text(size = 20),
  )

#
save.path<-"./test/comp_gpp_nt_and_dt/"
ggsave(paste0(save.path,"gpp_comparison.png"),p_comp)
