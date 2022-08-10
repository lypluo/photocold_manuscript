f_hardening <- function(temp, par){
  #changing the function according to the logistic function in weikipeidia:
  #https://en.wikipedia.org/wiki/Logistic_function
  xx_ori <- temp # * ppfd
  xx <- (-1)*par["b"] * (xx_ori - par["a"])
  yy <- 1 / (1 + exp(xx))
  plot(xx_ori,yy,xlab="Ta",ylab="par_stress")
  return(yy)
}

##
temp<-seq(-40,70,0.5)

#
par1<-c("a"=5,"b"=0.2)
f_hardening(temp,par1)

par1<-c("a"=5,"b"=0.1)
f_hardening(temp,par1)

par1<-c("a"=-5,"b"=0.1)
f_hardening(temp,par1)

#officially plotting:
library(ggplot2)
#-------------------------
#for "hardening"
#-------------------------
x=temp
par1<-c("a"=20,"b"=0.15)
y=f_hardening(temp,par = par1)
#
df1<-data.frame(x=x,y=y)
df_sel1<-df1[1:90,]
ggplot()+
  geom_line(data=df1,aes(x=x,y=y),col=adjustcolor("tomato",0.8),lty=2,size=1.05)+
  geom_line(data=df_sel1,aes(x=x,y=y),size=1.1,col="tomato")+
  theme_bw()
  # theme_void()

#-------------------------
#for "dehardening"
#-------------------------
x=temp
par2<-c("a"=12.5,"b"=0.3)
y=f_hardening(temp,par = par2)
#
df2<-data.frame(x=x,y=y)
df_sel2<-df2[91:nrow(df2),]
ggplot()+
  geom_line(data=df2,aes(x=x,y=y),col=adjustcolor("forestgreen",0.8),lty=2,size=1.05)+
  geom_line(data=df_sel2,aes(x=x,y=y),size=1.1,col="forestgreen")+
  theme_bw()
  # theme_void()

#--------------------------
#merge "hardenning" and "dehardening" period
#-------------------------

#
# p_plot<-ggplot()+
#   geom_line(data=df1,aes(x=x,y=y),col=adjustcolor("tomato",0.8),lty=2,size=1.05)+
#   geom_line(data=df_sel1,aes(x=x,y=y),size=1.1,col="tomato")+
#   geom_line(data=df2,aes(x=x,y=y),col=adjustcolor("forestgreen",0.8),lty=2,size=1.05)+
#   geom_line(data=df_sel2,aes(x=x,y=y),size=1.1,col="forestgreen")+
#   xlab("")+
#   # xlab(expression("T"[min]*"( °C )"))+
#   ylab(expression(italic("f")[stress]))+
#   theme_bw()+
#   theme(axis.title = element_text(size = 18),
#         axis.ticks.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.text = element_text(size=16))
#save the plot:
# path<-"./manuscript/figures/"
# ggsave(p_plot,filename = paste0(path,"f_stress.png"),height = 4,width = 5.5)


################################################
#re-make the plots:
#-------------------------
#for "hardening"
#-------------------------
temp<-seq(-50,70,0.5)
x=temp
par1<-c("a"=40,"b"=0.15)
y=f_hardening(temp,par = par1)
#
df1<-data.frame(x=x,y=y)
df_sel1<-df1[1:90,]
p_plot1<-ggplot()+
  geom_line(data=df1,aes(x=x,y=y),col=adjustcolor("tomato",0.8),lty=1,size=1.05)+
  # geom_line(data=df_sel1,aes(x=x,y=y),size=1.1,col="tomato")+
  xlab(expression("T"[min]*"( °C )"))+
  ylab(expression(italic("f")[1]))+
  # ggtitle("Hardening")+
  ggtitle("Phase 1")+
  theme_bw()+
  theme(axis.title = element_text(size = 18),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text = element_text(size=16),
        plot.title = element_text(hjust=0.5))
  
# theme_void()

#-------------------------
#for "dehardening"
#-------------------------
temp<-seq(-20,50,0.5)
par2<-c("a"=12.5,"b"=0.3)
y=f_hardening(temp,par = par2)
#
df2<-data.frame(x=x,y=y)
df_sel2<-df2[91:nrow(df2),]
p_plot2<-ggplot()+
  geom_line(data=df2,aes(x=x,y=y),col=adjustcolor("forestgreen",0.8),lty=1,size=1.05)+
  # geom_line(data=df_sel2,aes(x=x,y=y),size=1.1,col="forestgreen")+
  xlab(expression("GDD"*"( °C )"))+
  ylab(expression(italic("f")[2]))+
  # ggtitle("Dehardening")+
  ggtitle("Phase 2")+
  theme_bw()+
  theme(axis.title = element_text(size = 18),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text = element_text(size=16),
        plot.title = element_text(hjust=0.5)
        )
# theme_void()

##merge two plots:
library(cowplot)
library(grid)

path<-"./manuscript/figures/"
p_plot<-plot_grid(p_plot1,p_plot2,ncol=2)
ggsave(p_plot,filename = paste0(path,"f_stress.png"),height = 3,width = 6)
