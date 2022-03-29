#########################################################
#Site selection:Basic infos about the analysis-->the selected sites in Fluxnet2015
#########################################################
#refer Beni's code:https://rpubs.com/stineb/photocold_siteselection
#select sites with a winter-dormant cliamte (temperate and boreal) without dry season
#i.e. K?ppen climate: Cfa, Cfb, Cfc;Dfa, Dfb, Dfc, Dfd
#at the moment, only use forest-->classid:MF, ENF, DBF, and DNF
#---------------------
#selecting the sites fulfilling the above requirements
#---------------------
devtools::load_all("D:/Github/ingestr/")
library(ingestr)
#selecting the sites that in Climate zone "Cf.." and "Df..":
#The sites with no dry season in temperate and boreal region
df_sites_sel<-siteinfo_fluxnet2015 %>%
  filter(koeppen_code %in% c("Cfa", "Cfb", "Cfc", "Dfa", "Dfb", "Dfc", "Dfd") &
         classid %in% c("MF","ENF","DBF","DNF"))
#----------------
#save the data
#----------------
save.path<-"./data-raw/raw_data/sites_info/"
save(df_sites_sel,file=paste0(save.path,"Pre_selected_sites_info.RDA"))
