library(tidyverse)
library(rworldmap)
library(sf)
# Data 
library(spData)      
library(spDataLarge)

# Get map data
# worldMap <- map_data("world")

# Select only some countries and add values
europe <- 
  data.frame("country"=c("Austria",
                         "Belgium", 
                         "Germany",
                         "Spain", 
                         "Finland", 
                         "France", 
                         "Greece", 
                         "Ireland", 
                         "Italy", 
                         "Netherlands", 
                         "Portugal",
                         "Bulgaria","Croatia","Cyprus", "Czech Republic","Denmark","Estonia", "Hungary",
                         "Latvia", "Lithuania","Luxembourg","Malta", "Poland", "Romania","Slovakia",
                         "Slovenia","Sweden","UK", "Switzerland",
                         "Ukraine", "Turkey", "Macedonia", "Norway", "Slovakia", "Serbia", "Montenegro",
                         "Moldova", "Kosovo", "Georgia", "Bosnia and Herzegovina", "Belarus", 
                         "Armenia", "Albania", "Russia"),
             "Growth"=c(1.0, 0.5, 0.7, 5.2, 5.9, 2.1, 
                        1.4, 0.7, 5.9, 1.5, 2.2, rep(NA, 33)))

# Merge data and keep only Europe map data

data("world")

worldMap <- world

#match the countries:
worldMap<-worldMap %>%
  mutate(country=name_long,name_long=NULL)  ## change the "name_lone" to "country"
##merge the worldmap and europe datasets:
worldMap<-left_join(worldMap,europe)

#set the geometirc centre
centres <-
  worldMap %>%
  filter()%>%
  st_centroid()

worldMap <- worldMap
  # filter(!is.na(worldMap$Growth)) #remove the site do not have the growth data

# Plot it 

centroids <- 
  centres$geom %>% 
  purrr::map(.,.f = function(x){data.frame(long = x[1],lat = x[2])}) %>% 
  bind_rows %>% data.frame(country = centres$country)%>% 
  left_join(europe)


barwidth = 1
barheight = 0.75

ggplot()+ 
  geom_sf(data = worldMap, color = "black",fill = "lightgrey",
          colour = "white", size = 0.1)+
  coord_sf(xlim = c(-13, 35),  ylim = c(32, 71)) + 
  geom_rect(data = centroids,
            aes(xmin = long - barwidth,
                xmax = long + barwidth,
                ymin = lat,
                ymax = lat + Growth*barheight),fill="steelblue2") + 
  geom_text(data = centroids %>% filter(!is.na(Growth)),
            aes(x = long,
                y = lat + 0.5*Growth*0.75,
                label = paste0(Growth," %")),
            size = 2,col="blue")
#save the plot
  ggsave(file = "test.pdf",
         width = 10,
         height = 10)
