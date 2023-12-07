# map plot function adopt from "rbeni"(https://github.com/stineb/rbeni/blob/master/R/plot_map_simpl.R) 
#-->and revised by Yunpeng:updated in 2022-Oct
#' Empty global map
#'
#' Returns a ggplot object for an empty global map
#'
#' @param lonmin Left edge (longitude, in degrees), defaults to -180.
#' @param lonmax Right edge (longitude, in degrees), defaults to 180.
#' @param latmin Lower edge (latitude, in degrees), defaults to -90.
#' @param latmax Upper edge (latitude, in degrees), defaults to 90.
#' @return A ggplot object for a global map plot.
#' @export
#'
plot_map_simpl <- function(lonmin = -180, lonmax = 180, latmin =30, latmax = 85){
  
  library(rworldmap)
  library(ggplot2)
  library(rnaturalearth)
  library(sf)
  #
  lonmin= -30
  lonmax = 60
  latmin=30
  latmax=80 ##latmax should be higher 75
  
  # download global coastline data from naturalearth
  countries <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf")
  
  ##---------------------------------------------
  ## Projection
  ##---------------------------------------------
  # set coordinate systems
  # robinson <- CRS("+proj=robin +over")
  
  # create a bounding box for the robinson projection
  bb <- sf::st_union(sf::st_make_grid(
    st_bbox(c(xmin = lonmin,
              xmax = lonmax,
              ymax = latmax,
              ymin = latmin),
            crs = st_crs(4326)),
    n = c(12,6)))
  # bb_robinson <- st_transform(bb, as.character(robinson))
  
  # clip countries to bounding box
  # and transform
  # countries <- countries %>%
  #   st_buffer(0) %>%
  #   st_intersection(st_union(bb)) %>%
  #   st_transform(robinson)
  
  
  ## map theme
  ##---------------------------------------------
  theme_map <- theme_grey() +    # theme_minimal()
    theme(
      
      plot.title = element_text(hjust = 0, face="bold", size = 24),
      
      legend.position = "right", # c(0.07, 0.35), #"left"
      # legend.key.size = unit(c(5, 1), "mm"),
      legend.title=element_text(size=24),
      legend.text=element_text(size=20),
      
      # axis.line = element_blank(),
      axis.text = element_text(size = 20),
      # axis.title = element_blank(),
      
      # panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
      # plot.margin = unit( c(0, 0, 0, 5) , "mm")
    )
  
  ##---------------------------------------------
  ## Create ggplot object
  ##---------------------------------------------
  #update in October, 2022 by Yunpeng
  # create the breaks- and label vectors
  #update in March, 2022
  a<-expression("30°" ~ N, "60°" ~ N, "90°" ~ N)
  b<-expression("30°"  ~ W, "120°" ~ W, "60°" ~ W,"0°",
                "60°" ~ E, "120°" ~ E, "180°" ~ E)
  
  ewbrks <- seq(-120, 120, 60)
  nsbrks <- seq(45,75,15)
  ewlbls <- unlist(lapply(ewbrks, function(x) ifelse(x < 0, paste(-x, "°W"), ifelse(x > 0, paste(x, "°E"),x))))
  nslbls <- unlist(lapply(nsbrks, function(x) ifelse(x < 0, paste(-x, "°S"), ifelse(x > 0, paste(x, "°N"),x))))
  #find out the reason why I cannot add the maually label(long./lat.) into x and y axis
  #-->because the breaks cannot surpass the data range boundary-->for instance:
  #for range(-180,180), the breaks cannot be -180, and 180...
  
  #
  # sf_use_s2(FALSE)
  gg <- ggplot() +
    # background countries
    geom_sf(data = countries, color="black", fill='grey75', size = 0.1) +
    # geom_path(data=coastsCoarse, aes(long, lat, group=group), color='black',size=1.02)
    geom_sf(data = countries,
            colour = 'black',
            fill = 'grey75',  # NA for empty
            linetype = 'solid',
            size = 0.1) +
    # bounding box
    # geom_sf(data = bb,
    #         colour = 'black',
    #         linetype = 'solid',
    #         fill = NA,
    #         size = 0.1) +
    coord_sf(xlim = c(lonmin,lonmax),ylim=c(latmin,latmax),expand = FALSE)+
    # coord_sf(
    #   ylim = c(-60, 90)
    # ) +
    #manually add it-->!!and also the break should not outside the boundary, 
    #for instance:breaks cannot include -180 and 180 of longtitude
    scale_x_continuous(expand = c(0, 0),limits = c(lonmin,lonmax),
                       breaks =ewbrks,
                       labels = ewlbls)+
    scale_y_continuous(expand = c(0, 0),limits = c(latmin,latmax),
                       breaks =nsbrks,
                       labels = nslbls)+
    # scale_x_continuous(expand = c(0, 0),limits = c(-180,180),
    #                    breaks =c(-120,-60,0,60,120),
    #                    labels = c("120° W","60° W","0","60° E","120° E"))+  
    # scale_y_continuous(expand = c(0, 0),limits = c(30,90),
    #                    breaks =c(45,60,75),
    #                    labels = c("45° N","60° N","75° N"))+  
    # scale_x_continuous(expand = c(0, 0),breaks = ewbrks,labels = ewlbls) +
    # scale_y_continuous(expand = c(0, 0),breaks = nsbrks,labels = nslbls) +
    labs( x = "longtitude", y = "latitude")+
    theme_bw()+
    theme(axis.text = element_text(size=16),
          axis.title = element_text(size=18))
  # theme(axis.ticks.y.right = element_line(),
  #       axis.ticks.x.top = element_line(),
  #       panel.grid = element_blank())
  
  return(gg)
}
