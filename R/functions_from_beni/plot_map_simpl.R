# map plot function adopt from "rbeni"(https://github.com/stineb/rbeni/blob/master/R/plot_map_simpl.R) 
#-->but really does not work!!!-->because the unit of degree is square rather than circle
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
plot_map_simpl <- function(lonmin = -180, lonmax = 180, latmin = -60, latmax = 85){
  
  library(rworldmap)
  library(ggplot2)
  library(rnaturalearth)
  library(sf)
  
  # download global coastline data from naturalearth
  countries <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf")
  
  ##---------------------------------------------
  ## Projection
  ##---------------------------------------------
  # set coordinate systems
  # robinson <- CRS("+proj=robin +over")
  
  # create a bounding box for the robinson projection
  bb <- sf::st_union(sf::st_make_grid(
    st_bbox(c(xmin = -180,
              xmax = 180,
              ymax = 90,
              ymin = 30),
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
      
      plot.title = element_text(hjust = 0, face="bold", size = 18),
      
      legend.position = "right", # c(0.07, 0.35), #"left"
      # legend.key.size = unit(c(5, 1), "mm"),
      legend.title=element_text(size=12),
      legend.text=element_text(size=10),
      
      # axis.line = element_blank(),
      # axis.text = element_blank(),
      # axis.title = element_blank(),
      
      # panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
      # plot.margin = unit( c(0, 0, 0, 5) , "mm")
    )
  
  ##---------------------------------------------
  ## Create ggplot object
  ##---------------------------------------------
  #update in October, 2022 by Yunpeng
  lat.labels <- seq(30, 90, 30)
  lat.short  <- seq(30, 90, 10)
  lon.labels <- seq(-180, 180, 60)
  lon.short  <- seq(-180, 180, 10)

  # a<-expression("30°" ~ N,"40°" ~ N,"50°" ~ N,
  #               "60°" ~ N,"70°" ~ N,"80°" ~ N,"90°" ~ N)
  # a<-expression("30°" ~ N,"60°" ~ N,"90°" ~ N)
  # b<-expression("180°"  ~ W, "120°" ~ W, "60°" ~ W,"0°",
  #               "60°" ~ E, "120°" ~ E, "180°" ~ E)
  
  # lat.labels <- seq(30, 90, 30)
  # lon.labels <- seq(-180, 180, 60)
  # a<-c("30° N","60° N" ,"90° N")
  # b<-c("180° W","120° W", "60° W","0°","60° E","120° E","180° E")
  # create the breaks- and label vectors
  ewbrks <- seq(-180, 180, 60)
  nsbrks <- seq(30,90,30)
  ewlbls <- unlist(lapply(ewbrks, function(x) ifelse(x < 0, paste(-x, "°W"), ifelse(x > 0, paste(x, "°E"),x))))
  nslbls <- unlist(lapply(nsbrks, function(x) ifelse(x < 0, paste(-x, "°S"), ifelse(x > 0, paste(x, "°N"),x))))
  
  #
  # sf_use_s2(FALSE)
  gg <- ggplot() +
    # background countries
    # geom_sf(data = countries, color="black", fill='grey75', size = 0.1) +
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
    # scale_x_continuous(expand = c(0, 0)) +  #all those ways to add the labels dose not work
    # scale_y_continuous(expand = c(0, 0)) +
    # scale_x_continuous(expand = c(0, 0),breaks = ewbrks,labels = ewlbls) +
    # scale_y_continuous(expand = c(0, 0),breaks = nsbrks,labels = nslbls) +
    # scale_x_continuous(expand = c(0, 0), limits = c(-1+lonmin,lonmax+1),
    #                    breaks = lon.labels, labels = b) +
    # scale_y_continuous(expand = c(0, 0), limits = c(-1+latmin, latmax+1),
    #                    breaks = lat.labels, labels = a) +
    labs( x = "longtitude", y = "latitude")+
    theme_bw()
  # theme(axis.ticks.y.right = element_line(),
  #       axis.ticks.x.top = element_line(),
  #       panel.grid = element_blank())
  
  return(gg)
}
