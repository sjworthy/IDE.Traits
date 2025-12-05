# Code for making Figure 1: the map

# Youtube video used to help make the map
# https://www.youtube.com/watch?v=9Ex-f66qe2w

# load packages

library(pRecipe)
library(giscoR)
library(raster)
library(terra)
library(rayshader)
library(classInt)
library(tidyverse)
library(stringr)

# read in the lat/longs for the 83 IDE sites included
sites.63 = read.csv("./Formatted.Data/Revisions/IDE.sites.63.csv")

# get the whole world
world_sf = gisco_get_countries()

#download_data(dataset = "mswep",
              #path = getwd(),
              #timestep = "yearly")

mswep_data = terra::rast("./mswep_tp_mm_global_197902_202301_025_yearly.nc")
terra::ext(mswep_data) <- c(-180, 180, -90, 90)

names(mswep_data) = 1979:2023

mswep_df = mswep_data %>%
  as.data.frame(xy = TRUE)

# remove 2023, not complete data
mswep_df = mswep_df[,c(1:46)]
mswep_df_2 = mswep_df[,c(3:46)]
mswep_df$MAP = rowMeans(mswep_df_2)

terra::plot(mswep_data[[1]])
plot(sf::st_geometry(world_sf), add=TRUE)

theme_for_the_win = function(){
  theme_minimal() +
    theme(
      axis.line = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 11, color = "gray10"),
      legend.text = element_text(size = 10, color = "gray10"),
      legend.text.align = 0.5,
      panel.grid.major = element_line(color = NA),
      panel.grid.minor = element_line(color = NA),
      plot.background = element_rect(fill = NA, color = NA),
      legend.background = element_rect(fill = "white", color="black"),
      panel.border = element_rect(fill=NA, color=NA)
      #plot.margin = unit(c(t=0,r=0,b=0,l=0), "lines")
    )
}

breaks = classInt::classIntervals(mswep_df$MAP, n = 5, style = "pretty")$brks

colors = hcl.colors(n = length(breaks),rev = TRUE)

map1 = ggplot(data = mswep_df) +
  geom_raster(aes(x = x, y = y, fill = MAP)) +
  geom_sf(data = world_sf, fill = "transparent", color = "gray10",size = .5) +
  scale_fill_gradientn(name = "Mean Annual Precipitation (mm)",
                       colors = colors,
                       breaks = breaks,
                       limits = c(min(mswep_df$MAP),max(mswep_df$MAP)),
                       labels = round(breaks, 0),
                       trans = "log10")+
  guides(fill = guide_colorbar(direction = "vertical", barheight = unit(50, "mm"),
                               barwidth = unit(5, "mm"), title.position = "top", 
                               label.position = "right",title.hjust = .5, label.hjust = .5,
                               ncol=1, byrow=FALSE)) +
  theme_for_the_win()

map1

map1 = ggplot(data = mswep_df) +
  geom_raster(aes(x = x, y = y, fill = MAP)) +
  geom_sf(data = world_sf, fill = "transparent", color = "gray10",size = .5) +
  scale_fill_gradient2(name = str_wrap("Mean Annual Precipitation (mm) 1979-2022", 20),
                       low = "red3", mid = "#E2E2E2",high = "blue3",
                       trans = "log10",midpoint = log10(500))+
  geom_point(data = sites.63,
             aes(x = longitud, y = latitud), color = "black", size = 1.5,alpha = 0.7,
             position=position_jitter(h=0.50,w=0.50))+
  guides(fill = guide_colorbar(title.position = "top", 
                               label.position = "right",title.hjust = .5, label.hjust = .5)) +
  theme_for_the_win_legend()
map1  


theme_for_the_win_legend = function(){
  theme_minimal() +
    theme(
      axis.line = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "none",
      legend.title = element_text(size = 11, color = "gray10"),
      legend.text = element_text(size = 10, color = "gray10"),
      legend.text.align = 0.5,
      panel.grid.major = element_line(color = NA),
      panel.grid.minor = element_line(color = NA),
      plot.background = element_rect(fill = NA, color = NA),
      legend.background = element_rect(fill = "white", color="black"),
      panel.border = element_rect(fill=NA, color=NA)
      #plot.margin = unit(c(t=0,r=0,b=0,l=0), "lines")
    )
}