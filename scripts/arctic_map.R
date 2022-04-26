# Arctic map from stack over flow

library(rgdal)                                                                                                      
library(raster)
library(ggplot2)

data("wrld_simpl", package = "maptools") 
wm <- crop(wrld_simpl, extent(-180, 180, 45, 90))            
wm <- spTransform(wm, CRSobj = CRS(proj))

lat_long <- c("lat", "long")
qhi_point <- c(69.6, -138.9)

qhi_coords

# Defines the x axes required
x_lines <- seq(-120,180, by = 60)

ggplot() +
  geom_polygon(data = wm, aes(x = long, y = lat, group = group), fill = "grey", colour = "black", alpha = 0.8) +
  geom_point(data = qhi_coords, aes(x = long, y = lat), color = "red", alpha = 0.9, size = 5) +
  # Convert to polar coordinates
  coord_map("ortho", orientation = c(90, 0, 0)) +
  scale_y_continuous(breaks = seq(45, 90, by = 5), labels = NULL) +
  
  # Removes Axes and labels
  scale_x_continuous(breaks = NULL) +
  xlab("") + 
  ylab("") +
  # Adds labels
  #geom_text(aes(x = 180, y = seq(55, 85, by = 10), hjust = -0.2, label = paste0(seq(55, 85, by = 10), "°N"))) +
  #geom_text(aes(x = x_lines, y = 39, label = c("120°W", "60°W", "0°", "60°E", "120°E", "180°W"))) +
  
  # Adds axes
  #geom_hline(aes(yintercept = 45), size = 1)  +
  #geom_segment(aes(y = 45, yend = 90, x = x_lines, xend = x_lines), linetype = "dashed") +
  
  # Change theme to remove axes and ticks
  theme(panel.background = element_blank(),
       # panel.grid.major = element_line(size = 0.25, linetype = 'dashed',
        #                                colour = "black"),
        axis.ticks=element_blank()) 
  #labs(caption = "Designed by Mikey Harper")

## Yukon MAP

library(maps)

# Check all available geospatial objects:
# help(package='maps')

# Map of the world:
map('world',col="grey", fill=TRUE, bg="white", lwd=0.05, mar=rep(0,4), ylim=c(69.6549,69.5),xlim = c(-139.3698,-138.8186))
map("canada")
help(package='maps')
