
bats <- read.csv("input/mmpop_bats.csv", header=T)


library(raster)
library(rworldmap)
library(tidyverse)
library(ggspatial)

boundary <- extent(-180, 180, 0, 90)
boundary

map_outline <- getMap(resolution="high")
map_outline <- crop(map_outline, y=boundary) %>% fortify()

# bat colors
batColors <- setNames(c("#F0A0FF", "#0075DC", "#993F00", "#4C005C", "#191919", 
                        "#005C31", "#2BCE48", "#FFCC99", "#94FFB5", 
                        "#8F7C00", "#9DCC00", "#C20088", "#003380", "#FFA405", 
                        "#FFA8BB", "#426600", "#FF0010", "#5EF1F2"), levels(bats$species))


map <- ggplot()+  
  # Plot map land (fill), and outline (colour & size)
  geom_polygon(data=map_outline, aes(x=long, y=lat, group=group), fill="gray70")+#, colour="gray90", size=0.5)+
  # Remove extra space between map & axes
  coord_quickmap(expand=F)+
  # Axes labels
  xlab("Longitude")+
  ylab("Latitude")+
  # Adjust background theme/colour
  theme(
    axis.text = element_text(colour="black", size=8),
    # axis.title = element_text(colour="black", size=14),
    panel.background = element_rect(fill="white"),
    panel.border = element_rect(fill=NA, colour="black", size=0.5),
    legend.text=element_text(size=8),
    legend.title=element_blank(),
    legend.key.size=unit(0.7, "cm"),
    legend.position="none",
    panel.grid.minor=element_blank(),#line(colour="grey90", size=0.5, linetype="dashed"),
    panel.grid.major=element_blank()#ine(colour="grey90", size=0.5, linetype="dashed")
  )+
  scale_x_continuous(breaks=seq(-180,180,50))+
  geom_point(data=bats, aes(x=lon, y=lat, fill=species), size=2.5, shape=21,alpha=0.8)+
  scale_fill_manual(values=batColors)
# annotation_north_arrow(location="bl", which_north="true", pad_x=unit(0.25, "in"), pad_y=unit(0.25, "in"), style = north_arrow_fancy_orienteering)+
# annotation_scale(location="bl", width_hint=0.15)


map

ggsave("bat_map1.png", width=8, height=4, dpi=600)