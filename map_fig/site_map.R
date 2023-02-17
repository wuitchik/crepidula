# Followed the tutorial from https://www.zoology.ubc.ca/~kgilbert/mysite/Miscellaneous_files/R_MakingMaps.pdf

# libraries
library(sp)
library(maptools)
gpclibPermit()
library(maps)
library(mapdata)
library(mapproj)

sites = read.delim("sites.txt", sep = "\t", header = T)


map(database = "world",
    projection = "gilbert",
    xlim = c(-170, -20),
    ylim = c(20, 80),
    col = "gray80",
    fill = TRUE, 
    orientation = c(90,0,300))



plot = map("worldHires",
    xlim = c(-77, -59),
    ylim = c(35,50),
    col = "gray80", 
    fill = TRUE)

colours = c("Robbinston" = "#264D59", "Kettle Cove" = "#43978D","Beverly" = "#F9E07F", "Newport" = "#F9AD6A","Cape May" = "#D46C4E")


points(sites$Long, sites$Lat, pch = 19, col = colours)

ggsave("plot.pdf", plot)


