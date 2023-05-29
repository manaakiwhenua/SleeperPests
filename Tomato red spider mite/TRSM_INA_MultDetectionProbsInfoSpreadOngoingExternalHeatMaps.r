##########################################################################################################
###Implements custom function for mapping invasion probability (across simulation runs) in each node in each year
##########################################################################################################

library(ggplot2)
library(magick)
library(dplyr)
library(tidyr)
plot_population <- function(x, y, z, title="", units="", xlim=FALSE, ylim= FALSE, zero_colour=FALSE) {
  gg <- ggplot(data.frame(x=x, y=y, N=z), aes(x=x, y=y, fill=N)) +
    geom_tile(height=6, width=5) + 
    coord_fixed(ratio = 0.8) +                    # sets canvas proportions
    labs(x="", y="", title=title, fill=units) +   # removes xyz titles, adds title
    theme_void() +                                # removes background      
    theme(plot.title=element_text(hjust=0.5)) +   # centres the title
    guides(fill=guide_colourbar(barwidth=1, barheight=15))
  if (xlim && ylim) { 
    gg <- gg + coord_cartesian(xlim=xlim, ylim=ylim) 
  }
  if (zero_colour) {
    gg <- gg + scale_fill_gradientn(colours=c("cornsilk2", "darkseagreen3", "red"), values=c(0,0.05,1))
  } else {
    gg <- gg + scale_fill_gradientn(limits = c(0,1),colours=c("royalblue4", "yellow2", "red"), values=c(0,0.5,1))
  }
  print(gg)
}



input_file   <- "Tomato red spider mite data all NZ.csv"
d <- read.csv(input_file)
ResultsDir = "MultDetProbsInfoSpreadOngoingThreat/"

Nyears = 100

DetectionProbs = c(0.05,0.1,0.15,0.2,0.3,0.5,0.75,0.9)
if(DoClimateChange == TRUE)
	DetectionProbs = c(0.05,0.2,0.5,0.75,0.9,0.95,0.99)
for(detprob in 1:length(DetectionProbs))
{
ModelName = paste0("StandardCurrentClim_DetProb_",DetectionProbs[detprob])
if(DoClimateChange == TRUE)
ModelName = paste0("StandardFutureClim_DetProb_",DetectionProbs[detprob])
ImageDir = paste0(ResultsDir,"HeatMaps/",ModelName,"/")
dir.create(ImageDir,showWarnings = F,recursive = T)
Result = readRDS(paste0(ResultsDir,ModelName,"InvasionProb.rds"))
Nyears = ncol(Result)
for(year in 1:Nyears)
{
ImageFile = paste0(ImageDir, ModelName,"_year_", formatC( x = year, width = 3, flag = "0"),".png")
png(ImageFile)
plot_population(d$x, d$y, Result[,year], paste("Invasion probability year = ", year))
dev.off()
}
imgs = paste0(ImageDir, ModelName,"_year_", formatC( x = 1:Nyears, width = 3, flag = "0"),".png")
img_list <- lapply(imgs, image_read)

## join the images together
img_joined <- image_join(img_list)

## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = 2)

## save to disk
image_write(image = img_animated,
            path = paste0(ImageDir,ModelName,".gif"))

ModelName = paste0("OngoingCurrentClim_DetProb_",DetectionProbs[detprob])
if(DoClimateChange == TRUE)
ModelName = paste0("OngoingFutureClim_DetProb_",DetectionProbs[detprob])
ImageDir = paste0(ResultsDir,"HeatMaps/",ModelName,"/")
dir.create(ImageDir,showWarnings = F,recursive = T)
Result = readRDS(paste0(ResultsDir,ModelName,"InvasionProb.rds"))
for(year in 1:Nyears)
{
ImageFile = paste0(ImageDir, ModelName,"_year_", formatC( x = year, width = 3, flag = "0"),".png")
png(ImageFile)
plot_population(d$x, d$y, Result[,year], paste("Invasion probability year = ", year))
dev.off()
}
imgs = paste0(ImageDir, ModelName,"_year_", formatC( x = 1:Nyears, width = 3, flag = "0"),".png")
img_list <- lapply(imgs, image_read)

## join the images together
img_joined <- image_join(img_list)

## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = 2)

## save to disk
image_write(image = img_animated,
            path = paste0(ImageDir,ModelName,".gif"))

}
