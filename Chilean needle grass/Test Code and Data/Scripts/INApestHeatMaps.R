
library(magick)
INAheatmapShapeFile = function(ShapeObject, InvasionProbs, ImageDir)
{
dir.create(ImageDir)
for(year in 1:ncol(InvasionProbs))
{
ShapeObject$InvasionProb = InvasionProbs[,year]
Title = paste0(ModelName,"\nYear = ", year)
gg =  ggplot()+
geom_sf(data=ShapeObject, aes(fill=InvasionProb))+
scale_fill_gradientn(limits = c(0,1),colours=c("lightblue2","green4" ,"yellow2", "red"), values=c(0, 0.01,0.05,0.1,0.3,0.5,0.7,0.9,1))+
geom_tile(height=6, width=5) + 
labs(title=Title, fill="Invasion \nprobability")+
theme(plot.title=element_text(hjust=0.5)) +   # centres the title
guides(fill=guide_colourbar(barwidth=1, barheight=15))
ImageFile = paste0(ImageDir, ModelName,"_year_", formatC( x = year, width = 3, flag = "0"),".jpg")
ggsave(ImageFile,dpi = 90)
}

imgs = paste0(ImageDir, ModelName,"_year_", formatC( x = 1:ncol(InvasionProbs), width = 3, flag = "0"),".jpg")
img_list <- lapply(imgs, image_read)

## join the images together
img_joined <- image_join(img_list)

## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = 2)

## save to disk
image_write(image = img_animated,
            path = paste0(ImageDir,ModelName,".gif"))

}

INAheatmapXYpoints = function(XY, InvasionProbs, ImageDir)
{
dir.create(ImageDir)
for(year in 1:ncol(InvasionProbs))
{
XY$Z = InvasionProbs[,year]
Title = paste0(ModelName,"\nYear = ", year)
gg = ggplot(XY,aes(X,Y, color = Z))+
geom_point(size = 3) +
scale_colour_gradientn(limits = c(0,1),colours=c("lightblue2","green4" ,"yellow2", "red"), values=c(0, 0.01,0.05,0.1,0.3,0.5,0.7,0.9,1))+
coord_fixed(ratio = 1)+
labs(title = Title,color="Invasion \nprobability")+
theme(plot.title=element_text(hjust=0.5)) + 
#labs(fill="Invasion \nprobability")+
#    theme_void() +                                # removes background      
    theme(plot.title=element_text(hjust=0.5)) +   # centres the title
    guides(fill=guide_colourbar(barwidth=1, barheight=15))

ImageFile = paste0(ImageDir, ModelName,"_year_", formatC( x = year, width = 3, flag = "0"),".jpg")
ggsave(ImageFile,dpi = 90)
}
imgs = paste0(ImageDir, ModelName,"_year_", formatC( x = 1:ncol(InvasionProbs), width = 3, flag = "0"),".jpg")
img_list <- lapply(imgs, image_read)

## join the images together
img_joined <- image_join(img_list)

## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = 2)

## save to disk
image_write(image = img_animated,
            path = paste0(ImageDir,ModelName,".gif"))
}