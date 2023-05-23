
#Packages----
.libPaths()

if ( 1==2)
{
  install.packages("basemaps")
  #install.packages("ggmap")
  install.packages("colorblindcheck")
  install.packages("spacesXYZ")
}
library(INA)
library(tidyverse)
library(dplyr)
library(colorblindcheck)

memory.limit(size=600000) 

getwd()



source(r"[T:\Hamilton\Projects P-T\Spatlab\Development\SL1714_RLibrary_smap-dsm\R\raster-extents.R]")
source(r"[R:/Projects/SL2116_MPI_SleeperPests/Analysis/Chilean needle grass/Updated Analyses 2022/CurrentFutureClimate/Scripts/INASource_String_FindByCode.R]")

#Functions----
fun.plot.yaxis <- function()
{
  lst.scales <- c(1,2,5)
  scale.four <- (par("usr")[2] - par("usr")[1]) / 4
  scale.dist <- log10( scale.four)
  scale.round <- lst.scales[which.min(abs(lst.scales - scale.dist))]
  scale.step <-  as.numeric( scale.round * (10 **  ceiling(log10( scale.four)-1) ))
  xlab.step <- scale.step
  xlab.min <- floor( par("usr")[1]/xlab.step) * xlab.step
  xlab.max <- floor( par("usr")[2]/xlab.step) * xlab.step
  xlab.atx <- seq( from =  xlab.min, to = xlab.max, by = xlab.step)
  for ( i.loop in 1:length( xlab.atx))
  {
    axis( at = xlab.atx[i.loop], labels =  fun.string.exp(xlab.atx[i.loop]), side = 1, cex = 0.5)
  }
  
  ylab.step = scale.step
  ylab.min <- floor( par("usr")[3]/ylab.step) * ylab.step
  ylab.max <- floor( par("usr")[4]/ylab.step) * ylab.step
  ylab.aty <- seq( from =  ylab.min, to = ylab.max, by = ylab.step)
  for ( i.loop in 1:length( ylab.aty))
  {
    axis( at = ylab.aty[i.loop], labels =  fun.string.exp(ylab.aty[i.loop]), side = 2, cex = 0.75)
  }
}


#..Filenames----
#'@title fun.basename.DES
#'@name fun.basename.DES
#'@description Returns a filename specification from the values for DetectionProb, AnnualErradicationProb, SpreadReduction
#'@param d numeric DetectionProb
#'@param e numeric AnnualErradicationProb
#'@param s numeric SpreadReduction
#'@param short logical = T
#'@param norm logical = F
fun.basename.DES <- function( d, e, s, short = T, norm = F)
{
  #dir = ImageDir
  #DetectionProb = 0.05
  #AnnualErradicationProb = 0.196
  #SpreadReduction = 0.998
  #short <- T
  
  s.DP <- formatC( x = plyr::round_any(x = d, accuracy=.001, f=floor)
                   , format = "f", width = 3, flag = "0", digits = 3, zero.print = NULL, drop0trailing = F, decimal.mark = "." )
  s.EP <- formatC( x = plyr::round_any(x = e, accuracy=.001, f=floor)
                   , format = "f", width = 3, flag = "0", digits = 3, zero.print = NULL, drop0trailing = F, decimal.mark = "." )
  s.SR <- formatC( x = plyr::round_any(x = s, accuracy=.001, f=floor)
                   , format = "f", width = 3, flag = "0", digits = 3, zero.print = NULL, drop0trailing = F, decimal.mark = "." )
  
  if ( ! short)
  {
    if ( norm)
    {
      #--- norm's original file name algorithm
      val.return <- paste0("DetProb_", d,"_EradProb_", plyr::round_any(x = e, accuracy=.01, f=floor) ,"_SpreadReduction_",s)
    } else
    {
      #--- a file sorted algorithm
      val.return <- paste0("DetProb_", s.DP,"_EradProb_",s.EP, "_SpreadReduction_",SpreadReduction)
    }
  } else
  {
    #--- a file sorted algorithm with short names
    val.return <- paste0("DP", s.DP,"_EP",s.EP, "_SR",s.SR)
  }
  
  return( val.return)
}

#'@title fun.map.coastline
#'@name fun.map.coastline
#'@author Robbie Price
#'@description Returns a path for the coastline vector layer (currently returns the MWLR basedata coastline) allowing for checks to be added at a later date
#'@param x character Fullpath to the specified Coastline vector layer
#'@return character
fun.map.coastline <- function(x = r"[N:\Projects\BaseData\NZ\Coastline\LINZ\v20200803\nz-coastlines-and-islands-polygons-topo-150k.shp]")
{
  fn <- x
  
  if ( ! file.exists( fn))
  {
    stop( "Cannot locate Coastline dataset", fn)       
  }
  return( x)
}

#'@title fun.map.hillshade
#'@name fun.map.hillshade
#'@author Robbie Price
#'@description Returns a path for the hillshade raster (currently returns the MWLR basedata coastline) allowing for checks to be added at a later date
#'@param x character Fullpath to the specified Hillshade raster
#'@return character 
fun.map.hillshade <- function(x = r"[N:\Projects\BaseData\NZ\Landform\Hillshade\15mDEM\vsun359degat45az\nzhillshade_nztm_nzgd49_coast_359_45.tif]")
{
  fn <- x
  
  if ( ! file.exists( fn))
  {
    stop( "Cannot locate Hillshade dataset", fn)       
  }
  return( x)
}
#fun.map.hillshade(x = "")


fun.map.region <- function( x = r"[N:\Projects\BaseData\NZ\Boundaries\Loc_body\Region\StatsNZ_v20210209\regional-council-2020-generalised.shp]")
{
  fn <- x
  
  if ( ! file.exists( fn))
  {
    stop( "Cannot locate Hillshade dataset", fn)       
  }
  return( x)
}


#..Map Parameters----
#'@title fun.map.breaks
#'@name fun.map.breaks
#'@author Robbie Price
#'@description Returns breaks for a invasion probability map
#'@param
#'@return numeric
fun.map.breaks <- function()
{
  val.return <- c(0, 0.01,0.05,0.1,0.3,0.5,0.7,0.9,1)
  #class(val.return)
  
  return( val.return)
}

#'@title fun.map.colors.tealtored
#'@name fun.map.colors.tealtored
#'@author Robbie Price
#'@description Returns a list of colours to use when plotting invasion probability based on the breaks supplied
#'@param breaks numeric See  fun.map.breaks
#'@return character
#'@examples
#'fun.map.colors.tealtored( breaks = fun.map.breaks(), dir = getwd())
fun.map.colors.tealtored <- function ( breaks = fun.map.breaks(), dir = NA)
{
  colr <-  colorRampPalette(colors =  c( rgb( 0.2, 0.4, 0.7),  rgb( 0.80, 0.85, 0.00),  rgb( 0.55, 0.0, 0.01)))

  cols.probability <- colr(n = length(breaks))
  cols.invasion <- c( rgb( red = 0.01, green = 0.3, blue = 0.9, alpha = 0.1), add.alpha( col = cols.probability[2:length(cols.probability)], alpha = 0.8))
  #scales::show_col(cols.invasion)
  #colorblindcheck::palette_check(x = cols.invasion, plot = T)
  
  if ( ! is.na(dir))
  {
    if ( dir.exists(dir))
    {
      fn.check <- paste0(dir, "\\", "colorblindcheck.png")
      png( filename = fn.check, width = 18, height = 8, units = "cm", res = 150, pointsize = 6)
      par(mfrow = c(1, 2))
      scales::show_col(cols.invasion)
      colorblindcheck::palette_check(x = cols.invasion, plot = T)
      dev.off()
    }
  }
  
  #class(cols.invasion)
  
  return( cols.invasion)
}


#..Additional functions----
#'@title add.alpha
#'@name add.alpha
#'@author Markus Gesmann
#'@description Adds transparency to a vector of colours
#'@param col vector One or more colours to add transparency to.
#'@param alpha numeric 0-1 value representing the opacity of the colour.
#'@references 
#'Markus Gesmann (Apr 30, 2013) How to change the alpha value of colours in R. Retrieved from https://magesblog.com/post/2013-04-30-how-to-change-alpha-value-of-colours-in/
#'https://www.magesblog.com/post/2013-04-30-how-to-change-alpha-value-of-colours-in/
add.alpha <- function(col, alpha=1)
{
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}


#'@title fun.string.exp
#'@name fun.string.exp
#'@description Returns an expression formatting a number in to -BBB--- formatted numbers
#'@param chr character Eastingt or Northing value to format for axis label
#'chr = "1550000"
fun.string.exp <- function( chr)
{
  
  l3 <- substr(chr, start = nchar(chr)-2, stop = nchar(chr))
  l2 <- substr(chr, start = nchar(chr)-5, stop = nchar(chr)-3)
  l1 <- substr(chr, start = 1, stop = nchar(chr)-6)
  
  #bquote( "" ^ .(l1) ~ .(l2)  ^ .(l3))
  bquote( ""*.(l1)*bold(.(l2))*.(l3))
}


#..Main Plot Function----
#'@title fun.plot.map.detectionprobability
#'@name fun.plot.map.detectionprobability
#'@author Robbie Price
#'@description Creates map saved as PNG file for the specified data.  This creates consistent map layout for all region plots (hopefully).
#'@param v Vect Vector lay with column InvasionProbs for plotting
#'@param ras.hillshade.crop Rast Raster lay of hillshade
#'@param vct.coast Vect Coastline to play
#'@param vct.region vect Vector layer for region boundary. Default = NA
#'@param cols.invasion vector list of colours to use for v plot
#'@param brks.invasion vector list of breakpoints to use for v plot
#'@param outlines logical plot the outlines of the farm boundaries in solid line (or transparent)
#'@param fn.plot character Filename for output PNG file
#'@param ext.map bbox Bounding box extent for the plot
#'@param Title character Title for map
#'@param legend.title = "Invasion probability"
#'@param region character Name of region to plot (or "")
#'@param year numeric The year number to plot = 0
fun.plot.map.detectionprobability <- function ( v, ras.hillshade.crop, vct.coast, vct.region = NA
                                                , cols.invasion, brks.invasion, outlines, fn.plot, ext.map
                                                , Title, region, year, legend.title ="Invasion probability" )
{
  while( length( dev.list()) > 0)
  {
    dev.off()
  }

  #..plot an overview map as png----
#if ( (! is.na( vct.region)) & year == 1)

  if ( (! is.na( vct.region)))
  {
    fn.plot.overview <- paste0( dirname( path = fn.plot), r"[\]", "overview_", region, ".png")
    
    if (( ! file.exists(fn.plot.overview)) | year == 1)  
    {
      png( filename = fn.plot.overview, width = 4.5, height = 6, units = "cm", pointsize = 1/72, res = 300)
      ext.ovr <- ext.NZTM()
      
      terra::plot(x = vct.coast
                , mar=c(0.1, 0.1, 0.1, 0.1)
                , xlim = c( ext.ovr$xmin, ext.ovr$xmax)
                , ylim = c( ext.ovr$ymin, ext.ovr$ymax)
                , axes = FALSE
                , border = rgb(0.6,0.6,0.6)
                , col = rgb(0.9,0.9,0.9) 
                , lwd = 1)
      
      terra::plot(x = vct.region, col = NA, border = rgb(0.1,0.1,0.1), lwd = 3, lty = 1, legend = F, add = T )
      dev.off()
    } else 
    { 
      #warning( "Reusing overview file for", region)
    } 
  }
  
  
  #Open PNG device to create map image
  png( filename = fn.plot, width = 16.5, height = 15, units = "cm", pointsize = 1/72, res = 300)
  
  #..Initial Plot----
  terra::plot(x = ras.hillshade.crop
              , xlim =  c(ext.map$xmin,ext.map$xmax), ylim = c(ext.map$ymin,ext.map$ymax)
              , main = Title, mar=c(2.1, 1.3, 3.6, 7.3)
              , yaxt="n",yaxt="n" 
              , col = grey(level = c(0:100)/100, alpha = 0.5)
              , legend =F, axes = FALSE
              , asp = 1)
  
 
  
  #--- get everything required about the plot dimentsion
  lst.par <- par()
  plotpar <- list()
  plotpar$mar <- list( b = lst.par$mar, l = lst.par$mar[2], t = lst.par$mar[3], r = lst.par$mar[4])
  plotpar$win <- list( xl = lst.par$usr[1], xr = lst.par$usr[2], yb = lst.par$usr[3], yt = lst.par$usr[4], xc = (lst.par$usr[1]+lst.par$usr[2])/2)
  plotpar$pag <- list( w = dev.size("cm")[1], h = dev.size("cm")[2])
  plotpar$cm <- (plotpar$win$xr - plotpar$win$xl) / plotpar$pag$w - ( plotpar$mar$l + plotpar$mar$r)

  #.. scale bar
  lst.scales <- c(1,2,5)
  scale.four <- (par("usr")[2] - par("usr")[1]) / 4
  scale.dist <- log10( scale.four)
  scale.round <- lst.scales[which.min(abs(lst.scales - scale.dist))]
  scale.step <-  as.numeric( scale.round * (10 **  ceiling(log10( scale.four)-1) ))

  #..Frame around the map as ploted
  rect( xleft = plotpar$w$xl, xright = plotpar$w$xr, ybottom = plotpar$w$yb, ytop = plotpar$w$yt)

  #.. plot scale bar 1cm below frame
  p.s.xl <-  plotpar$win$xc - scale.step/2
  p.s.xr <-  plotpar$win$xc + scale.step/2
  
  lines(   x = c( p.s.xl, p.s.xr)
         , y = c( plotpar$win$yb - (plotpar$cm * 1.1) , plotpar$win$yb - (plotpar$cm * 1.1)), col = "black", xpd = T )
  lines(   x = c( p.s.xl, p.s.xl)
           , y = c( plotpar$win$yb - (plotpar$cm * 1.1) , plotpar$win$yb - (plotpar$cm * 0.8)), col = "black", xpd = T )
  lines(   x = c( p.s.xr, p.s.xr)
           , y = c( plotpar$win$yb - (plotpar$cm * 1.1) , plotpar$win$yb - (plotpar$cm * 0.8)), col = "black", xpd = T )
  text( x = p.s.xl, y = plotpar$win$yb - (plotpar$cm *  0.4)
        , label = "0", xpd = T, cex = 0.75, adj = (c(00.5,0.5 )) )
  text( x = p.s.xr, y = plotpar$win$yb - (plotpar$cm *  0.4)
        , label = scale.step/1000 , xpd = T, cex = 0.75, adj = (c(0.5,0.5 )) )
  text( x = p.s.xr + (plotpar$cm *  0.2), y = plotpar$win$yb - (plotpar$cm *  (1.1+0.8)/2)
        , label = "km" , xpd = T, cex = 0.75, adj = (c(0,0.5 )) )
  
  #..Axis labels----
  # Autodected best label range for plot base on a minimum of 4 lables

  
  #..Some map meta information----
  #lab.x <- par("usr")[2] - (par("usr")[2]-par("usr")[1])/100
  #lab.y <- par("usr")[3] + (par("usr")[4]-par("usr")[3])/100
  #text( x = lab.x, y = lab.y, label = "EPSG:2193 (NZTM)", xpd=NA, adj = (c(1,0 )), cex = 0.75 )
  text( x = plotpar$win$xr, y = plotpar$win$yb - (plotpar$cm *  0.4), label = "EPSG:2193 (NZTM)", xpd=NA, adj = (c(1,0.5 )), cex = 0.75 )

  #..Coastline----
  terra::plot(x = vct.coast, col = NA, border = rgb(0.4,0.4,0.4), lwd = 0.1, legend = F, add = T )
  
  #..Region (if requested)----
  if ( ! is.na( vct.region))
  {
    terra::plot(x = vct.region, col = NA, border = rgb(0.1,0.1,0.1), lwd = 0.8, lty = "3292", legend = F, add = T )
    ####lab.x <- par("usr")[1] + (par("usr")[2]-par("usr")[1])/2
    ####lab.y <- par("usr")[4] - (par("usr")[4]-par("usr")[3])/100
    ####text( x = lab.x, y = lab.y, label = region, xpd=NA, adj = (c(0.5,1 )), cex = 1 )
    ras.ovr <- png::readPNG(fn.plot.overview)
    rasterImage(image = ras.ovr
                , xleft = plotpar$win$xr+ (plotpar$cm * 0.1), xright = plotpar$win$xr + (plotpar$cm * 4.6)
                , ybottom = plotpar$win$yb, ytop = plotpar$win$yb + (plotpar$cm * 6)
                , xpd = T)
    
  }
  
  #..Outlines----
  #terra::plot(x = v,y = "InvasionProb" , lwd = 0.25, col = cols.invasion, breaks=c(0.01,0.05,0.1,0.3,0.5,0.7,0.9,1), legend = F , add = T)
  if ( outlines)
  {
    terra::plot(  x = v, y = "InvasionProb" 
                  , lwd = 0.2, border =  grey(level = 0.8, alpha = 1)
                  , type="interval"
                  , col = cols.invasion, breaks=brks.invasion
                  , legend = F , add = T)
    if ( 1==2)
    { ###DEBUG PLOT
      dev.off()
      dev.new()
      terra::plot(  x = v, y = "InvasionProb" 
                    , lwd = 0.2, border = NA
                    , col = cols.invasion
                    , type="interval"
                    , breaks= brks.invasion
                    , legend = T, all_levels =T)
      
    }
  } else
  {
    names(v)
    terra::plot(x = v,y = "InvasionProb" , lwd = 0.2
                , border =  grey(level = 0.3, alpha = 1), col = NA, legend = F , add = T)
    terra::plot(x = v, y = "InvasionProb"
                , type="interval"
                , col = cols.invasion, breaks=brks.invasion
                , border = NA 
                , legend = F , add = T)
  }
  
  par(xpd=TRUE)
  legend( "topright"
          , inset=c(-0.22,0)
          , bg  = "white"
          , title = legend.title 
          , legend =rev(brks.invasion), fill =rev(cols.invasion)
          , cex = 0.7,  title.adj = 0.1, adj = 0.1, xjust = 0.5, yjust = 0
          , bty = "n"
          , xpd=TRUE) #, lty=1, lwd=3)
  
  dev.off()
}




#Plot----
######################################################-
###Declare function for mapping farm level mean probability of  
##infestation (across realisations) in each year
######################################################-
#'@title InfestationHeatMaps
#'@name InfestationHeatMaps
#'@author Robbie Price (original by Norm Mason)
#'@description Creates an animated GIF of a time series heat map for pest invasion
#'@param Nyears numeric Simulation duration
#'@param Nperm numeric Number of permutations per paramteter combination
#'@param ImageDir character
#'@param FarmPolygons sf POlygon dataset that represents the Farm entities
#'@param InvasionProbResults array Three dimensional array of 0-1 probabilities of invasion?????? (ASK NORM)
#'@param TitleStem character  base for map title (will be augmented with year value for each time-dtep)
#'@param FileNameStem character basename of output files (will be augumented with a suffix by year and file extension )
#'@param Coastline character Fullpath to the Coastline Polygon dataset (might also be line?).  Default is function returning location = fun.map.coastline(),
#'@param Hillshade character Fullpath to the RASTER to use as background.  Default is function call to return location = fun.map.hillshade(),
#'@param breaks vector Numeric vector of break values. Default is functon call returning values fun.map.breaks(),
#'@param col vector List of colours to use for the Heat map.  Default is a function call returning colours fun.map.colors.tealtored(breaks = fun.map.breaks()) ,
#'@param outlines = F logical If true the outlines are draw for the proporty boundaries on top of fill colours, else lines are drawn first then fill overthe top
#'@param spf = 1 numeric Time in seconds for each frame (seconds per frame)
#'@param region character Name (or 4-letter code) of region to map. Defualt "" will ignore all code relating to regions  #2023/02/03
#'@export
InfestationHeatMaps  <- function(
    Nyears,
    InvasionProbResults,
    ImageDir,
    FarmPolygons,
    #InvasionResults,
    TitleStem,
    FileNameStem,
    Coastline = fun.map.coastline(),
    Hillshade = fun.map.hillshade(),
    Regions = fun.map.region(),
    breaks = fun.map.breaks(),
    cols.invasion = fun.map.colors.tealtored(breaks = fun.map.breaks(), dir = NA) ,
    outlines = F,
    spf = 1,
    region = "" #2023/02/03
)
{ 
  
  #Nyears = 3
  #warning( "NYears set to ", Nyears, " for debugging")
  
  #Create the output directory
  ImageDir.curated <- file.path( ImageDir)
  
  dir.create(path = ImageDir.curated, showWarnings = F, recursive = T)

  ImageDir.png <- paste0( ImageDir.curated, "\\", FileNameStem)
  dir.create(path = ImageDir.png, showWarnings = F, recursive = T)
  
  fn.plots <- paste0(ImageDir.png, "\\", FileNameStem,"_year_", formatC( x = 1:Nyears, width = 3, flag = "0"),".png")
  fn.gif <- paste0(ImageDir.curated, "\\", FileNameStem, "_animation.gif")
  
  #--- Assign the input parameters to the internal variables already in use (hard to refactor safely)
  fn.hillshade <- Hillshade #### r"[N:\Projects\BaseData\NZ\Landform\Hillshade\15mDEM\vsun359degat45az\nzhillshade_nztm_nzgd49_coast_359_45.tif]"
  fn.coastline <- Coastline #### r"[N:\Projects\BaseData\NZ\Coastline\LINZ\v20200803\nz-coastlines-and-islands-polygons-topo-150k.shp]"
  brks.invasion <- breaks
  #cols.invasion <- col
  ###fps <- ifelse( 100 %% fps == 0, fps, { warning( "fps needs to be a factor of 100, changed to 2." ); 2})

  #--- load the datasets
  ras.hillshade <- terra::rast(x = fn.hillshade)
  vct.coast <- terra::vect(x = fn.coastline) 

  #deal with possible CRS issues (WARNING: may not work for all transfroamtions)----
  #crs.maps <- sf::st_crs( sf::st_crs( x = 2193)$proj4string)$proj4string   
  #crs.hill <- sf::st_crs( x = sf::st_crs( x = terra::crs(ras.hillshade))$proj4string )$proj4string
  #crs.farm <- sf::st_crs(FarmPolygons)$proj4string
  #crs.coas <- sf::st_crs(x = terra::crs(vct.coast))$proj4string
  #Even doing the above is not reliable so will always project :(
  crs.value <- 2193
  crs.strin <- paste0( "epsg:", crs.value)

  #--- project the farm polygons  
  sf.FarmPolygons.NZMT <- sf::st_transform(x = FarmPolygons, crs =crs.strin )

  #Calculate map extent either by dataFarm or region extent ---- 
  vct.region <- NA
  if ( region == "" || ! file.exists(Regions))
  {
    ext.MAP <- sf::st_bbox(sf.FarmPolygons.NZMT)

    ext.map <- ext.round(x = ext.MAP, f = 5000, b = 5000)
    ext.map2 <- ext.round(x = ext.MAP, f = 10000, b = 10000)
    
    ext.map <- ext.square( x = ext.map)
    ext.map2 <- ext.square( x = ext.map2)
  
  } else
  {
    #....Use regions to define extents----
    #region = "MARL"
    #--- using the basedara regions file
    #prm.fn.regions <- fun.map.region()
    vct.regions <- terra::vect(x = Regions)
    #terra::plot( vct.regions)
    
    #--- fin dtrhe name of the region using the code provuded
    region.REGC2020_1 <- string.findbycode( x = region, y = vct.regions$REGC2020_1)
    region.REGC2020_2 <- string.findbycode( x = region, y = vct.regions$REGC2020_2)
    

    vct.region <- terra::subset(x = vct.regions, vct.regions$REGC2020_1 == region.REGC2020_1)    
    #terra::plot( vct.region)
    
    bbox.region <- terra::ext( vct.region)

    ext.MAP <- sf::st_bbox( obj = c(xmin = terra::xmin(bbox.region)
                                     , xmax = terra::xmax(bbox.region)
                                     , ymin = terra::ymin(bbox.region)
                                     , ymax = terra::ymax(bbox.region))
                          )
    
    ext.map <- ext.round(x = ext.MAP, f = 5000, b = 5000)
    ext.map2 <- ext.round(x = ext.MAP, f = 10000, b = 10000)
    
    ext.map <- ext.square( x = ext.map)
    ext.map2 <- ext.square( x = ext.map2)
    
  }
  
  ext.crop <- terra::ext(ext.map2$xmin-10000,  ext.map2$xmax+10000, ext.map2$ymin-10000, ext.map2$ymax+
                           10000)
  
  #ras.aoi <- terra::rast(extent = ext.crop, res = terra::res(ras.hillshade)) #, crs = terra::crs(ras.hillshade))
  ras.hillshade.crop <- terra::crop( x = ras.hillshade, y = ext.crop) 
  
  #--- because we can't compare crs's very reliably we will need to project regardless (STUPID R)
  vct.coast <- terra::project(x = vct.coast, y =  crs.strin)
  ras.hillshade.crop <- terra::project(x = ras.hillshade.crop, y = crs.strin)
  
  #--- now get the sf equivalents (terra doesn't plot properly so need to go back to sf objects...)
  vsp <- as(sf.FarmPolygons.NZMT, "Spatial")
  v <- terra::vect(vsp)

  year = 1
 
  ##Loop through realisations ( Nyears = 1)----
  year <- 1
  
  #InvasionProbResults$Mean <- rowMeans(x = InvasionProbResults, na.rm = T, dims = 1)
  while( length( dev.list()) > 0)
  {
    dev.off()
  }
  

  
  for(year in 1:Nyears)
  {
    cat("\r", "Creating map ", FileNameStem, ": ", year, " of ", Nyears, "...")
    Title = paste0(region.REGC2020_1, "\r\n", TitleStem," for Year ",year)  
    InvasionProb = InvasionProbResults[,year]
    v$InvasionProb <- InvasionProb
    
    fn.plot <- fn.plots[year] ### paste0(FileNameStem,"_year_", formatC( x = year, width = 3, flag = "0"),"_RJP_TOR.png")
    nchar(fn.plot)

    #....Calls single plot function as now also used for Mean plot (below)----
    fun.plot.map.detectionprobability(   v = v
                                       , ras.hillshade.crop = ras.hillshade.crop
                                       , vct.coast = vct.coast
                                       , vct.region = vct.region
                                       , cols.invasion = cols.invasion
                                       , brks.invasion = brks.invasion
                                       , fn.plot = fn.plot
                                       , outlines = outlines
                                       , ext.map = ext.map
                                       , Title = Title
                                       , region = region.REGC2020_2
                                       , year = year)

  }
  
  #################################################-
  #Create GIFF
  #..Create the GIFF----
  cat("\r\n", "Creating GIFF for ", FileNameStem, "\r\n")
  ##open the images
  img_list <- lapply(X = fn.plots, FUN = magick::image_read)
  
  ## join the images together
  img_joined <- magick::image_join(img_list)
  
  ## animate at 2 frames per second
  img_animated <- magick::image_animate(image = img_joined, optimize = T, delay = spf*100)
  
  ## view animated image
  #img_animated
  
  ## save to disk
  magick::image_write(image = img_animated, path = fn.gif)
  
  ##################################################-
  #..now plot the mean for the years----
  cat("\r", "Creating Mean values for ", FileNameStem, ": ", year, " of ", Nyears, "...")
  Title = paste0(region.REGC2020_1, "\r\n", TitleStem," Mean for Years 1 to ",year)  
  InvasionProbMean = rowMeans( x = InvasionProbResults[,1:Nyears], na.rm = T)
  v$InvasionProb <- InvasionProbMean
  #terra::writeVector(v, filename = r"[R:\Projects\SL2116_MPI_SleeperPests\Analysis\Test\Test.shp]",overwrite=T, options="ENCODING=UTF-8")
  #names(v)
  #scales::show_col(cols.invasion)
  #Plot---
  #fn.plot <- fn.plots[year] ### paste0(FileNameStem,"_year_", formatC( x = year, width = 3, flag = "0"),"_RJP_TOR.png")
  fn.plot.mean <- paste0(ImageDir.png, "\\", FileNameStem,"_Mean_001_", formatC( x = Nyears, width = 3, flag = "0"),".png")
#cat(  "\r\n". "DEBUG 647:  ", lst.Regions[region.idx], "\r\n")
  fun.plot.map.detectionprobability(     v = v
                                       , ras.hillshade.crop = ras.hillshade.crop
                                       , vct.coast = vct.coast
                                       , vct.region = vct.region
                                       , cols.invasion = cols.invasion
                                       , brks.invasion = brks.invasion
                                       , fn.plot = fn.plot.mean
                                       , outlines = outlines
                                       , ext.map = ext.map
                                       , Title = Title
                                       , legend.title = paste0( "Mean proportion", "\r\n", "of years invaded")
                                       , region = region.REGC2020_2
                                       , year = 0)  
  
  #dev.off()
  #warnings()   
}

