##################################################-
###Simulate spread from infested farms in Blind River
###Initial invasion approximates infestation documented here:
###https://www.tandfonline.com/doi/abs/10.1080/0028825X.1989.10414122
###Apply management as 1000-fold reduction in dispersal 
###And local erradication following 
###https://www.landcareresearch.co.nz/uploads/public/Discover-Our-Research/Biosecurity/Biocontrol-ecology-of-weeds/3-applications/CNG-review-SFF-project-2010.pdf
###Calls a core function for fixed initial infestations
##################################################-

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

###########################################################################-
###########################################################################-
###Read in fixed start function
###########################################################################-
###########################################################################-

#Source----
if ( ! exists("main.dir"))
{
  source(r"[R:\Projects\SL2116_MPI_SleeperPests\Analysis\Chilean needle grass\Updated Analyses 2022\CurrentFutureClimate\A01_MainlineDefaultSettings_Parameters.R]")
  #source(r"[R:\Projects\SL2116_MPI_SleeperPests\Analysis\Chilean needle grass\Updated Analyses 2022\CurrentFutureClimate\Mainline Default settings.R]")
  
}
if ( ! exists("INApestFixedStart"))
{
  source(r"[R:\Projects\SL2116_MPI_SleeperPests\Analysis\Chilean needle grass\Updated Analyses 2022\CurrentFutureClimate\INApestFixedStart.R]")
}

source(r"[T:\Hamilton\Projects P-T\Spatlab\Development\SL1714_RLibrary_smap-dsm\R\raster-extents.R]")




#Functions----

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
      val.return <- paste0("DetProb_", d,"_ErradProb_", plyr::round_any(x = e, accuracy=.01, f=floor) ,"_SpreadReduction_",s)
    } else
    {
      #--- a file sorted algorithm
      val.return <- paste0("DetProb_", s.DP,"_ErradProb_",s.EP, "_SpreadReduction_",SpreadReduction)
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
  scales::show_col(cols.invasion)
  colorblindcheck::palette_check(x = cols.invasion, plot = T)
  
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


#..Plot----
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
#'@param InvasionResults array Three dimensional array of 0-1 probabilities of invasion?????? (ASK NORM)
#'@param TitleStem character  base for map title (will be augmented with year value for each time-dtep)
#'@param FileNameStem character basename of output files (will be augumented with a suffix by year and file extension )
#'@param Coastline character Fullpath to the Coastline Polygon dataset (might also be line?).  Default is function returning location = fun.map.coastline(),
#'@param Hillshade character Fullpath to the RASTER to use as background.  Default is function call to return location = fun.map.hillshade(),
#'@param breaks vector Numeric vector of break values. Default is functon call returning values fun.map.breaks(),
#'@param col vector List of colours to use for the Heat map.  Default is a function call returning colours fun.map.colors.tealtored(breaks = fun.map.breaks()) ,
#'@param outlines = F logical If true the outlines are draw for the proporty boundaries on top of fill colours, else lines are drawn first then fill overthe top
#'@param spf = 1 numeric Time in seconds for each frame (seconds per frame)
#'@export
InfestationHeatMaps  <- function(
    Nyears,
    InvasionProbResults,
    ImageDir,
    FarmPolygons,
    InvasionResults,
    TitleStem,
    FileNameStem,
    Coastline = fun.map.coastline(),
    Hillshade = fun.map.hillshade(),
    breaks = fun.map.breaks(),
    col  = fun.map.colors.tealtored(breaks = fun.map.breaks(), dir = ImageDir) ,
    outlines = F,
    spf = 1
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
  cols.invasion <- col
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

  #--- prject the farm polygons  
  sf.FarmPolygons.NZMT <- sf::st_transform(x = FarmPolygons, crs =crs.strin )
  ext.NZTM <- sf::st_bbox(sf.FarmPolygons.NZMT)

  #Calculate map extent---- 
  ext.map <- ext.round(x = ext.NZTM, f = 5000, b = 5000)
  ext.map2 <- ext.round(x = ext.NZTM, f = 10000, b = 10000)
  
  ext.map <- ext.square( x = ext.map)
  ext.map2 <- ext.square( x = ext.map2)
  
  ext.crop <- terra::ext(ext.map2$xmin,  ext.map2$xmax, ext.map2$ymin, ext.map2$ymax)
  #ras.aoi <- terra::rast(extent = ext.crop, res = terra::res(ras.hillshade)) #, crs = terra::crs(ras.hillshade))
  ras.hillshade.crop <- terra::crop( x = ras.hillshade, y = ext.crop) 
  
  #--- because we can't compare crs's very reliably we will need to project regardless (STUPID R)
  vct.coast <- terra::project(x = vct.coast, y =  crs.strin)
  ras.hillshade.crop <- terra::project(x = ras.hillshade.crop, y = crs.strin)
  
  #--- now get the sf equivalents (terra doesn't plot properly so need to go back to sf objects...)
  vsp <- as(sf.FarmPolygons.NZMT, "Spatial")
  v <- terra::vect(vsp)

  year = 1
 
  ##Loop through realisations
  year <- 1
  for(year in 1:Nyears)
  {
    cat("\r", "Creating map ", FileNameStem, ": ", year, " of ", Nyears, "...")
    Title = paste0(TitleStem,"\nYear ",year)  
    InvasionProb = InvasionProbResults[,year]
    
    v$InvasionProb <- InvasionProb
    #terra::writeVector(v, filename = r"[R:\Projects\SL2116_MPI_SleeperPests\Analysis\Test\Test.shp]",overwrite=T, options="ENCODING=UTF-8")
    
    #scales::show_col(cols.invasion)
    #Plot---
    fn.plot <- fn.plots[year] ### paste0(FileNameStem,"_year_", formatC( x = year, width = 3, flag = "0"),"_RJP_TOR.png")
    nchar(fn.plot)
    
    png( filename = fn.plot, width = 16, height = 14, units = "cm", pointsize = 1/72, res = 300)
    #
    #dev.new()
    terra::plot(x = ras.hillshade.crop
                , xlim =  c(ext.map$xmin,ext.map$xmax), ylim = c(ext.map$ymin,ext.map$ymax)
                , main = Title, mar=c(2.1, 2.1, 4.1, 6.1)
                , yaxt="n",yaxt="n"
                , col = grey(level = c(0:100)/100, alpha = 0.5), legend =F, asp = 1)
    
    #terra::plot(x =ras.hillshade.crop, col = grey(level = c(0:100)/100, alpha = 0.5), legend =F, add = T)
    
    terra::plot(x = vct.coast, col = NA, border = rgb(0.4,0.4,0.4), lwd = 0.1, legend = F, add = T )
    
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
            , title = " Invasion probability"
            , legend =rev(brks.invasion), fill =rev(cols.invasion)
            , cex = 0.6,  title.adj = 0.1, adj = 0.1, xjust = 0, yjust = 0) #, lty=1, lwd=3)
    
    dev.off()

  }
  
  #Create GIFF
  #....Create the GIFF----
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
    
}

#Test Code----
fun.test.infestationHeatMaps <- function()
{
  warning("This function is for development and debugging only.", call. = T)
  # THIS CODE IS REQUIRED FOR DEBUGGING RUNS ONLY
  ######################################################-
  ######################################################-
  ###Estimate annual erradication probability required to 
  ###achieve 0.95 prob of achieving erradication
  ###following existing model of seedbank decline
  ###with management (using glyphosate) to prevent seed production
  ###use this as management efficay in INAscene function call
  ######################################################-
  ######################################################-
  ###Seeds per metre square
  Nseeds = 1500
  ###Seed loss from seed bank due to attrition
  SeedBankLoss = 0.61
  ###Seed loss from seed bank due to germination
  Germination = 0.064
  
  ###Model reduction in seed bank
  Nyears = 50
  Nseeds = vector(length = Nyears)
  N0 = 1500
  for(i in 1:Nyears)
  {
    Nseeds[i] = N0*(1-(SeedBankLoss+Germination))
    N0 = Nseeds[i] 
  }
  ###Find year where seed density <1 seed per hectare
  ###Consider this effective erradiation 
  ErradicationYear = length(Nseeds[Nseeds>1/10^4])
  
  ###Annual erradication prob to achieve 0.95 errodication
  ###after 14 (year when Nseeds <1 per hectare) years
  AnnualErradicationProb = 1-(1-0.95)^(1/ErradicationYear)
  
  #############################################################-
  #############################################################-

  ResultsDir= paste0(main.dir,"/HistoricExamplesCoreFunction/",ClimateScenarios[cs],"/")
  dir.create(ResultsDir,recursive = T)
  Regions = c("AUCK", "EBOP", "CANT", "OTAG", "GISB", "WAIK", "HBAY","MNWG", "WELL", "MARL", "STHL", "NRTH", "TNKI",
              "TASM", "NELS")
  region = 10
  RegionResultsDir= paste0(ResultsDir,Regions[region],"/")
  dir.create(RegionResultsDir)
  
  ########################################################-
  ###Read in simple XY points for INA
  ###Ecoclimatic index values and climate based extablishment for each farm
  ########################################################-
  
  data<-read_csv(paste0(ClimateScenarios[cs],"/Inputs/",Regions[region],"_Simple_points_for_INA.csv"))
  
  #Points for plotting in INA
  geocoords<-matrix(c(data$X, data$Y), byrow=F, ncol=2)
  
  ########################################################-
  ###Calculate distance-based dispersal probability
  ###and set long distance dispersal rate per farm
  #############################################################-
  
  ###Set mean annual number of incoming long distance events per farm year
  ###This ist the number of events is expected in a fully stocked matrix
  ###i.e all farms infested. Number of events will increase as more farms are infested
  LongDistProbPerFarm = 0.05
  
  ###Read in distance matrix
  dist_mat = readRDS(paste0(ClimateScenarios[cs],"/Inputs/",Regions[region],"_dist_mat_farm.rds"))
  
  ###Use custom distance kernel to obtain probabilities
  library(actuar)
  k= 7 #shape parameter
  b= 500 #scale parameter
  thresh= 1/20 #1 in 20 year spread risk between adjacent farms
  adj = thresh*actuar::ppareto(dist_mat, shape=k,  scale=b, lower.tail=F,log.p = F)
  
  ###Store names for nodes in network
  FarmNames = row.names(adj)

  ##############################################################-
  ###Assign climate based establishment prob based on ecoclimatic index
  ###Default function is logit
  ##############################################################-
  prob_est<-as.vector(data$Probability_Estab)
  if(EI_Prob_CurveType == "SplitLinear")
    prob_est<-as.vector(data$Probability_Estab_SplitLinear)
  
  ####################################################################-
  ###Read in historic invasion info
  ####################################################################-
  InputDir = r"[R:\Projects\SL2116_MPI_SleeperPests\Analysis\Chilean needle grass\Updated Analyses 2022\CurrentFutureClimate\Inputs\]"
  HistoricInfestation = read.csv(paste0(InputDir,"BlindRiverInfestedFarms.txt"),header = T, as.is =T)

  ####################################################-
  ###Input paramaters for INA
  ####################################################-
  ###Set number of realisations
  Nperm = 100
  ###Set simulation duration
  Nyears = 80
  
  ###Probability of management adoption when infestation detected
  ManageProb = 0.9
  
  ###Set spread reduction under management
  ###Assume 1000-fold reduction in spread prob under management
  SpreadReduction = 0.999
  
  ###Management efficacy defined as annual erradication probability
  ManEfficacy = AnnualErradicationProb
  
  ###Set probability of detection
  ###50% farmers confident of CNG ID 
  ###Default assumes those with ID knowledge have 50% chance of detection annually
  ###Declare vector to explore other probs
  DetectionProbs = c(0,0.05,0.1,0.15,0.2,0.25,0.3) 
  
  
  ###Declare initial infestation vector
  ###Start with infested farms from historic data
  InitBio = rep(0,times = nrow(adj))
  InitBio[FarmNames %in% HistoricInfestation$farm_id] = 1 
  
  ###Read in farm polygons
  RegionFarms = sf::st_read(paste0(ClimateScenarios[cs],"/Inputs/",Regions[region],"SheepBeefAtRisk",ClimateScenarios[cs], ".wgs84.shp"))
  ResultsDir= paste0(main.dir,"\\HistoricExamplesCoreFunction\\",ClimateScenarios[cs])
  dir.create(ResultsDir,recursive = T)
  Regions = c("AUCK", "EBOP", "CANT", "OTAG", "GISB", "WAIK", "HBAY","MNWG", "WELL", "MARL", "STHL", "NRTH", "TNKI",
              "TASM", "NELS")
  region = 10
  RegionResultsDir= paste0(ResultsDir, "\\", Regions[region])
  
  
  ImageDir = paste0(RegionResultsDir,"\\", "InvasionHeatMaps_V03")
  dir.create(path = ImageDir, showWarnings = F, recursive = T)
  #}
  detprob <- 1
  for(detprob in 1:length(DetectionProbs))
  {
    DetectionProb = DetectionProbs[detprob]
    TitleStem = paste0("Detection Prob. ",DetectionProb," Erradication Prob. ",round(AnnualErradicationProb,digits =3),
                       "\nSpread Reduction ",SpreadReduction)
    InvasionFileNameStem = paste0(RegionResultsDir, "\\", "DetProb_", DetectionProb,"_ErradProb_",round(AnnualErradicationProb,digits =2),
                                  "_SpreadReduction_",SpreadReduction)
    ##InvasionFileNameStem = paste0(RegionResultsDir, "\\", "DP_", DetectionProb,"_EP_",round(AnnualErradicationProb,digits =2),
    ##                              "_SR_",SpreadReduction)
    InvasionFileNameStem = paste0(RegionResultsDir, "\\"
                                   , fun.basename.DES( d = DetectionProb, e = AnnualErradicationProb, s = SpreadReduction, short = F, norm = T ))
    
    InvasionFileName <- paste0(InvasionFileNameStem,"InvasionProb.rds")
    file.exists(InvasionFileName)
    
    InvasionProbResults = readRDS(InvasionFileName)
    
    #ImageFileNameStem = paste0(ImageDir,"DetProb_",DetectionProb,"_ErradProb_",round(AnnualErradicationProb,digits =2),
    #                           "_SpreadReduction_",SpreadReduction,"_")
    #ImageFileNameStem = paste0("DetProb_",DetectionProb,"_ErradProb_",round(AnnualErradicationProb,digits =2), "_SpreadReduction_",SpreadReduction)
    ###ImageFileNameStem = paste0("DP_",DetectionProb,"_EP_",round(AnnualErradicationProb,digits =2), "_SR_",SpreadReduction)
                               
    ImageFileNameStem <- fun.basename.DES( d = DetectionProb, e = AnnualErradicationProb, s = SpreadReduction, short = T,norm = F )
    
    class(InvasionResults)
    
    #..Call----
    InfestationHeatMaps(
        Nyears = Nyears,
        InvasionProbResults = InvasionProbResults,
        ImageDir = ImageDir,
        FarmPolygons = RegionFarms,
        TitleStem = TitleStem,
        FileNameStem = ImageFileNameStem
        , Coastline = fun.map.coastline()
        , Hillshade = fun.map.hillshade()
        , breaks = fun.map.breaks()
        , col = fun.map.colors.tealtored(breaks = fun.map.breaks(), dir = ImageDir)
        , spf = 0.75
      )

    #break
  }
}


#fun.test.infestationHeatMaps()