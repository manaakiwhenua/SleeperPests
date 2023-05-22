#'@title string.findbycode
#'@name string.findbycode
#'@author Robbie Price
#'@description Parse a series of names for the closest matching value to a code.  Code can be any subset of the letters in the names 
#'@param x character Code or abbreviation to match against a lkist of names to find best match
#'@param y vector list of names to test code agains
#'@param s numeric Starting point in names to begin search (defaults to 1)
string.findbycode <- function( x, y, s = 1)
{
  val.return <- NA
  lst.index <- 1:length(y)

  i = 2
  for ( i in 1:length(y))  
  {
    lst.index[i] <- string.matchinorder( x = x, y = y[i], s = 1)  
  }
  
  if (any(! is.na(lst.index)))
  {
    val.return <- y[which( lst.index == min( lst.index, na.rm = T))]
  }
  
  return( val.return) 
}

#'@title string.matchinorder
#'@name string.matchinorder
#'@author Robbie Price
#'@description Scores the matching of a code against a string.  Low scores are better matches
#' attempts to find the sequential location of each chaarcter in "X" within "Y"
#' score is the sum of position locations.  Where a character x[j] occurs BEFORE x[i] then the length of x is added
#' CAT in CATERPILLAR C(1)A(2)T(3) 1+2+3 = 6
#' CAT in A catcher C(3)A(4)T(5) = 12
#' CAT in Actor C(2)A(4+1)T(4+3) = 14
#' 
#'@param x character Code or abbreviation to match against a lkist of names to find best match
#'@param y vector list of names to test code agains
#'@param s numeric Starting point in names to begin search (defaults to 1)
string.matchinorder <- function(x, y, s = 1)
{
  y <- toupper(y)
  x <- toupper(x)

  lst.matches <- rep( x = 0, times = nchar(x))
  i = 1
  s.fails <- 0
  #s = 1
  # i = 2
  
  for( i in 1:nchar(x) )  
  {
    c <- substr( x = x, start = i, stop = i)
    #cat( c, "\r\n")
    
    lst.finds <- unlist( gregexpr(pattern = c, text = y))
    lst.finds.s <- lst.finds[ which( lst.finds >= s)]

    s <- ifelse(length( lst.finds.s) > 0, lst.finds.s[1]+1, s)

    #cat(  lst.finds.s, sep = ", ")
    #cat( "\r\n")
    m <- i
    while( (length( lst.finds.s) == 0) && ( m > 0))
    {
      lst.finds <- unlist( gregexpr(pattern = c, text = y))
      lst.finds.s <- lst.finds[ which( lst.finds >= 1)]
      m <- m-1
      s.fails <- s.fails + nchar(y)
    }

    lst.matches[i] <- lst.finds.s[1] + s.fails
  }
  
  return( sum( lst.matches))
}

fun.test.string.findbycode <- function() 
{
  prm.fn.regions <- r"[N:\Projects\BaseData\NZ\Boundaries\Loc_body\Region\StatsNZ_v20210209\regional-council-2020-generalised.shp]"
  
  vct.regions <- terra::vect(x = prm.fn.regions)
  
  `Regions = c("AUCK", "EBOP", "CANT", "OTAG", "GISB", "WAIK", "HBAY","MNWG", "WELL", "MARL", "STHL", "NRTH", "TNKI",
              "TASM", "NELS")
  
  x = "OTAG"
  y = vct.regions$REGC2020_1[11]
  s = 1
  for ( i in Regions)
  {
    cat( i)
    s.region <- string.findbycode(x = i, y = vct.regions$REGC2020_1, s = 1)
    cat( " ", s.region, "r\n")
  }
`}
