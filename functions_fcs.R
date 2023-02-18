## This is a file with functions ##

## This function reads in the data from the .fcs-file and then returns the information as
## a dataframe that we can work with. For it to run, it needs the package flowCore, from
## bioconductor. We can then save the dataframes as a csv file, so we don't have to work
## with the fcs files. I have already run this, so no need to repeat it. I left it here so
## you know all the steps I took.
## Requires flowCore
getFCS <- function( fileName, path = ".", channels = "", tube = "" , changeH = TRUE ){
  
  ## get the raw data from the file with name 'fileName' that is located in 'path'
  fcsfile <- flowCore::read.FCS( if(path != "") file.path( path, fileName) else fileName) 
    
  ## extract the data into data.frame
  df <- dplyr::as_tibble( flowCore::exprs( fcsfile ) )
  
  ## assign all channels to be extracted, unless requested channels are given
  channels <- if(channels[1] == "") colnames(df) else channels
  
  # extract the data for the channels and write them as text
  df<- dplyr::select(df,  all_of( channels ) )
  
  # add a column to indicate sample/tube name
  df$tube <- if(tube == "") fileName else tube
  
  # if requested change all 'strange' signs to '_' using gsub
  if(changeH) names(df) <- gsub("([-/ \\|.()])","_", names(df))

  # R doesn't like variables starting with a number. If a variable does, then add "c_" in front of it
  names(df)[grepl("^[0-9]",names(df))]<- paste0("C_",names(df)[grepl("^[0-9]",names(df))])
  
  return(df)
}
  
  
### 
plot_oneGate <- function( df, x, y, gate, setX = "", setY = "" ){
  p1 <- df %>% 
    ggplot( aes_string(x = x, y = y ))   + # plot these channels
    geom_hex( bins = 50 , show.legend = FALSE )  +               # use colors to indicate number of counts
    geom_path(data = as.data.frame( gate ) %>% setNames( c( x, y ) ), size=0.2 ) +
    scale_fill_gradientn(colours = rainbow(6)) +
    theme( 
      axis.text = element_text( size = 5),
      axis.title = element_text(size=6)
    )
    
  
  if( length(setY) > 1 ) p1 <- p1 + scale_y_continuous(limits = setY)
  
  if( length(setX) > 1 ) p1 <- p1 + scale_x_continuous(limits = setX)
  
  return(p1)
}

### 
plot_noGate <- function( df, x, y, setX = "", setY = "", bins = 50 ){
  p1 <- df %>% 
    ggplot( aes_string(x = x, y = y ))   + # plot these channels
    geom_hex( bins = bins , show.legend = FALSE )  +               # use colors to indicate number of counts
    scale_fill_gradientn(colours = rainbow(6)) +
    theme( 
      axis.text = element_text( size = 5),
      axis.title = element_text(size=6)
    )
  
  
  if( length(setY) > 1 ) p1 <- p1 + scale_y_continuous(limits = setY)
  
  if( length(setX) > 1 ) p1 <- p1 + scale_x_continuous(limits = setX)
  
  return(p1)
}

### 
plot_twoGates <- function( df, x, y, gate1, gate2, setX = "", setY = "" ){
  p1 <- df %>% 
    ggplot( aes_string(x = x, y = y ))   + # plot these channels
    geom_hex( bins = 80 , show.legend = FALSE )  +               # use colors to indicate number of counts
    geom_path(data = as.data.frame( gate1 ) %>% setNames( c( x, y ) ), size=0.2  ) +
    geom_path(data = as.data.frame( gate2 ) %>% setNames( c( x, y ) ), size=0.2  ) +
    scale_fill_gradientn(colours = rainbow(6)) +
    theme( 
      axis.text = element_text( size = 5),
      axis.title = element_text(size=6)
    )
  
  if( length(setY) > 1 ) p1 <- p1 + scale_y_continuous(limits = setY)
  
  if( length(setX) > 1 ) p1 <- p1 + scale_x_continuous(limits = setX)
  
  return(p1)
}
 
gateCells_2d <- function( df, xCh, yCh, gate ){
  
  g_st <- st_polygon( x = list( gate )  )
  d_points <- st_as_sf( df, coords = c ( xCh, yCh ) )
  return( lengths( st_within( d_points, g_st ) ) == 1 )
}


biexpTrans <- function( x, a = 0.5, b = 1, c = 0.5, d = 1, f = 0, w = 0 ) { 
  return(a * exp( b * ( x - w ) ) - c * exp( -d * ( x - w ) ) + f )
  }

writeFCS <- function( fileName, overwrite = FALSE, path = ".", 
                      channels = "", tube = "" , changeH = TRUE ){
  if( !file.exists(gsub(".fcs$",replacement = ".csv", fileName)) | overwrite ){
  write_csv(
    file = gsub(".fcs",".csv", fileName),    # to the file in which we replaced .fcs for .csv
    x = format_numeric( getFCS( fileName, path = path, channels = channels, 
                                tube = tube , changeH = TRUE)  ) )
  } else print("Files already exists. Overwrite using writeFCS(fileName, overwrite = TRUE) ")
}

## change the format of column to numercic if it is numberic 
format_numeric <- function(x, ...) {
    numeric_cols <- vapply(x, is.numeric, logical(1))
    x[numeric_cols] <- lapply(x[numeric_cols], format, ...)
    x
}


# To find the peak we look in density dens$y which value is higher than the previous 
# one, and also higher that the next one:  this should be a local peak.
# Here we make a function does just that. It uses a complex lookup function (which etc) that 
# returns which value in the array dens$y is a local peak. It thens sorts the peaks, highest
# first, and then takes the value on the x axis, which should be the peak. 
# (Functions are bits of script we can use over and over again. We give it info, and it 
# gives us info back).
findPeak<- function (dens, tresh = 0, valley = FALSE, no = 1) {
  y = dens$y
  x = dens$x
  
  if(valley){
    pks <- which(diff(sign(diff(y, na.pad = FALSE)), na.pad = FALSE) > 0) + 2
  } else{
    pks <- which(diff(sign(diff(y, na.pad = FALSE)), na.pad = FALSE) < 0) + 2
  } 
  
  # keep only the peaks that are locally 'steep' enough (enough defined by treshhold 'tresh')
  if (tresh > 0 ) {
    pks <- pks[y[pks - 1] - y[pks] > tresh]
  }
  
  # sort the peaks, highest (or lowest) first
  if(valley){
    pks <- match( sort(y[pks],decreasing = F), y)
  } else{
    pks <- match( sort(y[pks],decreasing = T), y)
  } 
  if(length(pks) < no) no = length(pks)
  
  return( x[ pks[1:no] ] )
}

findValley <- function( d ){
  tot = length(d)
  least = ceiling(0.001*tot)
  if( tot< 10 ) return(3)
  rangeMin = min(sort(d)[least : (tot-least) ])
  rangeMax = max(sort(d)[least : (tot-least) ])
  #plot(density( d, bw = 0.05))
  #plot(density( d[d < (rangeMax)  & d> rangeMin], bw = 0.05))
  dens <- density( d[d < (rangeMax)  & d> rangeMin], bw = 0.05)
  twinPeaks <- findPeak( dens , no = 3 ) # get the three highest peaks
  
  p <- sort(twinPeaks) # sort the peaks by hights
  if(length(p)>2){  # if there are three peaks, see which ones are closest 
                    # together and keep the ones that are furthest apart.
   twinPeaks<- ifelse(c(p[2]-p[1]<p[3]-p[2],T), # check which are closest
                     ifelse(
                       c(dens$y[match(p[3],dens$x)]>0.1,T), # check if peak is high enough.
                       c(p[3],p[2]),
                       c(p[1],p[2])),
                     ifelse(c(dens$y[match(p[1],dens$x)]>0.1,T),c(p[1],p[2]),c(p[3],p[2])))
  } else twinPeaks <- p
  
  minPeak <- min(twinPeaks)
  maxPeak <- max(twinPeaks)
  # rangeMid <- (rangeMax + rangeMin)/2
  # 
  # minPeak <- ifelse( length( d[d < (rangeMid*1.2) & d> rangeMin])<10,
  #                     rangeMin,
  #                     findPeak( density( d[d < (rangeMid)  & d> rangeMin] ) ))
  # maxPeak <- ifelse( length( d[d < (rangeMid*0.8)  & d> rangeMin])<10,
  #                     rangeMax,
  #                     findPeak( density( d[d > (rangeMid) & d< rangeMax] ) ))
  # 
  if(length( d[d > minPeak & d< maxPeak])<10 ) return( 
    findPeak( density( d[d > rangeMin & d< rangeMax] ), 
              valley = TRUE ))
  
  Valley <- findPeak( density( d[d > minPeak & d< maxPeak] ), valley = TRUE )
  return( Valley )
}


plotFcsData <- function(df, channels = "", cutOf = "", timePlot = F, timelim = 100){ 
  require(tidyverse)
  
  theme_set(theme_minimal())
  
  leg <- theme(
    legend.title = element_text(color = "blue", size = 14),
    legend.text = element_text(color = "red", size = 10)
  )
  
  if( nrow( dplyr::filter(df, filt == 1)) < 250 ){ return( plot(x = 0, y = 0) ) }
  
  plots <- vector(mode = "list", 
                  length = 2 + 
                    if_else(sum(channels != "") > 0, 4* length(channels), 0) + 
                    if_else(sum(channels != "") > 1, 2, 0) + 
                    ifelse(timePlot,length(channels),0) )

  if(sum(channels != "") > 0){

      for(i in seq(1, length(channels) )){
        ## For each channel plot4 times, against SSC_A and distribution before and after filtering 
        plots[[  4*i - 3]] <- 
          df %>% 
          ggplot( aes_string(x = "SSC_A", y= channels[i]))   + # plot these channels
          geom_hex( bins = 50 , show.legend = FALSE) +               # use colors to indicate number of counts
          ylab(channels[i]) +
          scale_fill_gradientn(colours = rainbow(6)) +
          scale_y_continuous(limits = c(0,6)) +
          scale_x_continuous(limits = c(0000,200000))+ leg
        plots[[  4*i - 2]]<- 
          df %>% 
          dplyr::filter( filt == 1) %>%
          ggplot( aes_string(x = "SSC_A", y= channels[i]))   + # plot these channels
          geom_hex( bins = 50 , show.legend = FALSE) +               # use colors to indicate number of counts
          ylab(channels[i]) +
          scale_fill_gradientn(colours = rainbow(6)) +
          geom_hline( yintercept = cutOf[i], color = "blue") + # add line between wt and red
          scale_y_continuous(limits = c(0,6)) +
          scale_x_continuous(limits = c(0000,200000)) + leg
        plots[[  4*i - 1]] <-
          df %>%
          ggplot( aes_string(x = channels[i]))  + ##
          geom_density( fill = "black") +    ## make a density plot
          xlab(channels[i]) +
          scale_x_continuous(limits = c(0,6)) +
          geom_vline( xintercept = cutOf[i], color = "blue") # add line indicating below wt and above red
        plots[[  4*i ]] <-
          df %>%
          dplyr::filter( filt == 1) %>%
          ggplot( aes_string(x = channels[i]))  + ##
          geom_density( fill = "black") +    ## make a density plot
          xlab(channels[i]) +
          scale_x_continuous(limits = c(0,6)) +
          geom_vline( xintercept = cutOf[i], color = "blue") # add line below wt and above red
      }
    plots[[ 4*i + 1]] <-
      df %>%  
      ggplot( aes(x = SSC_A, y= FSC_A)) +   # plot these channels
      geom_hex( bins = 50 , show.legend = FALSE) +               # use colors to indicate number of counts
      scale_fill_gradientn(colours = rainbow(6)) +
      scale_y_continuous(limits = c(0000,200000)) + # make the x and y axis run from 0 to 250000
      scale_x_continuous(limits = c(0000,200000))
    plots[[ 4 * i + 2 ]]<-
      df %>%  
      dplyr::filter( filt == 1 ) %>% ## keep only the ones in our filter and give to next function
      ggplot( aes(x = SSC_A, y= FSC_A))   + # plot these channels
      geom_hex( bins = 50 , show.legend = FALSE) +               # use colors to indicate number of counts
      scale_fill_gradientn(colours = rainbow(6)) +
      scale_y_continuous(limits = c(0000,200000)) + # make the x and y axis run from 0 to 250000
      scale_x_continuous(limits = c(0000,200000))
    if(sum(channels != "") > 1){
      plots[[ 4 * i + 3 ]] <- 
        df %>% 
        ggplot( aes_string(x = channels[1], y= channels[2]))   + # plot these channels
        geom_hex( bins = 50 , show.legend = FALSE) +               # use colors to indicate number of counts
        ylab(channels[i]) +
        scale_fill_gradientn(colours = rainbow(6)) +
        scale_y_continuous(limits = c(0,6)) +
        scale_x_continuous(limits = c(0,6))+ leg
      plots[[ 4 * i + 4 ]]<- 
        df %>% 
        dplyr::filter( filt == 1) %>%
        ggplot( aes_string(x = channels[1], y= channels[2]))   + # plot these channels
        geom_hex( bins = 50 , show.legend = FALSE) +               # use colors to indicate number of counts
        ylab(channels[i]) +
        scale_fill_gradientn(colours = rainbow(6)) +
        scale_y_continuous(limits = c(0,6)) +
        scale_x_continuous(limits = c(0,6)) + leg
    }
      
      if(timePlot){
        for(i in seq(1, length(channels))){
          plots[[ length(plots)-(i-1) ]] <-
              df %>%  
              ggplot( aes_string(x = "Time", y= channels[i]) )   +
              geom_hex( bins = 50 ,show.legend = FALSE) +
              ylab(channels[i]) +   # use colors to indicate number of counts
              geom_vline( xintercept = timelim, color = "blue") +  
              scale_fill_gradientn(colours = rainbow(6)) + leg
          }
      }
    
  }
  #p <- do.call("arrangeGrob", c(plots, ncol= if_else(length(plots) > 4,4,2)) )
  
  # plot the figures together (for this we use grid.arrange from gridExtra)
  return(   do.call("arrangeGrob", c(plots, ncol= if_else(length(plots) > 4,4,2)) ) )

}

# remove spaces from a file name so it can be used in bash etc.
removeSymbolsFileName <- function( fileNames ){
  file.rename(fileNames, file.path(dirname(fileNames),gsub("([-/ \\|()])","_", basename(fileNames))))
}

