# Code used to generate output used in the paper,
#  "Diminishing CO2-driven gains in water use efficiency of global forests",
#  by Mark A Adams, Thomas N Buckley and Tarryn L Turnbull
# Created by Tom Buckley, 2019

##load libraries
library(doBy)
library(dplyr)
library(data.table)
library(reshape2)
library(plyr)

options("scipen"=100, "digits"=4)

##insert code to set working directory

#####################################################
#####################################################
#####################################################
#### First part of code converts input data (input.csv) to processed outputs,
#### and is included here for the sake of reproducibility. Users wishing to 
#### evaluate the processed iWUE (W) outputs can skip to the second part of the
#### code, below the line saving the file "processed.csv".
#####################################################
#####################################################
#####################################################



#load in collated input data
x <- read.table("input.csv", header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)

##partition based on data type (will recombine later after calculation of common Delta)
xu <- subset(x, datatype=="u") #uncorreced d13Cplant
xc <- subset(x, datatype=="c") #corrected d13C
xt <- subset(x, datatype=="t") #'triangle' Delta (capital Delta)
xw <- subset(x, datatype=="w") #iwUE
xi <- subset(x, datatype=="i") #ci/ca



######################################################
######################################################
##### deprocessing data
##### y=uncorrected d13Cplant
#####     (1) compute Delta

xu$Delta <- (xu$d13Catm - xu$y)/(1+0.001*xu$y)


##### y=corrected d13Cplant:
#####     (1) de-correct using d13Catm
#####     (2) compute Delta from decorrected value

xc$y = xc$y + xc$d13Catm + 6.4
xc$Delta <- (xc$d13Catm - xc$y)/(1+0.001*xc$y)


##### y=WUEi:
#####     (1) back out Delta as a + (b-a)*(1 - 1.6*WUEi/ca)

a = 4.4
xw$b <- 27

#study #34 used b = 29
z <- with(xw, study=="034") 
xw$b <- replace(xw$b, z, 29)

xw$Delta = a + (xw$b - a)*(1 - 1.6*xw$y/xw$ca)

xw <- subset(xw, select = -c(b)) #remove b col, no longer needed


##### y=ci:
#####     (1) back out Delta as a + (b-a)*ci/ca

xi$Delta = a + (27 - a)*xi$y/xi$ca 


##### y=Delta:
#####     (1) keep as is

xt$Delta = xt$y


######################################################
######################################################
#####
##### //done deprocessing data (all now Delta on a common basis)
##### //recombine and save

##recombine dfs 
x <- Reduce(function(x,y) merge(x,y,all=TRUE), list(xu, xc, xt, xw, xi))

##save deprocessed data (converted to common basis) to a file
write.csv(x, "converted.csv")



######################################################
######################################################
#####
##### calculations for iWUE



########
## (1) estimate A/ca needed for 'full' isotope model
##    assume doubling ratio (DR) = 1.45 from Keeling et al 2017

doublingratio=1.45
beta=log(doublingratio)/log(2)-1
A_ca280 = 9/280 #assume pre-industrial A=9 
x$A_ca = A_ca280*((x$ca/280)^beta)


########
## (2) set parameters for full isotope model

gm = 0.2 #mesophyll conductance (mol m-2 s-1)
f = 12 #photorespiration discrimination term (per mille)
gammastar = 43 #photorespiratory CO2 compensation point (ppm)
a = 4.4 #diffusive discrimination (per mille)
b = 30 #Rubisco discrimination (per mille)
am = 1.8 # (per mille)

########
##(3) compute ci/ca from Delta, and iWUE from ci/ca

##DR=1.45
x$ci_ca = (x$Delta - a + (b - am)*x$A_ca/gm + f*gammastar/x$ca)/(b - a)
x$iWUE = x$ca*(1-x$ci_ca)/1.6

##using simple isotope model
b_simple = 27
x$ci_ca_simple = (x$Delta - a)/(b_simple - a)
x$iWUE_simple = x$ca*(1-x$ci_ca_simple)/1.6

#save results
write.csv(x, file="processed.csv")




#####################################################
#####################################################
#####################################################
#### Second part of code, below, calculates rates of change of W
#### with respect to ca for chosen subsets of the overall dataset.
#####################################################
#####################################################
#####################################################




######################################################
######################################################
#####
##### calculations for dW/dca

#load in collated input data
x <- read.table("processed.csv", header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)

#create new column containing first and last year of each dataset
x$minyear <- ave(x$year, x$tree, FUN=min)
x$maxyear <- ave(x$year, x$tree, FUN=max)

#####
#####
## function get_dWdca(): 
##   for a given subset of full dataset (with given starting and ending years,
##   and optionally including only isotope series that span a certain period),
##   compute diWUE(W)/dca for each data series and store the results.
### NOTE: omit args spani and spanf to use all data rather than only series spanning a certain range
#
get_dWdca <- function(xin, startyear, endyear, spani, spanf) {

  x <- xin  

  #subset dataset for spanning constraints
  if(missing(spani)) {
    spani="alldata"; spanf=""
  }else{
    x <- subset(x, x$minyear<=spani & x$maxyear>=spanf)
  }

  #subset for period of interest
  x <- subset(x, x$year>=startyear & x$year<=endyear)
  
  ##omit rows with NAs
  x <- na.omit(x)
  
  ##get group-specific values of lifeform, etc, within the subsetted period 
  o=summaryBy(lifeform + group + lat + lon + minyear + maxyear ~ tree + site + study + author, data=x, FUN=mean)
  colnames(o)[5]<-"lifeform"
  colnames(o)[6]<-"group"
  colnames(o)[7]<-"lat"
  colnames(o)[9]<-"lon"
  colnames(o)[9]<-"minyear"
  colnames(o)[10]<-"maxyear"
  o$abslat=abs(round(o$lat,1)) #absolute value of latitude
  
  o$group[o$group==1] <- "angiosperm"
  o$group[o$group==2] <- "gymnosperm"
  o$lifeform[o$lifeform==1] <- "unknown"
  o$lifeform[o$lifeform==2] <- "deciduous"
  o$lifeform[o$lifeform==3] <- "evergreen"

  #create blank columns for dW/dca results
  o$dWdca <- NA
  o$intercept <- NA
  o$r2 <- NA
  o$p <- NA
  
  o$dWdca_simple <- NA
  o$r2_simple <- NA
  o$p_simple <- NA
  
  
  ########
  ## run regression on iWUE vs ca for each isotope series ("tree")
  #   to get dW/dca, record intercept, r^2 and p
  
  for(s in unique(x$tree)) {
    xsub=subset(x, tree==s)
    
    #model each version of WUEi calculation
    # detailed model with doubling ratio = 1.45
    m=lm(iWUE ~ ca, data=xsub)
    o$dWdca[[which(o$tree==s)]]= m$coefficients[2]
    o$intercept[[which(o$tree==s)]]= m$coefficients[1]
    o$r2[[which(o$tree==s)]]= summary(m)$r.squared
    o$p[[which(o$tree==s)]]= summary(m)$coefficients[,4][2]
    
    # simple isotope model
    m=lm(iWUE_simple ~ ca, data=xsub)
    o$dWdca_simple[[which(o$tree==s)]]= m$coefficients[2]
    o$r2_simple[[which(o$tree==s)]]= summary(m)$r.squared
    o$p_simple[[which(o$tree==s)]]= summary(m)$coefficients[,4][2]
    
  }
  
  #store output
  fname=paste("_dWdca_span_", spani, "_", spanf,
              "_years_", startyear, "_", endyear,
              ".csv", sep="")
  print(fname)
  write.csv(o, file=fname)
}


###################
###################
###################
### sample code to generate output for different subsets

### first decade of 20th century, using only series that span the 20th century
get_dWdca(x,1901,1910,1901,2000)

### first decade of 20th century, using all available data for that period
get_dWdca(x,1901,1910)

