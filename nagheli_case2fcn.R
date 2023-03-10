# nagheli_case2fcn.r
#
# Andrew J. Wiebe, 6 Dec 2022; modifies version from 3 June 2022, 6 Oct 2022
# 
# Use the Nagheli et al. (2020) method to find the stagnation point for a triangular aquifer (case 2)
#    with two stream boundaries and one semi-infinite boundary, and with multiple pumping wells.
# 
# NOTES (BEFORE CAP ZONE PTS EXTRACTION):
#    Qw is positive for extraction (pumping) and negative for injection (seems counter-intuitive)
#    dispot_str2_v3deriv_imagX(wells, consts, x) works for one well with fsolve
#    Both dispot_str2_v3deriv_imagX(wells, consts, x) and dispot_str2_v4deriv_imagX(wells, consts, x) work for two wells (Nagheli et al 2020 Fig 7b) with fsolve under certain x0 conditions
#    When there are multiple possibilities, need to weed out possible bad estimates (e.g., real(xf) < 0 and imag(xf) < 0 for alpha < pi/2)
#    Is there a way to estimate a priori the number of expected stagnation points based on the number of wells (and whether pumping or injection?)
#    Works OK (not perfectly for option 3 / fig7c)
#    Needed to add negative versions of the discharge pot contour magnitudes to draw the rest of the capture zone boundaries (Some of the discharge potential lines for the capture zone boundaries (through the stagnation points) get cut off at the cutoff lines)
#    Need to deal with case where stagnation point outside the aquifer (e.g., Whati case with no pf) and the entire aquifer is within the capture zone
#
# Notes (cap zone pts extraction)
#    option 1 / fig 7a: Works well with numpts = 200, threshAng = 1.0 rad - 9:45 am, 19 May 2022 (only issue: the set of points for the cap zone does not reach the x-axis -- need to extrapolate, or use more points for calculating the streamlines?
#    option 2 / fig 7b: Works well with numpts = 300; has a spurious line through the lower well for numpts = 200
#    option 3 / fig 7c: Almost works for numpts = 200 with modifications.
#                       Splitting curves works better when using two slopes over four points (outer two; skip the middle one) instead of three.
#                       Problem: extraction of points from contour lines fails to find half of capture zone through right stagnation point
#                       
# Based on the Octave code (Octave-5.1.0.0 Script) in nagheli_case2test4.m, 17 May 2022; modifies nagheli_case2test3.m, 16 May 2022, modifies nagheli_v11v2.m, 9 May 2022; and test_case2stagnatv7.m, 11 May 2022
#
# Provide an option to remove stagnation points outside the aquifer
#
# In an effort to evade branch cuts, return a matrix with columns for x,y,phi, then use raster stuff with slope and aspect to draw vectors
#
# 7 Sep 2022 - use points from 1/numpts to 2 in order to capture image well beyond permafrost boundary

# library(shiny)
# library(stats)
# library(nleqslv)
# library(signal) # for phase unwrap

# List function files
# https://www.r-bloggers.com/2013/01/storing-a-function-in-a-separate-file-in-r/
# source("C:\\RStudio_working\\dispot_str2_v4.R")
# source("C:\\RStudio_working\\dispot_str2_v3deriv_imagX.R")

nagheli_case2fcn <- function(numpts, option, wells, consts){

  # numpts <- 400 # one less than the number of points for rz and theta; try a larger number, based on correspondence with Dr. Nagheli.
  #### careful: note: numpts changes later for rectangular grid
  
  min_rz <- 0
  max_rz <- 1
  # max_rz <- 2 # try 2 here to get points out near image well
  step_rz <- (max_rz - min_rz) / numpts
  rD <- seq(min_rz, max_rz, by=step_rz)
  
  alpha <- consts[2]
  
  #### Wrapper function to return complex coordinates --------------------------------------------------------
  ffd <- function(x){
    zD <- complex(real = x[1], imaginary = x[2])
    zD
    x_comp <- dispot_str2_v3deriv_imagX(wells, consts, zD)
    x_comp
    x_return <- c(Re(x_comp), Im(x_comp))
    return(x_return)
  }
  
  #### Find Stagnation Points --------------------------------------------------------
  
  maxIter <- 200 ## use a grid instead of random points
  # numptsgrid <- 20
  # step <- 1/(numptsgrid)
  # xvect <- seq(step, 1, by=step)
  
  xf <- matrix(0.0, maxIter, 2) #xf = zeros(maxIter,1);
  # xf <- matrix(0.0, numptsgrid * numptsgrid, 2)
  colnames(xf) <- c("x","y")
  xfcount <- 1
  
  # set.seed(5)
  
  for (i in 1:maxIter){
  # for (i in 1:(numptsgrid)){
    # for (ii in 1:(numptsgrid)){
      ## NOTE: THE FOLLOWING MAY NOT BE IDEAL FOR CASES WHERE OTHER QUADRANTS AROUND THE ORIGIN ARE IN PLAY
      x0 <- c(runif(1, min = 0, max = 1), runif(1, min = 0, max = 1))
      # x0 <- c(xvect[i], xvect[ii])
      
      checkffd <- ffd(x0)
      
      if (!is.infinite(checkffd[1]) && !is.infinite(checkffd[2]) && !is.na(checkffd[1]) && !is.na(checkffd[2])){
        xftrial <- nleqslv(x0,ffd)
        
        # print(xftrial)
        
        if (abs(xftrial$x[1]) <= 2 && abs(xftrial$x[2]) <= 2){
          ## rounding: https://stackoverflow.com/questions/24694001/rounding-to-two-decimal-places-in-octave
          ## try rounding to 2 decimal places for comparison
          # check <- round(100*xftrial$x)/100 == round(100*xf)/100 # numbers are near zero
          check <- round(1000*xftrial$x)/1000 == round(1000*xf)/1000
          
          if (sum(check) == 0){
            xf[xfcount, c(1,2)] <- cbind(xftrial$x[1],xftrial$x[2])
            xfcount <- xfcount + 1
          }
        }
      }
    # }
  }
  
  xf = xf[c(1:xfcount - 1), c(1,2)]
  
  
  ## Note: if maxIter is high enough (e.g., 50 for five wells), all the stagnation points are found
  xy_stg <- matrix(0.0, (xfcount-1),2)
  
  x_stg <- xf[,"x"]
  y_stg <- xf[,"y"]
  xy_stg[,1] <- xf[,"x"]
  xy_stg[,2]<- xf[,"y"]
  
  # print(xy_stg)
  
  #### NOTE: AT THIS POINT, THE POINTS SHOULD BE ANALYZED AND ONLY THOSE WITHIN THE AQUIFER SHOULD BE SELECTED
  #### IDENTIFY STAGNATION POINTS LOCATED WITHIN AQUIFER (PUT THEM IN xy_stg2) ------
  xy_stg2 <- matrix(0.0,length(x_stg),2)
  counter <- 1
  
  for (i in 1:length(x_stg)){
  	### if radius < 1 and angle <= alpha, then OK
  	disttemp <- sqrt(x_stg[i]^2 + y_stg[i]^2)
  	thetatemp <- atan(y_stg[i] / x_stg[i])
  	
  	if (thetatemp < 0 && y_stg[i] > 0 &&  x_stg[i] < 0){ # NW quadrant
  	  thetatemp <- pi + thetatemp
  	}else if(thetatemp > 0 && y_stg[i] < 0 &&  x_stg[i] < 0){ # SW quadrant
  	  thetatemp <- pi + thetatemp
  	}else if (thetatemp < 0 && y_stg[i] < 0 && x_stg[i] > 0){ # SE quadrant
  	  thetatemp <- 2*pi + thetatemp
  	}
  	
  	if(thetatemp <= alpha){
  		xy_stg2[counter,] <- cbind(x_stg[i], y_stg[i])
  		counter <- counter + 1
  	}
  }
  
  xy_stg2 <- xy_stg2[c(1:counter-1), c(1,2),drop=FALSE] # remove zero rows, this time, prevent the matrix from becoming a vector if the number of rows drops down to 1 (https://stackoverflow.com/questions/61756409/why-does-nrow-return-a-null-value)
  
  #### find discharge potential magnitude through the stagnation point ------------------
  
  psiD_stg <- matrix(0.0, 2 * nrow(xy_stg2), 1)
  
  for (i in 1:nrow(xy_stg2)){
    x_temp <- cbind(xy_stg2[i,1], xy_stg2[i,2], 0)
    psiD_stg[i] <- Im(dispot_str2_v4(wells, consts, x_temp))
  }
  
  for (i in (nrow(xy_stg2)+1):(2*nrow(xy_stg2))){ # fixed this, but now nrow(indxs2) causing error # add the negative versions of the magnitudes for the other side of the cut lines? Seems to work
    psiD_stg[i] <- - psiD_stg[i - nrow(xy_stg2)]
  }
  
  
  #### --------------------------------------------------------------------------------------
  
  ##### CALCULATE EQUIPOTENTIAL LINES WITH RECTANGULAR GRID FOR FINDING CAPTURE ZONE BOUNDARIES ----
  
  # zDxy_sqr <- cbind((1:numpts) / numpts, (1:numpts) / numpts)  ### should this go from zero to numpts?
  # zDxy_sqr <- cbind((1:(2*numpts)) / numpts, (1:(2*numpts)) / numpts) # try going from 1/numpts up to 2
  
  ### from nagheli_case2test4.R
  if (option == 1 || option ==2 || option == 3){
    zDxy_sqr <- cbind((1:numpts) / numpts, (1:numpts) / numpts)  ### should this go from zero to numpts?????????????
  }else{
    # zDxy_sqr <- cbind(seq(-1,1,2/(numpts-1)), seq(0,2,(2/(numpts-1))))
    zDxy_sqr <- cbind(seq(-1,1,1/(numpts)), seq(0,2,(1/(numpts))))
  }
  
  # str(zDxy_sqr)
  
  counter <- 1
  
  xpts_sqr <- matrix(0.0,((2 * numpts + 1) * (2 * numpts + 1)),3)

  for (i in 1:(2*numpts+1)){
  	for (ii in 1:(2*numpts+1)){
  		xpts_sqr[counter,1] <- zDxy_sqr[i,1]
  		xpts_sqr[counter,2] <- zDxy_sqr[ii,2]
  		xpts_sqr[counter,3] <- 0 # phi0D
  		counter <- counter + 1
  	}
  }
  
  comOmegaD_sqr <- matrix(0.0,((2*numpts+1)*(2*numpts+1)),1)
  comOmegaD_sqr <- dispot_str2_v4(wells, consts, xpts_sqr)
  phiD_sqr <- Re(comOmegaD_sqr)
  psiD_sqr <- Im(comOmegaD_sqr)
  
  xy_phiD <- cbind(xpts_sqr[,1], xpts_sqr[,2], phiD_sqr[,1])
  
  phiD_sqr3 <- matrix(0.0,2*numpts+1, 2*numpts+1)
  psiD_sqr3 <- matrix(0.0,2*numpts+1, 2*numpts+1)
  
  counter <- 1
  
  for (i in 1:(2*numpts+1)){
  	for (ii in 1:(2*numpts+1)){
  	  phiD_sqr3[i,ii] <- phiD_sqr[counter]
  	  psiD_sqr3[i,ii] <- psiD_sqr[counter]
  		counter <- counter + 1
  	}
  }
  
  plotdata <- list()
  plotdata[[1]] <- xy_stg2
  plotdata[[2]] <- zDxy_sqr
  plotdata[[3]] <- phiD_sqr3
  plotdata[[4]] <- psiD_sqr3

  return(plotdata)
}
