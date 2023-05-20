# nagheli_case2fcn.r
#
# Andrew J. Wiebe, 25 Apr 2023; modifies version from 3 Jan 2023, 6 Dec 2022
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
#
# 7 Sep 2022 - use points from 1/numpts to 2 in order to capture image well beyond permafrost boundary
#
# Parameters:
#    numpts = the number of grid cells along x and y directions
#    option = option from Nagheli et al. (2020)
#    wells(:,1) = real coordinate for each well location (real(zwD))
#    wells(:,2) = complex coordinate for each well location (imag(zwD))
#    wells(:,3) = pumping rate for each well (>0 for injection, < 0 for extraction)
#    consts(1,1) - r_wedge, i.e., the length of the wedge in the actual system
#    consts(2,1) - alpha, i.e., the angle of the triangular aquifer
#    consts(3,1) - beta, i.e., the angle between the ambient aquifer flow direction and the x-axis
#    consts(4,1) - K, i.e., the hydraulic conductivity of the aquifer
#    consts(5,1) - b, i.e., the thickness of the aquifer
#    consts(6,1) - q0, i.e., regional uniform flow per unit width
# 
# Returns a list:
#    plotdata[[1]] <- x and y coordinates of stagnation points
#    plotdata[[2]] <- a two-column vector of x and y values corresponding to a grid for the phi (equipotential) and psi (flowline) values
#    plotdata[[3]] <- a two-dimensional matrix grid of phi (equipotential) values
#    plotdata[[4]] <- a two-dimensional matrix grid of psi (flowline) values
#


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
  
  maxIter <- 200  # random points
  # numptsgrid <- 20 # grid version
  # step <- 1/(numptsgrid)  # grid version
  # xvect <- seq(step, 1, by=step)  # grid version
  
  xf <- matrix(0.0, maxIter, 2)  # random points
  # xf <- matrix(0.0, numptsgrid * numptsgrid, 2) # grid version
  colnames(xf) <- c("x","y")
  xfcount <- 1
  
  set.seed(5) # random points
  
  for (i in 1:maxIter){ # random points
  # for (i in 1:(numptsgrid)){  # grid version
    # for (ii in 1:(numptsgrid)){  # grid version
      ## NOTE: THE FOLLOWING MAY NOT BE IDEAL FOR CASES WHERE OTHER QUADRANTS AROUND THE ORIGIN ARE IN PLAY
      x0 <- c(runif(1, min = 0, max = 1), runif(1, min = 0, max = 1)) # random points
      # x0 <- c(xvect[i], xvect[ii]) # grid version
      
      checkffd <- ffd(x0)
      
      if (!is.infinite(checkffd[1]) && !is.infinite(checkffd[2]) && !is.na(checkffd[1]) && !is.na(checkffd[2])){
        xftrial <- nleqslv(x0,ffd)
        
        # print(xftrial)
        
        if (abs(xftrial$x[1]) <= 2 && abs(xftrial$x[2]) <= 2){
          ## rounding: https://stackoverflow.com/questions/24694001/rounding-to-two-decimal-places-in-octave
          ## try rounding to 2 decimal places for comparison
          # check <- round(100*xftrial$x)/100 == round(100*xf)/100 # numbers are near zero
          # check <- round(1000*xftrial$x)/1000 == round(1000*xf)/1000
          # 
          # if (sum(check) == 0){
          #   xf[xfcount, c(1,2)] <- cbind(xftrial$x[1],xftrial$x[2])
          #   xfcount <- xfcount + 1
          # }
          
          ### updated from test23.R, 25 Apr 2023
          check1 <- round(100*xftrial$x[1])/100 == round(100*xf[,1])/100
          check2 <- round(100*xftrial$x[2])/100 == round(100*xf[,2])/100
          
          
          checklist1 <- which(check1 == TRUE)
          checklist2 <- which(check2 == TRUE)
          
          accept <- TRUE
          
          if(length(checklist1) > 0 && length(checklist2) > 0){
            for(i2 in 1:length(checklist1)){
              for(i3 in 1:length(checklist2)){
                if(checklist1[i2] == checklist2[i3]){
                  accept <- FALSE
                }
              }
            }
          }
          
          if(accept == TRUE){
            xf[xfcount, c(1,2)] <- cbind(xftrial$x[1],xftrial$x[2])
            xfcount <- xfcount + 1
          }
        }
      }
    # } # grid version
  }
  
  xf = xf[c(1:xfcount - 1), c(1,2), drop=FALSE] # xf = xf[c(1:xfcount - 1), c(1,2)]
  
  # print(xf)
  
  ## Note: if maxIter is high enough (e.g., 50 for five wells), all the stagnation points are found
  xy_stg <- matrix(0.0, (xfcount-1),2)
  
  x_stg <- xf[,"x"]
  y_stg <- xf[,"y"]
  xy_stg[,1] <- xf[,"x"]
  xy_stg[,2]<- xf[,"y"]
  
  # print(nrow(xy_stg))
  
  #### NOTE: AT THIS POINT, THE POINTS SHOULD BE ANALYZED AND ONLY THOSE WITHIN THE AQUIFER SHOULD BE SELECTED
  #### IDENTIFY STAGNATION POINTS LOCATED WITHIN AQUIFER (PUT THEM IN xy_stg2) ------
  xy_stg2 <- matrix(0.0,length(x_stg),2)
  counter <- 1
  
  for (i in 1:length(x_stg)){
  	### if radius < 1 and angle <= alpha, then OK
  	disttemp <- sqrt(x_stg[i]^2 + y_stg[i]^2)
  	thetatemp <- atan(y_stg[i] / x_stg[i])
  	
  	# print(thetatemp)
  	# print(y_stg[i])
  	# print(x_stg[i])
  	
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
  
  # print("nrow")
  # print(nrow(xy_stg2))
  
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
  if (option == 1 || option == 2 || option == 3){
    # zDxy_sqr <- cbind((1:numpts) / numpts, (1:numpts) / numpts)  ### should this go from zero to numpts?????????????
    zDxy_sqr <- cbind(seq(0,1,1/(numpts)), seq(0,1,(1/(numpts)))) # try this instead, for consistency with the option 4 line below
    xpts_sqr <- matrix(0.0,((numpts + 1) * (numpts + 1)),3) # revise this from below
    numpts_rev <- numpts + 1
  }else{
    # zDxy_sqr <- cbind(seq(-1,1,2/(numpts-1)), seq(0,2,(2/(numpts-1))))
    zDxy_sqr <- cbind(seq(-1,1,1/(numpts)), seq(0,2,(1/(numpts))))
    xpts_sqr <- matrix(0.0,((2 * numpts + 1) * (2 * numpts + 1)),3) # put this here instead of below
    numpts_rev <- 2 * numpts + 1
  }
  
  # str(zDxy_sqr)
  
  counter <- 1
  
  # xpts_sqr <- matrix(0.0,((2 * numpts + 1) * (2 * numpts + 1)),3)

  # for (i in 1:(2*numpts+1)){
  # 	for (ii in 1:(2*numpts+1)){
  for(i in 1:numpts_rev){
    for(ii in 1:numpts_rev){
  		xpts_sqr[counter,1] <- zDxy_sqr[i,1]
  		xpts_sqr[counter,2] <- zDxy_sqr[ii,2]
  		xpts_sqr[counter,3] <- 0 # phi0D
  		counter <- counter + 1
  	}
  }
  
  # comOmegaD_sqr <- matrix(0.0,((2*numpts+1)*(2*numpts+1)),1) # change this
  comOmegaD_sqr <- matrix(0.0,(numpts_rev * numpts_rev),1) # changed 3 Jan 2023
  comOmegaD_sqr <- dispot_str2_v4(wells, consts, xpts_sqr)
  phiD_sqr <- Re(comOmegaD_sqr)
  psiD_sqr <- Im(comOmegaD_sqr)
  
  xy_phiD <- cbind(xpts_sqr[,1], xpts_sqr[,2], phiD_sqr[,1])
  
  # phiD_sqr3 <- matrix(0.0,2*numpts+1, 2*numpts+1)
  # psiD_sqr3 <- matrix(0.0,2*numpts+1, 2*numpts+1)
  phiD_sqr3 <- matrix(0.0,numpts_rev, numpts_rev) # 3 Jan 2023
  psiD_sqr3 <- matrix(0.0,numpts_rev, numpts_rev) # 3 Jan 2023
  
  counter <- 1
  
  # for (i in 1:(2*numpts+1)){ # 3 Jan 2023
  	# for (ii in 1:(2*numpts+1)){ # 3 Jan 2023
  	for (i in 1:numpts_rev){
  	 for (ii in 1:numpts_rev){
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
