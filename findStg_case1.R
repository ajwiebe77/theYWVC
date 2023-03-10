# nagheli_case1test1.m
# 
# Andrew J. Wiebe, 3 June 2022
# 
# Objective: Find stagnation points, discharge potential, and streamlines for case 1 of Nagheli et al 2020 method.
# 
# Notes:
#    This version assumes confined conditions.
#    This version assumes a homogeneous aquifer
#    Dr. Nagheli advises a small step for dx and dy (e.g., 0.001) in order to have the flow lines go through the well point.
#    The fig6b case does not look exactly the same - maybe more than 200 points along each direction is not enough? There is one extra stagnation point
#    There are cut lines in both figures, with an offset in the flowlines at the cutlines
#
# Based on Octave code from nagheli_case1test1.m, 12 May 2022; modifies nagheli_v11v2.m, 9 May 2022
# Coded for R 1 June 2022
#
# Notes:
# Problem: There is the possibility that the code will not find all stagnation points (a grid search might be better than random numbers)
# There are some equipotential lines along the bottom of the map
# The stagnation points are not exactly the same as in the Nagheli et al 2020 paper... increase the number of points?

# library(stats)
# library(nleqslv)

# FIND STAGNATION POINTS --------------------------------------------------------------

findStg_case1 <- function(wells, consts){

  ## wrapper function that uses only one variable but calls dispot_str1_v1deriv that needs three variables
  ffd1 <- function(x){
    zD <- complex(real = x[1], imaginary = x[2]) # be careful if there is more than one row!
    ## use dispot_str1_v1deriv.R
    x_comp <- dispot_str1_v1deriv(wells, consts, zD)
    x_return <- c(Re(x_comp), Im(x_comp))
    return(x_return)
  }
  
  maxIter <- 200 ## FUTURE: use a grid instead of random points
  
  xf <- matrix(0.0, maxIter, 2) #xf = zeros(maxIter,1);
  colnames(xf) <- c("x","y")
  xfcount <- 1
  
  set.seed(5)
  
  for (i in 1:maxIter){
    ## NOTE: THE FOLLOWING MAY NOT BE IDEAL FOR CASES WHERE OTHER QUADRANTS AROUND THE ORIGIN ARE IN PLAY
    x0 <- c(runif(1, min = 0, max = 1), runif(1, min = 0, max = 1))
    #xftrial <- fsolve(ffd1, x0, tol=0.01)#tol = .Machine$double.eps^(0.5))
    xftrial <- nleqslv(x0,ffd1)
    
    if (abs(xftrial$x[1]) <= 2 && abs(xftrial$x[2]) <= 2){
      ## rounding: https://stackoverflow.com/questions/24694001/rounding-to-two-decimal-places-in-octave
      ## try rounding to 2 decimal places for comparison
      check <- round(100*xftrial$x)/100 == round(100*xf)/100
      
      if (sum(check) == 0){
        #xf[xfcount, c(1,2)] <- xftrial$x # works but can't read the elements as expected
        xf[xfcount, c(1,2)] <- cbind(xftrial$x[1],xftrial$x[2])
        xfcount <- xfcount + 1
      }
    }
  }
  
  xf = xf[c(1:xfcount - 1), c(1,2)]
  
  ## Note: if maxIter is high enough (e.g., 50 for five wells), all the stagnation points are found
  xy_stg <- matrix(0.0, (xfcount-1),2)
  
  # x_stg <- xf[,"x"]
  # y_stg <- xf[,"y"]
  xy_stg[,1] <- xf[,"x"]
  xy_stg[,2]<- xf[,"y"]

  return(xy_stg)
}