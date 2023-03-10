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

# List function files
# https://www.r-bloggers.com/2013/01/storing-a-function-in-a-separate-file-in-r/
#source("C:\\RStudio_working\\dispot_str1_v1.R")

nagheli_case1fcn <- function(numpts, wells, consts){
  
  min_rz <- 0
  max_rz <- 1
  step_rz <- (max_rz - min_rz) / numpts
  rD <- seq(min_rz, max_rz, by=step_rz)
  
  #### --------------------------------------------------------------------------------------
  
  ## Use a different function to find the stagnation points (findStg_case1.R)
  # xy_stg <- findStg_case1(wells, consts)
  
  # OR keep the code here (nagheli_case1fcn.R and nagheli_case2fcn.R should use the same approach)
  
  ## wrapper function that uses only one variable but calls dispot_str1_v1deriv that needs three variables
  ffd1 <- function(x){
    zD <- complex(real = x[1], imaginary = x[2]) # be careful if there is more than one row!
    ## use dispot_str1_v1deriv.R
    x_comp <- dispot_str1_v1deriv(wells, consts, zD)
    x_return <- c(Re(x_comp), Im(x_comp))
    return(x_return)
  }
  
  # maxIter <- 200 ## use a grid instead of random points
  numptsgrid <- 20
  step <- 1/(numptsgrid)
  xvect <- seq(step, 1, by=step)
  
  # xf <- matrix(0.0, maxIter, 2) #xf = zeros(maxIter,1);
  xf <- matrix(0.0, numptsgrid * numptsgrid, 2)
  colnames(xf) <- c("x","y")
  xfcount <- 1
  
  # set.seed(5)
  
  # for (i in 1:maxIter){
  for (i in 1:(numptsgrid)){
    for (ii in 1:(numptsgrid)){
      ## NOTE: THE FOLLOWING MAY NOT BE IDEAL FOR CASES WHERE OTHER QUADRANTS AROUND THE ORIGIN ARE IN PLAY
      # x0 <- c(runif(1, min = 0, max = 1), runif(1, min = 0, max = 1))
      x0 <- c(xvect[i], xvect[ii])
      
      xftrial <- nleqslv(x0,ffd1)
      
      if (abs(xftrial$x[1]) <= 2 && abs(xftrial$x[2]) <= 2){
        ## rounding: https://stackoverflow.com/questions/24694001/rounding-to-two-decimal-places-in-octave
        ## try rounding to 2 decimal places for comparison
        check <- round(100*xftrial$x)/100 == round(100*xf)/100
        
        if (sum(check) == 0){
          xf[xfcount, c(1,2)] <- cbind(xftrial$x[1],xftrial$x[2])
          xfcount <- xfcount + 1
        }
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
  
  #### --------------------------------------------------------------------------------------
    
  ### Prepare variables for discharge and potential functions
  
  zD_xy <- matrix(0.0, length(rD) * length(rD),2);
  colnames(zD_xy) <- c("x","y")
  zD_x3 <- matrix(0.0, length(rD), length(rD))
  zD_y3 <- matrix(0.0, length(rD), length(rD))
  
  counter <- 1;
  
  for (i in 1:length(rD)){
    for (ii in 1:length(rD)){
      zD_xy[counter,"x"] <- (i-1)/length(rD)
      zD_xy[counter,"y"] <- (ii-1)/length(rD)
      zD_x3[i,ii] <- zD_xy[counter,"x"]
      zD_y3[i,ii] <- zD_xy[counter,"y"]
  
      counter <- counter + 1
    }
  }
  
  phi0D <- matrix(0.0, nrow(zD_xy),1)
  
  ##### FIND DIMENSIONLESS DISCHARGE POTENTIAL --------------------------------------------------------------
    
  xpts <- cbind(zD_xy[,"x"], zD_xy[,"y"], phi0D)
  comOmegaD <- matrix(0.0, nrow(zD_xy),1)
  comOmegaD <- dispot_str1_v1(wells, consts, xpts)
  phiD <- Re(comOmegaD)
  
  ##### FIND STREAMLINES (STREAM FUNCTION) --------------------------------------------------------------
  psiD <- Im(comOmegaD)
  
  #### find streamline magnitude through the stagnation point
  #### revise
  # psiD_stg = imag(dispot_str1_v1(wells, consts, [xf(1,1), xf(2,1), 0]));
  
  #### 3D GRID -----------------------------------------------------------------------------------------
  
  phiD3 <- matrix(0.0, length(rD), length(rD))
  psiD3 <- matrix(0.0, length(rD), length(rD))
  counter <- 1;
  
  for (i in 1:length(rD)){
    for (ii in 1:length(rD)){
      phiD3[i,ii] <- phiD[counter]
      psiD3[i,ii] <- psiD[counter]
      counter <- counter + 1;
    }
  }
  
  ##### PLOT
  #### Compare with Nagheli et al 2020 Figure 6a or 6b
  # if (option == 1){
  #   png("Nagheli_et_al_2020_fig6a.png",height=1500,width=1500,res=300)  
  #   #pdf("Nagheli_et_al_2020_fig6a.pdf")  
  # }else{
  #   png("Nagheli_et_al_2020_fig6b.png",height=1500,width=1500,res=300)
  #   #pdf("Nagheli_et_al_2020_fig6b.pdf")
  # }
  # 
  # ptA <- c(0, 0)
  # ptB <- c(1, 0)
  # boundary <- rbind(ptA, ptB)
  # 
  # #numContoursDis = 10
  # #numContoursStr = 20
  # 
  # plot(boundary[,1], boundary[,2], type = "l", cex = 2, xlab="Dimensionless x-axis", ylab="Dimensionless y-axis")
  # 
  # #contour(zD_x3, zD_y3, phiD3, numContoursDis, 'b--'); # discharge potential
  # contour(zD_xy[1:(numpts+1),2], zD_xy[1:(numpts+1),2], phiD3, nlevels=10, col="blue")
  # #contour(zD_x3, zD_y3, abs(psiD3), numContoursStr, 'r-'); # streamlines
  # par(new=TRUE)
  # contour(zD_xy[1:(numpts+1),2], zD_xy[1:(numpts+1),2], abs(psiD3), nlevels=20, col="red")
  # ## contour(zD_x3, zD_y3, psiD3, [psiD_stg psiD_stg], 'k-'); # streamlines
  # 
  # extwells <- subset(wells, wells[,3]>0)
  # injwells <- subset(wells, wells[,3]<0)
  # 
  # points(extwells[,1], extwells[,2],pch=18, col="black")
  # points(injwells[,1], injwells[,2],pch=5)
  # points(x_stg, y_stg, pch=18, col="green")
  # arrowLen <- 0.2
  # arrows(0, 0, arrowLen * cos(beta), arrowLen * sin(beta))
  # 
  # dev.off()
  
  # print(counter)
  
  # plotdata <- data.frame(
  #   # zD_xyd = zD_xy,
  #   phiD3d = phiD3,
  #   psiD3d = psiD3
  # )
  
  plotdata <- list()
  plotdata[[1]] <- xy_stg # added this to the list - 21 June 2022
  # plotdata[[2]] <- zD_xy[1:(numpts+1),2]
  plotdata[[2]] <- cbind(zD_xy[1:(numpts+1),2], zD_xy[1:(numpts+1),2]) ## changed this on 12 Dec 2022 for consistency with case 2
  plotdata[[3]] <- phiD3
  plotdata[[4]] <- psiD3
  
  return(plotdata)
}