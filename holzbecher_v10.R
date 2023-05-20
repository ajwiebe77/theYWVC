### holzbecher_v10.R
### Andrew J. Wiebe, 17 May 2023
#
# Objective: Remove branch cuts from contour lines by using the sequential contouring method by Holzbecher (2018)
#
# 17 May 2023 - works for Whitehorse example (one flowline near the stagnation point is incomplete)
# 17 May 2023 - verification: works for case 2 option 4 (one segment seems to be missing outside the boundary around (-0.5,0.4), but looks correct within the aquifer
# 17 May 2023 - verification: works for case 2 option 3
# 17 May 2023 - verification: works for case 2 option 2
# 17 May 2023 - verification: works for case 2 option 1
 
# Reference:
#    Holzbecher, E. 2018.Streamline visualization of potential flow with branch cuts, with applications to groundwater. J. Flow Vis. Image Process. 25(2):119-144 (2018). DOI: 10.1615/JFlowVisImageProc.2018025918.
#
# Revisions since v1:
#    Fixed code for contourLines so that negative x values would be accommodated (index calculations)
#    in xsequence commands: changed the zD_xy index for x
#    Works if the wells are not immediately adjacent (i.e., within 7 or so cells in the matrix)
#    v4: added code to ignore stagnation points outside the aquifer
#    fixed problem with some stagnation points outside domain (x > maxx, y > maxy)
#    Added a way to deal with inf or nan values in psi based on https://stackoverflow.com/questions/35974779/finding-the-max-of-a-r-dataframe-column-ignoring-inf-and-na
#    Modified so that case 2 option 1 works (keep single row matrices as matrices); deal with sequential contouring for the case of only one well
#    Modified so that case 2 option 2 works adjust maxz update using a temporary variable since part of a matrix must be selected
#    
#    
# Parameters:
#    wells(:,1) = real coordinate for each well location (real(zwD))
#    wells(:,2) = complex coordinate for each well location (imag(zwD))
#    wells(:,3) = pumping rate for each well (>0 for injection, < 0 for extraction)
#    consts(1,1) = r, i.e., the length of the constant head boundary
#    consts(2,1) = beta, i.e., the angle between the ambient aquifer flow direction and the x-axis
#    consts(3,1) = K, i.e., the hydraulic conductivity of the aquifer
#    consts(4,1) = b, i.e., the thickness of the aquifer
#    consts(5,1) = q0, i.e., regional uniform flow per unit width
#    x(:,1) = zDx, list of dimensionless x coordinates at which to evaluate the equation
#    x(:,2) = zDy, list of dimensionless y coordinates at which to evaluate the equation
#    x(:,3) = phiD0 - the initial dimensionless discharge potential with no pumping, for each (real(zD)=zDx,imag(zD)=zDy) pair
# 
# Return values:
#    returnList[[1]] <- psi3, a list of the streamlines that avoids branch cuts, as a set of x and y points for each
#    returnList[[2]] <- capzones, a list of the capture zone boundary lines, as a set of points for each

holzbecher_v10 <- function(psi, zD_xy, wells, consts, xy_stg){
  
  distthreshstg <- 0.01
  
  if(length(consts) == 5){ # 12 Dec 2022 - add this check to address both cases 1 and 2
    r <- consts[1]
    K <- consts[3]
    b <- consts[4]
  }else{
    r <- consts[1]
    K <- consts[4]
    b <- consts[5]
  }
  
  wells2 <- wells
  wells2[,3] <- wells2[,3] / (K * b * r) # dimensionless Q values
  
  wells3 <- wells2[which(wells2[,3] > 0),,drop=FALSE] # 11 Apr 2023 - pumping wells
  wells4 <- wells2[which(wells2[,3] < 0),,drop=FALSE] # 19 Apr 2023 - injection wells
  
  numContoursStr <- 30 # 15
  
  N <- nrow(wells2)
  n <- nrow(psi)
  
  minx <- min(zD_xy[,1])
  maxx <- max(zD_xy[,1])
  xstep <- ((maxx - minx)/(nrow(zD_xy) - 1))
  
  miny <- min(zD_xy[,2])
  maxy <- max(zD_xy[,2])
  ystep <- ((maxy - miny)/(nrow(zD_xy) - 1))
  
  ### 12 Dec 2022 - add a check for xy_stg so that only stagnation points within the aquifer are addressed (e.g., points outside box for case 1):
  xy_stg_aqf <- xy_stg
  checkx <- which(xy_stg[,1] < minx)

  if(length(checkx) > 0){
    xy_stg_aqf <- xy_stg[-checkx,, drop=FALSE] # https://stackoverflow.com/questions/21977894/conditional-deleting-rows-in-a-matrix-in-r
  }
  
  checkx2 <- which(xy_stg[,1] > maxx) # add a check of the maximum values, too - 19 Apr 2023
  
  if(length(checkx2) > 0){
    xy_stg_aqf <- xy_stg[-checkx2,, drop=FALSE] # https://stackoverflow.com/questions/21977894/conditional-deleting-rows-in-a-matrix-in-r
  }
  
  checky <- which(xy_stg_aqf[,2] < miny)

  if(length(checky) > 0){
    xy_stg_aqf <- xy_stg_aqf[-checky,,drop=FALSE]
  }
  
  checky2 <- which(xy_stg_aqf[,2] > maxy) # add a check of the maximum values, too - 19 Apr 2023
  
  if(length(checky2) > 0){
    xy_stg_aqf <- xy_stg_aqf[-checky2,,drop=FALSE]
  }
  
  ### 25 Jan 2023 Add a check to ensure that stagnation points are not too close to each other
  if(nrow(xy_stg_aqf) > 1){
    checkdist <- matrix(0, nrow(xy_stg_aqf), nrow(xy_stg_aqf))
    skiprows <- matrix(0,1,nrow(xy_stg_aqf))
    xy_stg_aqf_temp <- xy_stg_aqf
    for(i in 1:(nrow(xy_stg_aqf_temp) - 1)){
      if(skiprows[1,i] == 0){
        for(ii in (i+1):nrow(xy_stg_aqf_temp)){
          diststg <- sqrt((xy_stg_aqf[i,1] - xy_stg_aqf[ii,1])^2 + (xy_stg_aqf[i,2] - xy_stg_aqf[ii,2])^2)
          if(diststg < distthreshstg){
            xy_stg_aqf <- xy_stg_aqf[-ii,,drop=FALSE] # remove duplicate
            skiprows[1,ii] <- 1
          }
        }
      }
    }
  }
  
  num_stg_pts <- nrow(xy_stg_aqf)
  xy_stg <- matrix(xy_stg_aqf, num_stg_pts,2)
  
  # Holzbecher 2018 Algorithm, adapted for direction (branch cuts go down instead of left)
  # Divide domain into blocks based on the well locations
  ybreaks <- sort(wells2[,2], decreasing=FALSE)

  if(nrow(wells2) > 1){ # add this 3 Jan 2023
    wells2 <- wells2[order(wells2[,2], decreasing=FALSE),] # problem here for Whati (case 2 opt 1)
  }

  xbreaks <- wells2[,1]
  
  # Method:
  # Set capphi[N+1] <- capphi
  # for region with y < ybreaks[1], do nothing (i.e., no changes to psi)
  
  # Draw stream function contours for capphi for region zy > zN (i.e., unaffected region)
  # change the contour interval in the different sections this time
  
  minz <- min(psi[is.finite(psi)]) # https://stackoverflow.com/questions/35974779/finding-the-max-of-a-r-dataframe-column-ignoring-inf-and-na
  maxz <- max(psi[is.finite(psi)])
  
  spacing <- (maxz - minz) / numContoursStr
  
  if(nrow(wells2) > 1){ # add this 3 Jan 2023
    nspace <- ceiling(wells2[which(wells2[,2] == max(wells2[,2])),3] / spacing) * 2
    nspace1 <- nspace[1] # just in case there is a tie for the maximum y value # 27 Jan 2023
    spacing <- wells2[which(wells2[,2] == max(wells2[,2])),3] / nspace
    spacing <- spacing[1] # just in case there is a tie for the maximum y value # 27 Jan 2023
  }else{ # add this 3 Jan 2023
    nspace <- ceiling(wells2[1,3] / spacing) * 2
    nspace1 <- nspace
    spacing <- wells2[1,3] / nspace
  }
  
  ### psi values at stagnation points ------------------------------------------
  ### previously used the revised psi matrix to determine what the additional contour values should be for cap zone boundaries going through the stagnation points
  ### this time: use given psi matrix
  
  temp_psi <- matrix(0.0, nrow(xy_stg), 1)
  
  stgcounter <- 1
  for(i in 1:nrow(xy_stg)){
    if(ceiling((xy_stg[i,1] - minx)/xstep) <= n && ceiling((xy_stg[i,2] - miny)/ystep) <= n){ # add this 3 Jan 2023; ensure stagnation points in bounds # 25 Jan 2023: fixed problems: should use n not N and there was a missing "<= n" expression
      temp_psi[stgcounter,1] <- psi[ceiling((xy_stg[i,1] - minx)/xstep), ceiling((xy_stg[i,2] - miny)/ystep)] # problem here: if stagnation point is outside outer arc for wedge aquifer, then subscript out of bounds
      stgcounter <- stgcounter + 1
    }
  }
  
  psiD_stg <- temp_psi # try this to remove some of the extra lines
  
  ### Calculate combinations of sums of pairs of well pumping rates
  
  shift <- wells2[,3]
  
  if(nrow(wells2) > 1){
    for(wi1 in 1:(nrow(wells2) - 1)){
      for(wi2 in (wi1 + 1):nrow(wells2)){
        shift <- c(shift, abs(wells2[wi1,3]) + abs(wells2[wi2,3])) # 24 Apr 2023
        shift <- unique(shift)
      }
    }
    
    # 24 Apr 2023 - adjust shift to remove duplicates  
    shift <- sort(shift)
    shiftindx <- 1
    indx2 <- 1
    for(wi1 in 2:length(shift)){
      if(abs(shift[wi1] - shift[indx2]) > 1e-10){ # don't accept if difference is very small and the numbers are essentially the same
        shiftindx <- c(shiftindx, wi1)
        indx2 <- wi1
      }
    }
    
    shift <- shift[shiftindx]
  }
  
  ### Create a list of stagnation points psi values and adjusted stagnation points' psi values 
  ### by adding (to positive psi values) or subtracting (from negative psi values) the combinations of sums of well pairs
  # try the following for case 1 option 2:
  
  if(length(wells4) > 1){
    psiD_stg4 <- psiD_stg
    
    for(wi1 in 1:length(shift)){
      psiD_stg4 <- rbind(psiD_stg4, t(t(psiD_stg4[which(psiD_stg4 > 0)] + shift[wi1]))) # try to fix finding the last segment - case 1 opt 2 - 24 Apr 2023
      psiD_stg4 <- rbind(psiD_stg4, t(t(psiD_stg4[which(psiD_stg4 < 0)] - shift[wi1]))) # try to fix finding the last segment - case 1 opt 2 - 24 Apr 2023
      psiD_stg4 <- rbind(psiD_stg4, t(t(psiD_stg4[which(psiD_stg4 < 0)] + shift[wi1]))) # try to fix finding the last segment - case 1 opt 2 - 25 Apr 2023
      ### psiD_stg4 <- rbind(psiD_stg4, t(t(psiD_stg4 + shift[wi1]))) # 24 Apr 2023 # does not work for case 1 option 1 but works for case 1 option 2
    }
  
    psiD_stg4 <- unique(psiD_stg4)
    psiD_stg4 <- sort(psiD_stg4)
  }else{
    ### revise list based on checklevel (later) - works so far for all cases except case 1 option 2
    psiD_stg5 <- psiD_stg
    for(wi1 in 1:length(psiD_stg)){
      psiD_stg5 <- c(psiD_stg5, unique(c(psiD_stg[wi1] + c(0, shift), psiD_stg[wi1] - shift))) # 15 May 2023
    }

    psiD_stg5 <- unique(psiD_stg5)
    psiD_stg5 <- sort(psiD_stg5)

    psiD_stg4 <- psiD_stg5
  }
  
  ### ----------------------------------------------------------------------
  
  # Method:
  # For j = seq(N,1, by=-1)
  #   Draw stream function contours in region Ij for capphi_j, partitioning based on whether psi(z) <= 0 (do nothing) or psi(z) > 0 (add Qj to psi(z))
  # End
  
  conLines <- list()
  
  if(nrow(wells2) > 1){
  
    for(i in seq(N,1,by=-1)){
      if(i > 1){
        startIndx <- floor((ybreaks[(i-1)] - miny)/ystep)
        stopIndx <- floor((ybreaks[i] - miny)/ystep)
      }else{
        startIndx <- 1
        stopIndx <- floor((min(ybreaks) - miny)/ystep)
      }
      
      # cycle all psi values in the section, modify if psi[x,y] > 0

      if(i < N){
        for(x in 1:n){
          for(y in startIndx:stopIndx){
            if(is.finite(psi[x,y])){ # helpful for case 2 option 3
              if(psi[x,y] > 0){
                psi[x,y] <- psi[x,y] + abs(wells2[i,3]) # try absolute value; something goes wrong when the negative (injection well) is processed - 24 Apr 2023
              }
            }
          }
        }
      }
      
      tempz <- psi[,startIndx:stopIndx]
      maxz <- max(maxz, tempz[is.finite(tempz)]) # deal with possible NaN or Inf values
      
      ### try to generalize:
      if(nrow(wells2) > 1){ # add this 3 Jan 2023
        nspace <- nspace1 * abs(wells2[i,3])/wells2[which(wells2[,2] == max(wells2[,2])),3] # try absolute value; something goes wrong when the negative (injection well) is processed - 24 Apr 2023
        nspace <- nspace[1] # in case there are ties for the maximum y value # 27 Jan 2023
      }else{
        nspace <- nspace1 # add this 3 Jan 2023
      }
      
      spacing <- abs(wells2[i,3]) / nspace # try absolute value; something goes wrong when the negative (injection well) is processed - 24 Apr 2023
      
      ### Loop over x values and contour psi in segment blocks
      
      for(ii in 1:(N + 2 - i)){
        xbreaks2 <- xbreaks[i:N]
        xbreaks2 <- sort(xbreaks2,decreasing=FALSE)
        tempLines <- list()
        
        if(ii == 1){
          if(i == N){
            psi3 <- contourLines(zD_xy[1:n,1], zD_xy[floor((ybreaks[N] - miny)/ystep):n,2], psi[,floor((ybreaks[N] - miny)/ystep):n], levels=seq(minz,maxz,by=spacing)) # try accounting for xbreaks here
            tempLines <- contourLines(zD_xy[1:n,1], zD_xy[floor((ybreaks[N] - miny)/ystep):n,2], psi[,floor((ybreaks[N] - miny)/ystep):n], levels=psiD_stg4)
            if(length(tempLines) > 0){
              conLLength <- length(conLines)
              conLines[(conLLength + 1):(conLLength + length(tempLines))] <- tempLines
            }
          }
          
          if((floor(((xbreaks2[ii] - minx)/xstep) - 2) - 1) > 2 && (stopIndx - 1) - (startIndx + 1) > 2){
            tempPsi <- contourLines(zD_xy[1:floor(((xbreaks2[ii] - minx)/xstep) - 2),1], zD_xy[(startIndx + 1):(stopIndx - 1),2], psi[1:floor(((xbreaks2[ii] - minx)/xstep) - 2),(startIndx + 1):(stopIndx - 1)], levels=seq(minz,maxz,by=spacing)) # revise to deal with both negative and positive x values; changed the zD_xy index for x
            
            if(length(tempPsi) > 0){# 11 Apr 2023
              psi3[(length(psi3) + 1):(length(psi3) + length(tempPsi))] <- tempPsi
            }
            
            tempLines <- contourLines(zD_xy[1:floor(((xbreaks2[ii] - minx)/xstep) - 2),1], zD_xy[(startIndx + 1):(stopIndx - 1),2], psi[1:floor(((xbreaks2[ii] - minx)/xstep) - 2),(startIndx + 1):(stopIndx - 1)], levels=psiD_stg4) # changed the zD_xy index for x
            
            if(length(tempLines) > 0){
              conLLength <- length(conLines)
              conLines[(conLLength + 1):(conLLength + length(tempLines))] <- tempLines
            }
            
          }
        }else if(ii == (N + 2 - i)){
          
          if(n - floor(((xbreaks2[(ii - 1)] - minx)/xstep) + 2) > 2 && (stopIndx - 1) - (startIndx +1) > 2){
            tempPsi <- contourLines(zD_xy[floor(((xbreaks2[(ii - 1)] - minx)/xstep) + 2):n,1], zD_xy[(startIndx +1):(stopIndx - 1),2], psi[floor(((xbreaks2[(ii - 1)] - minx)/xstep) + 2):n,(startIndx + 1):(stopIndx - 1)], levels=seq(minz,maxz,by=spacing)) # changed the zD_xy index for x
            
            if(length(tempPsi) > 0){ # add this check for mirror case
              psi3[(length(psi3) + 1):(length(psi3) + length(tempPsi))] <- tempPsi
            }
            
            tempLines <- contourLines(zD_xy[floor(((xbreaks2[(ii - 1)] - minx)/xstep) + 2):n,1], zD_xy[(startIndx + 1):(stopIndx - 1),2], psi[floor(((xbreaks2[(ii - 1)] - minx)/xstep) + 2):n,(startIndx + 1):(stopIndx - 1)], levels=psiD_stg4) # changed the zD_xy index for x
              
            if(length(tempLines) > 0){
              conLLength <- length(conLines)
              conLines[(conLLength + 1):(conLLength + length(tempLines))] <- tempLines
            }
            
          }
        }else{
          
          if(floor((xbreaks2[ii] - minx)/xstep) - floor((xbreaks2[(ii - 1)] - minx)/xstep) > 5 && (stopIndx - 1) - (startIndx +1) > 2){ # otherwise, there is a selection of one column or less
            tempPsi <- contourLines(zD_xy[floor(((xbreaks2[(ii - 1)] - minx)/xstep) + 2):floor(((xbreaks2[ii] - minx)/xstep) - 2),1], zD_xy[(startIndx +1):(stopIndx - 1),2], psi[floor(((xbreaks2[(ii - 1)] - minx)/xstep) + 2):floor(((xbreaks2[ii] - minx)/xstep) - 2),(startIndx + 1):(stopIndx - 1)], levels=seq(minz,maxz,by=spacing)) # changed the zD_xy index for x
            
            if(length(tempPsi) > 0){# 11 Apr 2023
              psi3[(length(psi3) + 1):(length(psi3) + length(tempPsi))] <- tempPsi
            }
                
            tempLines <- contourLines(zD_xy[floor(((xbreaks2[(ii - 1)] - minx)/xstep) + 2):floor(((xbreaks2[ii] - minx)/xstep) - 2),1], zD_xy[(startIndx +1):(stopIndx - 1),2], psi[floor(((xbreaks2[(ii - 1)] - minx)/xstep) + 2):floor(((xbreaks2[ii] - minx)/xstep) - 2),(startIndx + 1):(stopIndx - 1)], levels=psiD_stg4) # changed the zD_xy index for x
            
            if(length(tempLines) > 0){
              conLLength <- length(conLines)
              conLines[(conLLength + 1):(conLLength + length(tempLines))] <- tempLines
            }
          }
        }
        
      }
      
    }
  }else{ # 25 Apr 2023
    # Only one well  
    startIndx <- 1
    stopIndx <- floor((min(ybreaks) - miny)/ystep)
    
    xbreaks2 <- xbreaks[1]
    tempLines <- list()
    
    psi3 <- contourLines(zD_xy[1:n,1], zD_xy[floor((ybreaks[1] - miny)/ystep):n,2], psi[,floor((ybreaks[1] - miny)/ystep):n], levels=seq(minz,maxz,by=spacing))
    tempLines <- contourLines(zD_xy[1:n,1], zD_xy[floor((ybreaks[1] - miny)/ystep):n,2], psi[,floor((ybreaks[1] - miny)/ystep):n], levels=psiD_stg4)
    if(length(tempLines) > 0){
      conLLength <- length(conLines)
      conLines[(conLLength + 1):(conLLength + length(tempLines))] <- tempLines
    }
    
    tempPsi <- contourLines(zD_xy[c(1:floor(((xbreaks2[1] - minx)/xstep) - 2), floor(((xbreaks2[1] - minx)/xstep) + 1):n),1], zD_xy[(startIndx + 1):(stopIndx - 1),2], psi[c(1:floor(((xbreaks2[1] - minx)/xstep) - 2), floor(((xbreaks2[1] - minx)/xstep) + 1):n),(startIndx + 1):(stopIndx - 1)], levels=seq(minz,maxz,by=spacing)) # revise to deal with both negative and positive x values; changed the zD_xy index for x
    
    if(length(tempPsi) > 0){# 11 Apr 2023
      psi3[(length(psi3) + 1):(length(psi3) + length(tempPsi))] <- tempPsi
    }
    
    tempLines <- contourLines(zD_xy[c(1:floor(((xbreaks2[1] - minx)/xstep) - 2), floor(((xbreaks2[1] - minx)/xstep) + 1):n),1], zD_xy[(startIndx + 1):(stopIndx - 1),2], psi[c(1:floor(((xbreaks2[1] - minx)/xstep) - 2), floor(((xbreaks2[1] - minx)/xstep) + 1):n),(startIndx + 1):(stopIndx - 1)], levels=psiD_stg4) # changed the zD_xy index for x
    
    if(length(tempLines) > 0){
      conLLength <- length(conLines)
      conLines[(conLLength + 1):(conLLength + length(tempLines))] <- tempLines
    }
    
    
  }
  
  # plot(c(0,1), c(0,0), xlim=c(0,1), ylim=c(0,1))
  # plot(c(0,1), c(0,0), xlim=c(-0.5,1), ylim=c(0,1))
  
  # for(i2 in 1:length(psi3)){ # plot the streamlines
    # lines(psi3[[i2]]$x, psi3[[i2]]$y, col="red")
  # }
  
  ### Find and group the contour segments related to lines through the stagnation points-----------------
  
  if(length(conLines) > 1){ # 11 Apr 2023
    ### at this point, the segments need to be cleaned up (sharp angles removed, lines going through wells removed)
    ### then the adjacent band needs to be addressed
    
    ### identify sharp breaks in a segment and split
    toldistedge <- 0.01
    toldist <- 0.025 # 20 Mar 2023
    conLLength2 <- length(conLines)
    splitflags <- c(0,0) # 10 May 2023 - a pair in a row has been split and should not be re-joined
    removeflag <- matrix(0,1,length(conLines))
    rfcount <- 1
    
    # 9 May 2023
    threshAng <- 0.2
    indxskip <- 3 # 27 Jan 2023
    
    for(i in 1:conLLength2){
      temp_pts <- cbind(conLines[[i]]$x, conLines[[i]]$y)
      angmaxdiff <- 0
      indxAngMaxD <- 0
      splitLine <- FALSE # 9 May 2023
      if(nrow(temp_pts) > 10){ # 9 May 2023
        for(ii in (indxskip + 2):(nrow(temp_pts) - 5)){ # 9 May 2023
          mergeResult <- mergecheck(temp_pts[1:ii,], temp_pts[(ii+1):nrow(temp_pts),], toldist, indxskip, 0.6) # 10 May 2023
          
          if(mergeResult != 2 && mergeResult != 3){ # 12 May 2023 # removing some unnecessarily?
            splitLine <- TRUE
            if(indxAngMaxD[1] == 0){ # 9 May 2023
              indxAngMaxD <- ii
            }else if(indxAngMaxD[1] > 0){
              indxAngMaxD <- c(indxAngMaxD, ii)
            }
          }
        }
        
        if(splitLine == TRUE){ # 9 May 2023
            ### split the segment into two
          if(indxAngMaxD[1] > 0){ # 9 May 2023 -- revised to account for multiple breaks
            
            # ignore a cluster of breaks for now
            
            if(length(indxAngMaxD) == 1){
              list1 <- list(level = conLines[[i]]$level, x = conLines[[i]]$x[1:indxAngMaxD], y = conLines[[i]]$y[1:indxAngMaxD])
              list2 <- list(level = conLines[[i]]$level, x = conLines[[i]]$x[(indxAngMaxD+1):nrow(temp_pts)], y = conLines[[i]]$y[(indxAngMaxD+1):nrow(temp_pts)])
              
              conLLength <- length(conLines)
              
              conLines[[1 + conLLength]] <- list1
              conLines[[2 + conLLength]] <- list2
              
              splitflags <- rbind(splitflags, c(1 + conLLength, 2 + conLLength)) # 10 May 2023
              removeflag[rfcount] <- i
              rfcount <- rfcount + 1
            }else if(length(indxAngMaxD) > 1){ # 9 May 2023
              if(indxAngMaxD[1] > 5){
                list1 <- list(level = conLines[[i]]$level, x = conLines[[i]]$x[1:indxAngMaxD[1]], y = conLines[[i]]$y[1:indxAngMaxD[1]])
                conLLength <- length(conLines)
                conLines[[1 + conLLength]] <- list1
                
                splitflags <- rbind(splitflags, c(conLLength, 1 + conLLength)) # 10 May 2023
              }
              
              for(i2 in 1:(length(indxAngMaxD) - 1)){
                if(indxAngMaxD[i2 + 1] - indxAngMaxD[i2] > 5){
                  list1 <- list(level = conLines[[i]]$level, x = conLines[[i]]$x[(indxAngMaxD[i2]+1):indxAngMaxD[i2 + 1]], y = conLines[[i]]$y[(indxAngMaxD[i2]+1):indxAngMaxD[i2 + 1]])
                  conLLength <- length(conLines)
                  conLines[[1 + conLLength]] <- list1

                  splitflags <- rbind(splitflags, c(conLLength, 1 + conLLength)) # 10 May 2023
                }
              }
              
              list1 <- list(level = conLines[[i]]$level, x = conLines[[i]]$x[(indxAngMaxD[length(indxAngMaxD)] + 1):nrow(temp_pts)], y = conLines[[i]]$y[(indxAngMaxD[length(indxAngMaxD)]+1):nrow(temp_pts)])
              conLLength <- length(conLines)
              conLines[[1 + conLLength]] <- list1
              splitflags <- rbind(splitflags, c(conLLength, 1 + conLLength)) # 10 May 2023
              removeflag[rfcount] <- i
              rfcount <- rfcount + 1
            }
          }
        }
      }
    }
    
    # now there are rfcount - 1 indices in conLines that are flagged for removal
    
    
    ### ---------------------------
    ### remove cross-cutting lines? - Next version
    # conLines <- removeCrossCuts(conLines, removeflag, wells, toldist) # 28 Apr 2023
    
    ### remove lines through wells - Next version
    # removeflag_rev <- getLinesExtWells(conLines, removeflag, wells, toldist) # 3 May 2023
    # rfcount <- rfcount + length(removeflag_rev)
    # removeflag <- c(removeflag, removeflag_rev)
    
    ### ---------------------------
    
    threshAng <- 0.1
    indxskip <- 3 # 27 Jan 2023 - the correct segments may not get merged because points at the end of the segment are too close together and do not give a representative angle
    flags <- matrix(0,1,length(conLines))
    capzones <- list()
    
    for(i in 1:length(psiD_stg)){
      ### for each stagnation point, identify list of possible segments
      nseg <- 0
      seg <- 0
      segchosen <- 0
      
      checklevel <- unique(c(psiD_stg[i] + c(0, shift), psiD_stg[i] - shift))
      
      ### Compile segments in vector 'seg' that have a psi value that is a potential match with the stagnation point value i
      for(ii in 1:length(conLines)){
        
        if(length(which(removeflag == ii)) == 0 && length(which(abs(checklevel - conLines[[ii]]$level) == 0)) > 0){
          nseg <- nseg + 1
          if(nseg == 1){
            seg <- ii
          }else{
            seg <- c(seg, ii)
          }
        }
      }
      
      if(nseg > 0){ ### if some segments found...
        for(ii in 1:nseg){ # for the number of possible segments
          temp_pts <- cbind(conLines[[seg[ii]]]$x, conLines[[seg[ii]]]$y)
          x1 <- temp_pts[1,1]
          y1 <- temp_pts[1,2]
          x2 <- temp_pts[nrow(temp_pts),1]
          y2 <- temp_pts[nrow(temp_pts),2]
          
          #### identify segments with endpoints closest to stagnation point related to psiD_stg[i]
          flagClose <- FALSE
          
          dist1 <- sqrt((xy_stg[i,1] - x1)^2 + (xy_stg[i,2] - y1)^2) # compare endpoint 1 with stagnation point
          if(dist1 < toldist){
            flagClose <- TRUE
          }
          
          dist2 <- sqrt((xy_stg[i,1] - x2)^2 + (xy_stg[i,2] - y2)^2) # compare endpoint 2 with stagnation point
          if(dist2 < toldist){
            flagClose <- TRUE
          }
          
          if(flagClose == TRUE && dist1 < toldist && dist2 >= toldist){
            
            # if closest endpoint goes near a pumping well, ignore
            ### verify that neither endpoint terminates near a well
            dist <- sqrt((wells3[,1] - x1)^2 + (wells3[,2] - y1)^2) # compare tip of first segment with all extraction wells # 11 Apr 2023
            if(min(dist) < toldist){
              flags[1,seg[ii]] <- seg[ii]
              flagClose <- FALSE
            }
            
            dist <- sqrt((wells3[,1] - x2)^2 + (wells3[,2] - y2)^2) # compare tail of first segment with all extraction wells # 11 Apr 2023
            if(min(dist) < toldist){
              flags[1,seg[ii]] <- seg[ii]
              flagClose <- FALSE
            }
          }else if(flagClose == TRUE && dist2 < toldist && dist1 >= toldist){
            # if closest endpoint goes near the well, ignore
            ### verify that neither endpoint terminates near a well
            dist <- sqrt((wells3[,1] - x1)^2 + (wells3[,2] - y1)^2) # compare tip of first segment with all wells # 11 Apr 2023
            if(min(dist) < toldist){
              flags[1,seg[ii]] <- seg[ii]
              flagClose <- FALSE
            }
            
            dist <- sqrt((wells3[,1] - x2)^2 + (wells3[,2] - y2)^2) # compare tail of first segment with all wells # 11 Apr 2023
            if(min(dist) < toldist){
              flags[1,seg[ii]] <- seg[ii]
              flagClose <- FALSE
            }
            
          }else if(dist1 < toldist && dist2 < toldist){
            # if closest endpoint goes near the well, ignore
            ### verify that neither endpoint terminates near a well
            dist <- sqrt((wells3[,1] - x1)^2 + (wells3[,2] - y1)^2) # compare tip of first segment with all wells # 11 Apr 2023
            if(min(dist) < toldist){
              flags[1,seg[ii]] <- seg[ii]
              flagClose <- FALSE
            }
            
            dist <- sqrt((wells3[,1] - x2)^2 + (wells3[,2] - y2)^2) # compare tail of first segment with all wells # 11 Apr 2023
            if(min(dist) < toldist){
              flags[1,seg[ii]] <- seg[ii]
              flagClose <- FALSE
            }
          }
          
          ### Estimate whether the segment is directed toward the well rather than tangentially
          #### Cosine Law
          ### find well closest to the stagnation point
          dist <- sqrt((xy_stg[i,1] - wells3[,1])^2 + (xy_stg[i,2] - wells3[,2])^2) # 11 Apr 2023
          indx <- which(dist == min(dist))
          a1 <- sqrt((x1 - x2)^2 + (y1 - y2)^2)
          b1 <- sqrt((wells3[indx,1] - x2)^2 + (wells3[indx,2] - y2)^2) # 11 Apr 2023
          c1 <- sqrt((x1 - wells3[indx,1])^2 + (y1 - wells3[indx,2])^2) # 11 Apr 2023
          theta <- acos((a1*a1 - b1*b1 - c1*c1) / (-2*b1*c1))
          
          dist1 <- sqrt((x1 - x2)^2 + (y1 - y2)^2) # distance between the two endpoints of the segment
          dist2 <- sqrt(2 * min(dist) * min(dist) *(1 - cos(theta))) # expected approx distance based on a triangle within a circle around the well with a radius of the well-stagnation point distance
          
          if((theta < 0.5 * threshAng && dist1 > 1.2 * dist2) || (theta < 0.5 * threshAng && dist1 < 0.8 * dist2) || (dist1 < 0.0001 && dist2 < 0.0001)){ # 26 April 2023 - # use 20% as a metric and verify that the distances are not both essentially zero (overlapping points); works for case 2 opt3 and case 1 opt 1; not for case 1 opt 2 (crashes)
          #   ### segment directed toward well
             flagClose <- FALSE
          }
          
          ### Conclusion - if flagged, then remove
          if(flags[1,seg[ii]] > 0){
            removeflag[rfcount] <- seg[ii]
            rfcount <- rfcount + 1
          }
          
          if(flagClose == TRUE && length(which(removeflag == seg[ii])) == 0){
            if(segchosen[1] == 0){
              segchosen <- seg[ii]
            }else if(segchosen[1] != seg[ii]){ # avoid duplicates # 8 May 2023
              segchosen <- c(segchosen, seg[ii])
            }
            
          }
        }
        
        #### Continue with segment(s) chosen, find next segment that seems to link at endpoint -----
        seglist <- 0 # track list of segments
        segorder <- 0 # 1 means order OK; 2 means order should be reversed
        mergecodelist <- 0 # track mergetype
        threshAngTrial <- threshAng
        
        if(length(segchosen) == 4){ # 17 May 2023 for case 1 option 2
          temp_pts <- cbind(conLines[[segchosen[1]]]$x, conLines[[segchosen[1]]]$y)
          temp_pts2 <- cbind(conLines[[segchosen[2]]]$x, conLines[[segchosen[2]]]$y)
          temp_pts3 <- cbind(conLines[[segchosen[3]]]$x, conLines[[segchosen[3]]]$y)
          temp_pts4 <- cbind(conLines[[segchosen[4]]]$x, conLines[[segchosen[4]]]$y)
          
          ### Does one segment go near a well? (If so, don't use that one)
          
          mergecode1 <- mergecheck(temp_pts, temp_pts2, toldist, indxskip, threshAngTrial * 2)
          mergecode2 <- mergecheck(temp_pts, temp_pts3, toldist, indxskip, threshAngTrial * 2)
          mergecode3 <- mergecheck(temp_pts, temp_pts4, toldist, indxskip, threshAngTrial * 2)
          mergecode4 <- mergecheck(temp_pts2, temp_pts3, toldist, indxskip, threshAngTrial * 2)
          mergecode5 <- mergecheck(temp_pts2, temp_pts4, toldist, indxskip, threshAngTrial * 2)
          mergecode6 <- mergecheck(temp_pts3, temp_pts4, toldist, indxskip, threshAngTrial * 2)
          
          if(mergecode1 > 0 && mergecode2 == 0 && mergecode3 == 0 && checkSplits(splitflags, segchosen[1], segchosen[2]) == FALSE){
            # keep 1 and 2
            segchosen <- segchosen[1:2]
          }else if(mergecode2 > 0 && mergecode1 == 0 && mergecode3 == 0 && checkSplits(splitflags, segchosen[2], segchosen[3]) == FALSE){
            # keep 2 and 3
            segchosen <- segchosen[2:3]
          }else if(mergecode3 > 0 && mergecode1 == 0 && mergecode2 == 0 && checkSplits(splitflags, segchosen[1], segchosen[3]) == FALSE){
            # keep 1 and 3
            segchosen <- c(segchosen[1], segchosen[3])
          }
          
        }else if(length(segchosen) == 3){ # remove one
          temp_pts <- cbind(conLines[[segchosen[1]]]$x, conLines[[segchosen[1]]]$y)
          temp_pts2 <- cbind(conLines[[segchosen[2]]]$x, conLines[[segchosen[2]]]$y)
          temp_pts3 <- cbind(conLines[[segchosen[3]]]$x, conLines[[segchosen[3]]]$y)
          mergecode1 <- mergecheck(temp_pts, temp_pts2, toldist, indxskip, threshAngTrial * 2)
          mergecode2 <- mergecheck(temp_pts2, temp_pts3, toldist, indxskip, threshAngTrial * 2)
          mergecode3 <- mergecheck(temp_pts, temp_pts3, toldist, indxskip, threshAngTrial * 2)
          
          if(mergecode1 > 0 && mergecode2 == 0 && mergecode3 == 0 && checkSplits(splitflags, segchosen[1], segchosen[2]) == FALSE){
            # keep 1 and 2
            segchosen <- segchosen[1:2]
          }else if(mergecode2 > 0 && mergecode1 == 0 && mergecode3 == 0 && checkSplits(splitflags, segchosen[2], segchosen[3]) == FALSE){
            # keep 2 and 3
            segchosen <- segchosen[2:3]
          }else if(mergecode3 > 0 && mergecode1 == 0 && mergecode2 == 0 && checkSplits(splitflags, segchosen[1], segchosen[3]) == FALSE){
            # keep 1 and 3
            segchosen <- c(segchosen[1], segchosen[3])
          }
          
          if(length(segchosen) == 3){
            temp_pts <- cbind(conLines[[segchosen[1]]]$x, conLines[[segchosen[1]]]$y)
            temp_pts2 <- cbind(conLines[[segchosen[2]]]$x, conLines[[segchosen[2]]]$y)
            temp_pts3 <- cbind(conLines[[segchosen[3]]]$x, conLines[[segchosen[3]]]$y)
            mergecode1 <- mergecheck(temp_pts, temp_pts2, toldist, indxskip, threshAngTrial * 4)
            mergecode2 <- mergecheck(temp_pts2, temp_pts3, toldist, indxskip, threshAngTrial * 4)
            mergecode3 <- mergecheck(temp_pts, temp_pts3, toldist, indxskip, threshAngTrial * 4)
            
            if(mergecode1 > 0 && mergecode2 == 0 && mergecode3 == 0 && checkSplits(splitflags, segchosen[1], segchosen[2]) == FALSE){
              # keep 1 and 2
              segchosen <- segchosen[1:2]
            }else if(mergecode2 > 0 && mergecode1 == 0 && mergecode3 == 0 && checkSplits(splitflags, segchosen[2], segchosen[3]) == FALSE){
              # keep 2 and 3
              segchosen <- segchosen[2:3]
            }else if(mergecode3 > 0 && mergecode1 == 0 && mergecode2 == 0 && checkSplits(splitflags, segchosen[1], segchosen[3]) == FALSE){
              # keep 1 and 3
              segchosen <- c(segchosen[1], segchosen[3])
            }else{
              # simply keep 1 and 2
              segchosen <- segchosen[1:2]
            }
          }
        }
        
        if(length(segchosen) == 1){ ### if segchosen has one element
          segchosen[2] <- segchosen[1]
          
          seglist <- segchosen[1]
          segorder <- 2 # test initial reverse
          
        }else if(segchosen[1] != segchosen[2]){ ### if segchosen has two elements or 3 - ignore the third for now
          temp_pts <- cbind(conLines[[segchosen[1]]]$x, conLines[[segchosen[1]]]$y)
          temp_pts2 <- cbind(conLines[[segchosen[2]]]$x, conLines[[segchosen[2]]]$y)
          
          mergecode <- 0
          
          while(mergecode == 0 && threshAngTrial <= 5 * threshAng){ # 3 May 2023 # adjust based on threshAng?
            mergecode <- mergecheck(temp_pts, temp_pts2, toldist, indxskip, threshAngTrial) # call function to evaluate whether to merge; no segments are yet merged, so no need for the mergecheck2 function
            
            if(mergecode == 0){
              threshAngTrial <- threshAngTrial + threshAng
            }
          }
          
          if(mergecode > 0){
            seglist <- c(segchosen[1], segchosen[2])
            if(mergecode == 1){
              segorder <- c(2, 1)
            }else if(mergecode == 2){
              segorder <- c(2, 2)
            }else if(mergecode == 3){
              segorder <- c(1, 1)
            }else if(mergecode == 4){
              segorder <- c(1, 2)
            }
            
            mergecodelist <- mergecode
          }
        }
        
        ### loop while not finished deal with segchosen[1], then deal with segchosen[2] in the same way
        
        for(segi in 1:2){
          if(segchosen[segi] > 0){
            done <- FALSE
            iii <- 1
            
            toldisttrial <- toldist
            threshAngTrial <- threshAng
            
            while(iii <= nseg && done == FALSE){
              foundseg <- FALSE
              flagexttail <- FALSE  #if flagexttail is FALSE, attach next segment to tip; if TRUE, attach to tail
              check1 <- FALSE
              check2 <- FALSE
              
              if(segi == 1){
                temp_pts <- cbind(conLines[[seglist[1]]]$x, conLines[[seglist[1]]]$y)
                
                if(segorder[1] == 1){
                  flagexttail <- FALSE
                }else{
                  flagexttail <- TRUE
                }
                
              }else if(segi == 2){
                temp_pts <- cbind(conLines[[seglist[length(seglist)]]]$x, conLines[[seglist[length(seglist)]]]$y)
                
                if(segorder[length(seglist)] == 1){
                  flagexttail <- TRUE
                }else{
                  flagexttail <- FALSE
                }
              }
              
              
              if(length(seglist) > 0){
                check1 <- checkSplits(splitflags, seglist[1], seg[iii])
                check2 <- checkSplits(splitflags, seglist[length(seglist)], seg[iii])
              }
              
              if(segchosen[1] != seg[iii] && seg[iii] != segchosen[2] && checkSplits(splitflags, segchosen[1], seg[iii]) == FALSE && checkSplits(splitflags, segchosen[2], seg[iii]) == FALSE && length(which(removeflag == seg[iii])) == 0 && length(which(seglist == seg[iii])) == 0 && check1 == FALSE && check2 == FALSE){
                temp_pts2 <- cbind(conLines[[seg[iii]]]$x, conLines[[seg[iii]]]$y)
                mergecode <- 0
                
                if(nrow(temp_pts) > indxskip && nrow(temp_pts2) > indxskip){
                  mergecode <- mergecheck2(temp_pts,temp_pts2, toldisttrial, indxskip, threshAngTrial, flagexttail) # call function to evaluate whether to merge
                }
                
                if(mergecode > 0){
                  if(mergecodelist[1] == 0){
                    mergecodelist <- mergecode
                  }else if(segi == 1){
                    mergecodelist <- c(mergecode, mergecodelist) # added these cases and switched order here on 3 May 2023; mergecodelist had incorrect order
                  }else if(segi == 2){
                    mergecodelist <- c(mergecodelist, mergecode) # original line for either segi == 1 or segi == 2
                  }
                  
                  if(segi == 1){
                    seglist <- c(seg[iii], seglist)
                    
                    if(segorder[1] == 1 && mergecode == 3){
                      segorder <- c(1, segorder)
                    }else if(segorder[1] == 1 && mergecode == 1){
                      segorder <- c(2, segorder)
                    }else if(segorder[1] == 2 && mergecode == 2){
                      segorder <- c(2, segorder)
                    }else if(segorder[1] == 2 && mergecode == 4){
                      segorder <- c(1, segorder)
                    }else if(segorder[1] == 1 && mergecode == 2){ # add additional cases - filled in 26 Apr 2023
                      segorder <- c(1, segorder)
                    }else if(segorder[1] == 1 && mergecode == 4){
                      segorder <- c(2, segorder)
                    }else if(segorder[1] == 2 && mergecode == 1){
                      segorder <- c(1, segorder)
                    }else if(segorder[1] == 2 && mergecode == 3){
                      segorder <- c(2, segorder)
                    }
                  }else if(segi == 2){
                    seglist <- c(seglist, seg[iii])
                    
                    if(segorder[length(segorder)] == 1 && mergecode == 3){
                      segorder <- c(segorder, 1)
                    }else if(segorder[length(segorder)] == 1 && mergecode == 4){
                      segorder <- c(segorder, 2)
                    }else if(segorder[length(segorder)] == 2 && mergecode == 1){
                      segorder <- c(segorder, 1)
                    }else if(segorder[length(segorder)] == 2 && mergecode == 2){
                      segorder <- c(segorder, 2)
                    }else if(segorder[length(segorder)] == 1 && mergecode == 2){ # add additional cases
                      segorder <- c(segorder, 1)
                    }else if(segorder[length(segorder)] == 1 && mergecode == 1){
                      segorder <- c(segorder, 2)
                    }else if(segorder[length(segorder)] == 2 && mergecode == 3){
                      segorder <- c(segorder, 2)
                    }else if(segorder[length(segorder)] == 2 && mergecode == 4){
                      segorder <- c(segorder, 1)
                    }
                  }
                  
                  # update current segment
                  foundseg <- TRUE  # need to go through the list again for the new segment
                  
                  # print(paste(foundseg, seg[iii], toldisttrial, threshAngTrial))
                  
                  x1 <- temp_pts2[1,1]
                  y1 <- temp_pts2[1,2]
                  x2 <- temp_pts2[nrow(temp_pts2),1]
                  y2 <- temp_pts2[nrow(temp_pts2),2]
                  
                  # rgb1 <- runif(1)
                  # rgb2 <- runif(1)
                  # rgb3 <- runif(1)
                  # points(conLines[[seg[iii]]]$x, conLines[[seg[iii]]]$y,col=rgb(rgb1,rgb2,rgb3))
                  # 
                  # if(toldisttrial > toldist){
                    # print(paste(seg[iii], "special case  ///////////////////////////////////////////"))
                    # points(conLines[[seg[iii]]]$x, conLines[[seg[iii]]]$y,col="red")
                  # }
                   
                  toldisttrial <- toldist # reset the toldisttrial in case it was extended
                  threshAngTrial <- threshAng # reset # 3 May 2023
                }
              } # otherwise skip this iii value
              
              #### stop when you reach an outer boundary or run out of segments to check
              if(foundseg == FALSE){
                iii <- iii + 1
                
                x1 <- temp_pts[1,1]
                y1 <- temp_pts[1,2]
                x2 <- temp_pts[nrow(temp_pts),1]
                y2 <- temp_pts[nrow(temp_pts),2]
              }else{
                iii <- 1
              }
              
              # if at edge, done = true
              if(x1 < (minx + toldistedge) || y1 < (miny + toldistedge) || x1 > (maxx - toldistedge) || y1 > (maxy - toldistedge) || x2 < (minx + toldistedge) || y2 < (miny + toldistedge) || x2 > (maxx - toldistedge) || y2 > (maxy - toldistedge)){
                done <- TRUE
              }
              
              # done if near boundary (above case) or if at an injection well
              if(length(wells4) > 3){ # more than one row
                dist <- c(sqrt((x1 - wells4[,1])^2 + (y1 - wells4[,2])^2), sqrt((x2 - wells4[,1])^2 + (y2 - wells4[,2])^2)) # consider both endpoints of the last segment - 19 Apr 2023
                indx <- which(dist == min(dist))
                
                if(dist[indx] < toldistedge){
                  done <- TRUE # near injection well
                  # print("endpoint near injection well ------------------------------------------------------")
                }
              
              }else if(length(wells4) == 3){
                dist <- c(sqrt((x1 - wells4[1])^2 + (y1 - wells4[2])^2), sqrt((x2 - wells4[1])^2 + (y2 - wells4[2])^2)) # 19 Apr 2023
                indx <- which(dist == min(dist))
                
                if(dist[indx] < 2 * toldistedge){
                  done <- TRUE # near injection well
                  # print("endpoint near injection well ------------------------------------------------------")
                }
              }
              
              if(iii == nseg && done == FALSE && toldisttrial < 4 * toldist){ # we've gone through the list and the segments have reached neither a domain boundary nor an injection well for the current segment - 19 Apr 2023
                toldisttrial <- toldisttrial + toldist # extend the valid match distance to try to find more segments
                iii <- 1 # repeat the search for the current segment endpoint
                # print(paste(x1,y2,x2,y2, "update toldist *************************"))
              }else if(iii == nseg && done == FALSE && toldisttrial >= 5 * toldist && threshAngTrial < 6 * threshAng){
                threshAngTrial <- threshAngTrial + threshAng # 3 May 2023
                iii <- 1 # repeat the search for the current segment endpoint
                # print(paste(x1,y2,x2,y2, "update angle *************************"))
              }
            }
          }
          
        }
        # print(seglist)
        # print(segorder)
        # print(mergecodelist)
        
      }
      
      ### check and fix the order of the segments for continuity -- minimize the gaps between segments
      ### for each pair of segments, there are four possible gap distances between the endpoints; choose the minimum one and set the forward and reverse accordingly
      
      if(seglist[1] != 0){ # 17 Apr 2023
        
        for(ii in 1:(length(seglist) - 1)){
          
          d1 <- sqrt((conLines[[seglist[ii]]]$x[1] - conLines[[seglist[ii+1]]]$x[1])^2 + (conLines[[seglist[ii]]]$y[1] - conLines[[seglist[ii+1]]]$y[1])^2) # tip to tip
          d2 <- sqrt((conLines[[seglist[ii]]]$x[1] - conLines[[seglist[ii+1]]]$x[length(conLines[[seglist[ii+1]]]$x)])^2 + (conLines[[seglist[ii]]]$y[1] - conLines[[seglist[ii+1]]]$y[length(conLines[[seglist[ii+1]]]$y)])^2) # tip to tail
          d3 <- sqrt((conLines[[seglist[ii]]]$x[length(conLines[[seglist[ii]]]$x)] - conLines[[seglist[ii+1]]]$x[1])^2 + (conLines[[seglist[ii]]]$y[length(conLines[[seglist[ii]]]$y)] - conLines[[seglist[ii+1]]]$y[1])^2) # tail to tip
          d4 <- sqrt((conLines[[seglist[ii]]]$x[length(conLines[[seglist[ii]]]$x)] - conLines[[seglist[ii+1]]]$x[length(conLines[[seglist[ii+1]]]$x)])^2 + (conLines[[seglist[ii]]]$y[length(conLines[[seglist[ii]]]$y)] - conLines[[seglist[ii+1]]]$y[length(conLines[[seglist[ii+1]]]$y)])^2) # tail to tail
          
          if(d1 < d2 && d1 < d3 && d1 < d4){
            segorder[ii] <- 2 # ensure that directions are 2 then 1
            segorder[ii + 1] <- 1
          }else if(d2 < d1 && d2 < d3 && d2 < d4){
            segorder[ii] <- 2 # ensure that directions are 2 then 2
            segorder[ii + 1] <- 2
          }else if(d3 < d1 && d3 < d2 && d3 < d4){
            segorder[ii] <- 1# ensure that directions are 1 then 1
            segorder[ii + 1] <- 1
          }else if(d4 < d1 && d4 < d2 && d4 < d3){
            segorder[ii] <- 1# ensure that directions are 1 then 2
            segorder[ii + 1] <- 2
          }
        }
        
        ### compiles the capture zone boundary segments in a reasonable order
        if(length(seglist) > 0){
          czlength <- length(capzones)
          
          for(ii in 1:length(seglist)){
            
            if(ii == 1){
              if(segorder[ii] == 1){
                tempx <- conLines[[seglist[ii]]]$x
                tempy <- conLines[[seglist[ii]]]$y
              }else if(segorder[ii] == 2){ # 17 Apr 2023
                xseq <- seq(length(conLines[[seglist[ii]]]$x),1, by=-1)
                tempx <- conLines[[seglist[ii]]]$x[xseq]
                tempy <- conLines[[seglist[ii]]]$y[xseq]
              }
            }else{
              if(segorder[ii] == 1){
                tempx <- c(tempx, conLines[[seglist[ii]]]$x)
                tempy <- c(tempy, conLines[[seglist[ii]]]$y)
              }else if(segorder[ii] == 2){
                xseq <- seq(length(conLines[[seglist[ii]]]$x),1, by=-1)
                tempx <- c(tempx, conLines[[seglist[ii]]]$x[xseq])
                tempy <- c(tempy, conLines[[seglist[ii]]]$y[xseq])
              }
            }
          }
          
          # print("append capzones segments")
          capzones[[czlength + 1]] <- list(x = tempx, y = tempy)
        }
      }
    }
    
    # for(i3 in 1:length(capzones)){
    #   points(capzones[[i3]]$x, capzones[[i3]]$y)
    # }
    
    returnList <- list()
    returnList[[1]] <- psi3
    returnList[[2]] <- capzones
  }else{ # no stagnation points in bounds
    returnList <- list()
    returnList[[1]] <- psi3
    returnList[[2]] <- list()
  }
  
  return(returnList)
}
