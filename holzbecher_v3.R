### holzbecher_v3.R
### Andrew J. Wiebe, 9 Dec 2022; modifes holzbecher_v2.R, 6 Dec 2022; modifies holzbecher_v1.R, 2 Dec 2022
#
# based on test18.R, 2 Dec 2022
# Uses code from nagheli_case1test1.R
#
# PROBLEMS: may need to remove stagnation points outside of the expected range (x = (-1,1), y = (0,2) for case 2 option 4, etc., for different cases)
#
# Revisions since v1:
#    Replaced x_stg with xy_stg; replaced length(x_stg) with nrow(xy_stg);
#    Fixed code for contourLines so that negative x values would be accommodated (index calculations)
#    in xsequence commands: changed the zD_xy index for x
#    (Need to fix code so that y values can be negative) -- adjusted, not tested
#    Works if the wells are not immediately adjacent (i.e., within 7 or so cells in the matrix)
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

holzbecher_v3 <- function(psi, zD_xy, wells, consts, xy_stg){

  if(length(consts) == 5){ # 12 Dec 2022 - add this check to address both cases 1 and 2
    r <- consts[1]
    # beta <- consts[2]
    K <- consts[3]
    b <- consts[4]
  }else{
    r <- consts[1]
    # beta <- consts[3]
    K <- consts[4]
    b <- consts[5]
  }
  
  # q0 <- consts[5]
  # q0D <- q0 / (K * b)
  # 
  # Qw <- wells[,3]
  # QwD <- Qw / (K * b * r)
  
  wells2 <- wells
  wells2[,3] <- wells2[,3] / (K * b * r) # dimensionless Q values
  
  numContoursStr <- 15
  
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
  
  # print(xy_stg_aqf)
  
  checky <- which(xy_stg_aqf[,2] < miny)

  if(length(checky) > 0){
    xy_stg_aqf <- xy_stg_aqf[-checky,,drop=FALSE]
  }

  # print(xy_stg_aqf)
  
  xy_stg <- xy_stg_aqf
  num_stg_pts <- nrow(xy_stg)
  # print(num_stg_pts)
  
  # plot(boundary[,1], boundary[,2], type = "l", cex = 2, xlab="X", ylab="Y",xlim=c(0,1),ylim=c(0,1))
  # lines(c(0,1),c(0,0), lwd = 3, col="black")
  # contour(zD_xy[1:n,2], zD_xy[1:n,2], phiD3, drawlabels=FALSE, col="blue", xlim=c(0,1),ylim=c(0,1))
  
  # Holzbecher 2018 Algorithm, adapted for direction (branch cuts go down instead of left)
  # Divide domain into blocks based on the well locations
  ybreaks <- sort(wells2[,2], decreasing=FALSE)
  # xbreaks <- sort(wells2[,1], decreasing=FALSE)
  wells2 <- wells2[order(wells2[,2], decreasing=FALSE),]
  xbreaks <- wells2[,1]
  # xbreaks2a <- sort(xbreaks,decreasing=FALSE)
  
  # Set capphi[N+1] <- capphi
  # for region with y < ybreaks[1], do nothing (i.e., no changes to psi)
  
  # Draw stream function contours for capphi for region zy > zN (i.e., unaffected region)
  # change the contour interval in the different sections this time
  
  minz <- min(psi)
  maxz <- max(psi)
  spacing <- (maxz - minz) / numContoursStr
  
  nspace <- ceiling(wells2[which(wells2[,2] == max(wells2[,2])),3] / spacing) * 2
  nspace1 <- nspace
  spacing <- wells2[which(wells2[,2] == max(wells2[,2])),3] / nspace
  
  # par(new=TRUE)
  # contour(zD_xy[1:n,2], zD_xy[floor(ybreaks[N] * n):n,2], psi[,floor(ybreaks[N] * n):n], levels=seq(minz,maxz,by=spacing), drawlabels=FALSE, col="red", xlim=c(0,1),ylim=c(0,1))
  # psi3 <- contourLines(zD_xy[1:n,2], zD_xy[floor(ybreaks[N] * n):n,2], psi[,floor(ybreaks[N] * n):n], levels=seq(minz,maxz,by=spacing)) # fix this y index calc
  # psi3 <- contourLines(zD_xy[1:n,1], zD_xy[floor(ybreaks[N] * n):n,2], psi[,floor(ybreaks[N] * n):n], levels=seq(minz,maxz,by=spacing)) # fix this y index calc
  psi3 <- contourLines(zD_xy[1:n,1], zD_xy[floor((ybreaks[N] - miny)/ystep):n,2], psi[,floor((ybreaks[N] - miny)/ystep):n], levels=seq(minz,maxz,by=spacing)) # try accounting for xbreaks here
  
  # psi3 <- contourLines(zD_xy[1:floor((xbreaks2a[1] - minx)/xstep),1], zD_xy[floor((ybreaks[N] - miny)/ystep):n,2], psi[1:floor((xbreaks2a[1] - minx)/xstep),floor((ybreaks[N] - miny)/ystep):n], levels=seq(minz,maxz,by=spacing)) # not needed
  
  # for(i in 1:length(psi3)){ # plot the streamlines
  #   lines(psi3[[i]]$x, psi3[[i]]$y, col="green")
  # }
  
  # for(i in 1:length(xbreaks2a)){ # don't need this stuff; just make sure block between y=0.5 and y=1  does not get duplicated
  #   psi3 <- contourLines(zD_xy[1:n,1], zD_xy[floor((ybreaks[N] - miny)/ystep):n,2], psi[,floor((ybreaks[N] - miny)/ystep):n], levels=seq(minz,maxz,by=spacing))
  # }
  
  # For j = seq(N,1, by=-1)
  #   Draw stream function contours in region Ij for capphi_j, partitioning based on whether psi(z) <= 0 (do nothing) or psi(z) > 0 (add Qj to psi(z))
  # End
  
  for(i in seq(N,1,by=-1)){
    if(i > 1){
      # startIndx <- floor(ybreaks[(i-1)] * n)
      # stopIndx <- floor(ybreaks[i] * n)
      startIndx <- floor((ybreaks[(i-1)] - miny)/ystep)
      stopIndx <- floor((ybreaks[i] - miny)/ystep)
    }else{
      startIndx <- 1
      # stopIndx <- floor(min(ybreaks) * n)
      stopIndx <- floor((min(ybreaks) - miny)/ystep)
    }
    
    # cycle all psi values in the section, modify if psi[x,y] > 0
    for(x in 1:n){
      for(y in startIndx:stopIndx){
        if(psi[x,y] > 0){
          psi[x,y] <- psi[x,y] + wells2[i,3]
        }
      }
    }
    
    maxz <- max(maxz, psi[,startIndx:stopIndx])
    
    # try to generalize:
    nspace <- nspace1 * wells2[i,3]/wells2[which(wells2[,2] == max(wells2[,2])),3]
    
    spacing <- wells2[i,3] / nspace
    
    for(ii in 1:(N + 2 - i)){
      xbreaks2 <- xbreaks[i:N]
      xbreaks2 <- sort(xbreaks2,decreasing=FALSE)
      
      # print(xbreaks2)
      
      # par(new=TRUE)
      if(ii == 1){
        # contour(zD_xy[1:floor(xbreaks2[ii] * n - 2),2], zD_xy[(startIndx + 1):(stopIndx - 1),2], psi[1:floor(xbreaks2[ii] * n - 2),(startIndx + 1):(stopIndx - 1)], levels=seq(minz,maxz,by=spacing), drawlabels=FALSE, col="red", xlim=c(0,1),ylim=c(0,1))
        # tempPsi <- contourLines(zD_xy[1:floor(xbreaks2[ii] * n - 2),2], zD_xy[(startIndx + 1):(stopIndx - 1),2], psi[1:floor(xbreaks2[ii] * n - 2),(startIndx + 1):(stopIndx - 1)], levels=seq(minz,maxz,by=spacing))
        
        if((floor(((xbreaks2[ii] - minx)/xstep) - 2) - 1) > 2 && (stopIndx - 1) - (startIndx + 1) > 2){
          tempPsi <- contourLines(zD_xy[1:floor(((xbreaks2[ii] - minx)/xstep) - 2),1], zD_xy[(startIndx + 1):(stopIndx - 1),2], psi[1:floor(((xbreaks2[ii] - minx)/xstep) - 2),(startIndx + 1):(stopIndx - 1)], levels=seq(minz,maxz,by=spacing)) # revise to deal with both negative and positive x values; changed the zD_xy index for x
          psi3[(length(psi3) + 1):(length(psi3) + length(tempPsi))] <- tempPsi
        }
      }else if(ii == N + 2 - i){
        # contour(zD_xy[floor(xbreaks2[(ii - 1)] * n + 2):n,2], zD_xy[(startIndx +1):(stopIndx - 1),2], psi[floor(xbreaks2[(ii-1)] * n + 2):n,(startIndx + 1):(stopIndx - 1)], levels=seq(minz,maxz,by=spacing), drawlabels=FALSE, col="red", xlim=c(0,1),ylim=c(0,1))
        # tempPsi <- contourLines(zD_xy[floor(xbreaks2[(ii - 1)] * n + 2):n,2], zD_xy[(startIndx +1):(stopIndx - 1),2], psi[floor(xbreaks2[(ii-1)] * n + 2):n,(startIndx + 1):(stopIndx - 1)], levels=seq(minz,maxz,by=spacing))
        if(n - floor(((xbreaks2[(ii - 1)] - minx)/xstep) + 2) > 2 && (stopIndx - 1) - (startIndx +1) > 2){
          tempPsi <- contourLines(zD_xy[floor(((xbreaks2[(ii - 1)] - minx)/xstep) + 2):n,1], zD_xy[(startIndx +1):(stopIndx - 1),2], psi[floor(((xbreaks2[(ii - 1)] - minx)/xstep) + 2):n,(startIndx + 1):(stopIndx - 1)], levels=seq(minz,maxz,by=spacing)) # changed the zD_xy index for x
          psi3[(length(psi3) + 1):(length(psi3) + length(tempPsi))] <- tempPsi
        }
      }else{
        # contour(zD_xy[floor(xbreaks2[(ii - 1)] * n + 2):floor(xbreaks2[ii] * n - 2),2], zD_xy[(startIndx +1):(stopIndx - 1),2], psi[floor(xbreaks2[(ii - 1)] * n + 2):floor(xbreaks2[ii] * n - 2),(startIndx + 1):(stopIndx - 1)], levels=seq(minz,maxz,by=spacing), drawlabels=FALSE, col="red", xlim=c(0,1),ylim=c(0,1))
        # tempPsi <- contourLines(zD_xy[floor(xbreaks2[(ii - 1)] * n + 2):floor(xbreaks2[ii] * n - 2),2], zD_xy[(startIndx +1):(stopIndx - 1),2], psi[floor(xbreaks2[(ii - 1)] * n + 2):floor(xbreaks2[ii] * n - 2),(startIndx + 1):(stopIndx - 1)], levels=seq(minz,maxz,by=spacing))
        if(floor((xbreaks2[ii] - minx)/xstep) - floor((xbreaks2[(ii - 1)] - minx)/xstep) > 5 && (stopIndx - 1) - (startIndx +1) > 2){ # otherwise, there is a selection of one column or less
          tempPsi <- contourLines(zD_xy[floor(((xbreaks2[(ii - 1)] - minx)/xstep) + 2):floor(((xbreaks2[ii] - minx)/xstep) - 2),1], zD_xy[(startIndx +1):(stopIndx - 1),2], psi[floor(((xbreaks2[(ii - 1)] - minx)/xstep) + 2):floor(((xbreaks2[ii] - minx)/xstep) - 2),(startIndx + 1):(stopIndx - 1)], levels=seq(minz,maxz,by=spacing)) # changed the zD_xy index for x
          psi3[(length(psi3) + 1):(length(psi3) + length(tempPsi))] <- tempPsi
        }
      }
      
    }
  }
  
  # for(i in 1:length(psi3)){ # plot the streamlines
  #   lines(psi3[[i]]$x, psi3[[i]]$y, col="green")
  # }
  
  ### Use the revised psi matrix to determine what the additional contour values should be for cap zone boundaries going through the stagnation points
  
  temp_psi <- matrix(0.0, nrow(xy_stg), 1)
  counter <- 1
  for(i in 1:nrow(xy_stg)){
    # temp_psi[counter,1] <- psi[ceiling(xy_stg[i,1] * n), ceiling(xy_stg[i,2] * n)]
    temp_psi[counter,1] <- psi[ceiling((xy_stg[i,1] - minx)/xstep), ceiling((xy_stg[i,2] - miny)/ystep)]
    # temp_psi[counter,1] <- psi[floor((xy_stg[i,1] - minx)/xstep), floor((xy_stg[i,2] - miny)/ystep)] # try floor here -- leads to zero capzone boundaries
    counter <- counter + 1
  }
  
  psiD_stg <- temp_psi
  
  ### Follow each stagnation points contour until boundary or source/sink -----------------------------------
  
  conLines <- list()
  
  ybreaks <- sort(ybreaks, decreasing = TRUE)
  wells3 <- wells2[order(wells2[,2], decreasing=TRUE),] # order by y values
  xbreaks <- sort(xbreaks, decreasing = FALSE)
  # xy_stg <- cbind(x_stg, y_stg)
  # xy_stg_coord <- ceiling(xy_stg * n)
  # xy_stg_coord <- c(ceiling((xy_stg[i,1] - minx)/xstep), ceiling((xy_stg[i,2] - miny)/ystep))
  
  for(i in 1:nrow(xy_stg)){
    if(xy_stg[i,2] > ybreaks[1]){
      # ysequence <- seq(floor(ybreaks[1] * n + 2),n,by=1) # fix this
      ysequence <- seq(floor(((ybreaks[1] - miny)/ystep) + 2),n,by=1)
    }else if(xy_stg[i,2] > ybreaks[length(ybreaks)]){
      # ystart <- floor(ybreaks[min(which(ybreaks < xy_stg[i,2]))] * n) + 2 # fix this
      # yend <- floor(ybreaks[max(which(ybreaks > xy_stg[i,2]))] * n) - 1 # fix this
      ystart <- floor((ybreaks[min(which(ybreaks < xy_stg[i,2]))] - miny)/ystep) + 2
      yend <- floor((ybreaks[max(which(ybreaks > xy_stg[i,2]))] - miny)/ystep) - 1
      ysequence <- seq(ystart,yend,by=1)
    # }else if(xy_stg[i,2] > 0){
    }else if(xy_stg[i,2] > miny){
      ystart <- 1
      # yend <- (ybreaks[length(ybreaks)] * n) - 1 # fix this
      yend <- floor((ybreaks[length(ybreaks)] - miny)/ystep) - 1
      ysequence <- seq(ystart,yend,by=1)
    }
    
    if(xy_stg[i,1] > xbreaks[length(xbreaks)]){
      # xsequence <- seq(floor(xbreaks[length(xbreaks)] * n + 2),n,by=1) # fix this
      xsequence <- seq(floor((xbreaks[length(xbreaks)] - minx)/xstep + 2),n,by=1)
    }else if(xy_stg[i,1] > xbreaks[1]){
      # xstart <- floor(xbreaks[max(which(xbreaks < xy_stg[i,1]))] * n) + 2 # fix this
      # xend <- floor(xbreaks[min(which(xbreaks > xy_stg[i,1]))] * n) - 1 # fix this
      xstart <- floor((xbreaks[max(which(xbreaks < xy_stg[i,1]))] - minx)/xstep + 2)
      xend <- floor((xbreaks[min(which(xbreaks > xy_stg[i,1]))] - minx)/xstep - 1)
      xsequence <- seq(xstart,xend,by=1)
    # }else if(xy_stg[i,1] > 0){
    }else if(xy_stg[i,1] > minx){
      # xsequence <- seq(1,floor(xbreaks[1] * n - 1),by=1) # fix this
      xsequence <- seq(1,floor((xbreaks[1] - minx)/xstep - 1),by=1)
    }
    
    tempLines <- contourLines(zD_xy[xsequence,1], zD_xy[ysequence,2], psi[xsequence,ysequence], levels=psiD_stg[i]) # changed the zD_xy index for x
    
    if(length(tempLines) > 0){
      conLLength <- length(conLines)
      conLines[(conLLength + 1):(conLLength + length(tempLines))] <- tempLines
    }
  }
  
  ### at this point, the segments need to be cleaned up (sharp angles removed, lines going through wells removed)
  ### then the adjacent band needs to be addressed
  
  ### identify sharp breaks in a segment and split
  conLLength2 <- length(conLines)
  threshSlp <- 2.5
  threshAngPerp <- 0.65 # 0.6
  threshAngPara <- 6.0
  splitflags <- matrix(0,length(conLines)*3,length(conLines)*3) # "split matrix" - a 1 implies that the segments were joined but underwent a split operation; do not re-join them!
  removeflag <- matrix(0,1,length(conLines))
  rfcount <- 1
  for(i in 1:conLLength2){
    temp_pts <- cbind(conLines[[i]]$x, conLines[[i]]$y)
    angmaxdiff <- 0
    indxAngMaxD <- 0
    if(nrow(temp_pts) > 5){
      for(ii in 2:(nrow(temp_pts) - 1)){
        i1 <- ii + 1
        i2 <- ii
        slope1 <- (temp_pts[i1,2] - temp_pts[(i1-1),2])/(temp_pts[i1,1] - temp_pts[(i1-1),1])
        slope2 <- (temp_pts[i2,2] - temp_pts[(i2-1),2])/(temp_pts[i2,1] - temp_pts[(i2-1),1])
        # theta1 <- atan(slope1)
        # theta2 <- atan(slope2)
        theta1 <- atan2((temp_pts[i1,2] - temp_pts[(i1-1),2]), (temp_pts[i1,1] - temp_pts[(i1-1),1]))
        theta2 <- atan2((temp_pts[i2,2] - temp_pts[(i2-1),2]), (temp_pts[i2,1] - temp_pts[(i2-1),1]))
        
        # if(abs(theta1 - theta2) > angmaxdiff){
        if(abs(theta1 - theta2) > angmaxdiff && abs(theta1 - theta2) < threshAngPara){
          angmaxdiff <- abs(theta1 - theta2)
          indxAngMaxD <- ii
        }
        # print(paste(ii, theta1, theta2, slope1, slope2))
        # print(paste(ii, abs(theta1 - theta2), abs(slope1 - slope2)))
      }
      
      
      
      if(angmaxdiff > threshAngPerp){
        ### split the segment into two
        list1 <- list(level = conLines[[i]]$level, x = conLines[[i]]$x[1:indxAngMaxD], y = conLines[[i]]$y[1:indxAngMaxD])
        list2 <- list(level = conLines[[i]]$level, x = conLines[[i]]$x[(indxAngMaxD+1):nrow(temp_pts)], y = conLines[[i]]$y[(indxAngMaxD+1):nrow(temp_pts)])
        conLLength <- length(conLines)

        conLines[[1 + conLLength]] <- list1
        conLines[[2 + conLLength]] <- list2
        splitflags[(1 + conLLength),(2 + conLLength)] <- 1
        splitflags[(2 + conLLength),(1 + conLLength)] <- 1
        removeflag[rfcount] <- i
        rfcount <- rfcount + 1
      }
    }
    
  }
  
  # splitflags
  # splitflags2 <- splitflags
  
  # if(rfcount > 1){
    # conLines <- conLines[- removeflag[1:(rfcount - 1)]]
    # splitflags2 <- splitflags2[-removeflag[1:(rfcount - 1)], -removeflag[1:(rfcount - 1)]]
  # }
  
  # splitflags <- splitflags2
  
  ### Problem: after splitting and removing, the indices of the ones to merge are ???
  ### problem persists...
  ### Instead: Try flagging the list items to be removed and then ignore them based on the flags instead of actually removing them
  
  ### from test12.R
  ### Now merge the segments according to similarity of the tips and tails
  ### First: identify which segments should be merged and use the merge matrix to reconstruct what the mergers should be afterward
  
  thresholdNumPts <- 2
  # threshAng <- 1.5 # radians
  threshAng <- 0.6 # radians
  threshSlp <- 2.5
  toldist <- 0.05
  flags <- matrix(0,1,length(conLines))
  mergeflags <- matrix(0,length(conLines),length(conLines)) # "merge matrix"; 1 means merge; the row gives the starting segment, the column gives the second segment; 
  mergetype <- matrix(0,length(conLines),length(conLines)) # the row gives the starting segment, the column gives the second segment; 1 for tip to tip, 2 for tip to tail, 3 for tail to tip, and 4 for tail to tail
  
  for(i in 1:(length(conLines) - 1)){ ## compare endpoints of each set of points to all other endpoints of sets of points
    temp_pts <- cbind(conLines[[i]]$x, conLines[[i]]$y)
    # points(temp_pts[,1],temp_pts[,2],col=rgb(runif(1),runif(1),runif(1)))
    
    if(nrow(temp_pts) > 1){ #use nrow instead
      
      for(ii in (i+1):length(conLines)){
        # if(splitflags[i,ii] == 0){
        if(length(which(removeflag == i)) == 0 && length(which(removeflag == ii)) == 0){
          temp_pts2 <- cbind(conLines[[ii]]$x, conLines[[ii]]$y)
          
          if(nrow(temp_pts2) > 1){ # use nrow instead
            x1 <- temp_pts[1,1]
            y1 <- temp_pts[1,2]
            x2 <- temp_pts[nrow(temp_pts),1]
            y2 <- temp_pts[nrow(temp_pts),2]
            x3 <- temp_pts2[1,1]
            y3 <- temp_pts2[1,2]
            x4 <- temp_pts2[nrow(temp_pts2),1]
            y4 <- temp_pts2[nrow(temp_pts2),2]
            
            ### verify that neither set of points terminates near a well
            dist <- sqrt((wells[,1] - x1)^2 + (wells[,2] - y1)^2) # compare tip of first segment with all wells
            if(min(dist) < toldist){
              flags[1,i] <- i
            }
            
            dist <- sqrt((wells[,1] - x2)^2 + (wells[,2] - y2)^2) # compare tail of first segment with all wells
            if(min(dist) < toldist){
              flags[1,i] <- i
            }
            
            dist <- sqrt((wells[,1] - x3)^2 + (wells[,2] - y3)^2) # compare tip of second segment with all wells
            if(min(dist) < toldist){
              flags[1,ii] <- ii
            }
            
            dist <- sqrt((wells[,1] - x4)^2 + (wells[,2] - y4)^2) # compare tail of second segment with all wells
            if(min(dist) < toldist){
              flags[1,ii] <- ii
            }
            
            # print(paste("flags", flags[1,i], flags[1,ii]))
            
            if(flags[1,i] == 0 && flags[1,ii] == 0){
              ### compare segment endpoints
              dist <- sqrt((x3 - x1)^2 + (y3 - y1)^2) # tip to tip
              # print(dist)
              
              if(dist < toldist){ # merge tip to tip
                # ### check if slope is continuous
                i1 <- 2
                i2 <- 2
                slope1 <- (temp_pts[i1,2] - temp_pts[(i1-1),2])/(temp_pts[i1,1] - temp_pts[(i1-1),1])
                slope2 <- (temp_pts2[i2,2] - temp_pts2[(i2-1),2])/(temp_pts2[i2,1] - temp_pts2[(i2-1),1])
                # theta1 <- atan(slope1)
                # theta2 <- atan(slope2)
                theta1 <- atan2((temp_pts[i1,2] - temp_pts[(i1-1),2]),(temp_pts[i1,1] - temp_pts[(i1-1),1]))
                theta2 <- atan2((temp_pts2[i2,2] - temp_pts2[(i2-1),2]),(temp_pts2[i2,1] - temp_pts2[(i2-1),1]))
                
                # print(abs(theta1 - theta2))
                
                # if(abs(slope1 - slope2) < threshSlp && abs(theta1 - theta2) < threshAng){
                if(abs(theta1 - theta2) < threshAng || abs(theta1 - theta2) > threshAngPara){
                  mergeflags[i,ii] <- 1
                  mergetype[i,ii] <- 1
                }
                
              }else{
                dist <- sqrt((x4 - x1)^2 + (y4 - y1)^2) # tip to tail
                # print(dist)
                
                if(dist < toldist){ # merge tip to tail
                  ### check if slope is continuous
                  i1 <- 2
                  i2 <- nrow(temp_pts2)
                  slope1 <- (temp_pts[i1,2] - temp_pts[(i1-1),2])/(temp_pts[i1,1] - temp_pts[(i1-1),1])
                  slope2 <- (temp_pts2[i2,2] - temp_pts2[(i2-1),2])/(temp_pts2[i2,1] - temp_pts2[(i2-1),1])
                  # theta1 <- atan(slope1)
                  # theta2 <- atan(slope2)
                  theta1 <- atan2((temp_pts[i1,2] - temp_pts[(i1-1),2]),(temp_pts[i1,1] - temp_pts[(i1-1),1]))
                  theta2 <- atan2((temp_pts2[i2,2] - temp_pts2[(i2-1),2]),(temp_pts2[i2,1] - temp_pts2[(i2-1),1]))
                  
                  # print(abs(theta1 - theta2))
                  
                  # if(abs(slope1 - slope2) < threshSlp && abs(theta1 - theta2) < threshAng){
                  if(abs(theta1 - theta2) < threshAng || abs(theta1 - theta2) > threshAngPara){
                    mergeflags[i,ii] <- 1
                    mergetype[i,ii] <- 2
                  }
                }else{
                  dist <- sqrt((x3 - x2)^2 + (y3 - y2)^2) # tail to tip
                  # print(dist)
                  
                  if(dist < toldist){ # merge tip to tail
                    ### check if slope is continuous
                    i1 <- nrow(temp_pts)
                    i2 <- 2
                    slope1 <- (temp_pts[i1,2] - temp_pts[(i1-1),2])/(temp_pts[i1,1] - temp_pts[(i1-1),1])
                    slope2 <- (temp_pts2[i2,2] - temp_pts2[(i2-1),2])/(temp_pts2[i2,1] - temp_pts2[(i2-1),1])
                    # theta1 <- atan(slope1)
                    # theta2 <- atan(slope2)
                    theta1 <- atan2((temp_pts[i1,2] - temp_pts[(i1-1),2]), (temp_pts[i1,1] - temp_pts[(i1-1),1]))
                    theta2 <- atan2((temp_pts2[i2,2] - temp_pts2[(i2-1),2]), (temp_pts2[i2,1] - temp_pts2[(i2-1),1]))
                    
                    # print(abs(theta1 - theta2))
                    
                    # if(abs(slope1 - slope2) < threshSlp && abs(theta1 - theta2) < threshAng){
                    if(abs(theta1 - theta2) < threshAng || abs(theta1 - theta2) > threshAngPara){
                      mergeflags[i,ii] <- 1
                      mergetype[i,ii] <- 3
                    }
                  }else{
                    dist <- sqrt((x4 - x2)^2 + (y4 - y2)^2) # tail to tail
                    # print(dist)
                    
                    if(dist < toldist){ # merge tip to tail
                      ### check if slope is continuous
                      i1 <- nrow(temp_pts)
                      i2 <- nrow(temp_pts2)
                      slope1 <- (temp_pts[i1,2] - temp_pts[(i1-1),2])/(temp_pts[i1,1] - temp_pts[(i1-1),1])
                      slope2 <- (temp_pts2[i2,2] - temp_pts2[(i2-1),2])/(temp_pts2[i2,1] - temp_pts2[(i2-1),1])
                      # theta1 <- atan(slope1)
                      # theta2 <- atan(slope2)
                      theta1 <- atan2((temp_pts[i1,2] - temp_pts[(i1-1),2]), (temp_pts[i1,1] - temp_pts[(i1-1),1]))
                      theta2 <- atan2((temp_pts2[i2,2] - temp_pts2[(i2-1),2]), (temp_pts2[i2,1] - temp_pts2[(i2-1),1]))
                      
                      # print(abs(theta1 - theta2))
                      
                      # if(abs(slope1 - slope2) < threshSlp && abs(theta1 - theta2) < threshAng){
                      if(abs(theta1 - theta2) < threshAng || abs(theta1 - theta2) > threshAngPara){
                        mergeflags[i,ii] <- 1
                        mergetype[i,ii] <- 4
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  # print(mergeflags)
  
  ### from test12.R
  ### Now: find the mergers in mergeflags and connect the sets of points based on the cases specified in mergetype
  ### make this a function later
  con_xy <- list()
  counter_xy <- 1
  mergekey <- matrix(0,sum(mergeflags),7) # format: con_xy element, conLines element A, order of tip (1) and tail (2), conLines element B, order of tip (1) and tail (2)
  mergekeyIndx <- 1
  removeflag2 <- matrix(0,1,length(conLines))
  rfcount2 <- 1
  
  for(i in 1:(nrow(mergeflags) - 1)){
    for(ii in (i + 1):ncol(mergeflags)){
      # if(mergeflags[i,ii] == 1){
      if(mergeflags[i,ii] == 1 &&  length(which(removeflag == i)) == 0 && length(which(removeflag == ii)) == 0){
        ### determine case
        if(length(which(mergekey[,2] == i)) == 0 && length(which(mergekey[,5] == i)) == 0 && length(which(mergekey[,2] == ii)) == 0 && length(which(mergekey[,5] == ii)) == 0){
          kcase <- 1 # neither segment is in con_xy yet
        }else if(length(which(mergekey[,2] == i)) == 1 && length(which(mergekey[,5] == i)) == 0 && length(which(mergekey[,2] == ii)) == 0 && length(which(mergekey[,5] == ii)) == 0){
          kcase <- 2 # segment i is in list and is found in first slot of mergekey
        }else if(length(which(mergekey[,2] == i)) == 0 && length(which(mergekey[,5] == i)) == 1 && length(which(mergekey[,2] == ii)) == 0 && length(which(mergekey[,5] == ii)) == 0){
          kcase <- 3 # segment i is in list and is found in second slot of mergekey
        }else if(length(which(mergekey[,2] == i)) == 0 && length(which(mergekey[,5] == i)) == 0 && length(which(mergekey[,2] == ii)) == 1 && length(which(mergekey[,5] == ii)) == 0){
          kcase <- 4 # segment ii is in list and is found in first slot of mergekey
        }else if(length(which(mergekey[,2] == i)) == 0 && length(which(mergekey[,5] == i)) == 0 && length(which(mergekey[,2] == ii)) == 0 && length(which(mergekey[,5] == ii)) == 1){
          kcase <- 5 # segment ii is in list and is found in second slot of mergekey
        }else if(length(which(mergekey[,2] == i)) == 1 && length(which(mergekey[,5] == i)) == 0 && length(which(mergekey[,2] == ii)) == 0 && length(which(mergekey[,5] == ii)) == 1){
          kcase <- 6 # both segments are in the list; segment i is in first slot and segment ii is in the second slot
        }else if(length(which(mergekey[,2] == i)) == 0 && length(which(mergekey[,5] == i)) == 1 && length(which(mergekey[,2] == ii)) == 1 && length(which(mergekey[,5] == ii)) == 0){
          kcase <- 7 # both segments are in the list; segment ii is in first slot and segment i is in the second slot
        }else if(length(which(mergekey[,2] == i)) == 1 && length(which(mergekey[,5] == i)) == 0 && length(which(mergekey[,2] == ii)) == 1 && length(which(mergekey[,5] == ii)) == 0){
          kcase <- 8 # both segments are in the list and both are in the first slot
        }else if(length(which(mergekey[,2] == i)) == 0 && length(which(mergekey[,5] == i)) == 1 && length(which(mergekey[,2] == ii)) == 0 && length(which(mergekey[,5] == ii)) == 1){
          kcase <- 9 # both segments are in the list and both are in the second slot
        }
        
        if(mergetype[i,ii] == 1){ # tip to tip -----------------------------------
          if(kcase == 1){
            temp_pts <- cbind(conLines[[i]]$x, conLines[[i]]$y)
            temp_pts2 <- cbind(conLines[[ii]]$x, conLines[[ii]]$y)
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts) + nrow(temp_pts2),2)
            con_xy[[counter_xy]][,1] <- c(conLines[[i]]$x[seq(nrow(temp_pts),1,by=-1)], conLines[[ii]]$x)
            con_xy[[counter_xy]][,2] <- c(conLines[[i]]$y[seq(nrow(temp_pts),1,by=-1)], conLines[[ii]]$y)
            mergekey[mergekeyIndx,] <- c(counter_xy, i, 2, 1, ii, 1, 2)
            mergekeyIndx <- mergekeyIndx + 1
            counter_xy <- counter_xy + 1
          }else if(kcase == 2){
            # need to merge with existing merger and then flag the existing one for deletion; overwrite the entry in mergekey
            mkIndx <- which(mergekey[,2] == i) # row in mergekey
            # reverse direction of pts in ii and append to start of row
            temp_pts <- cbind(conLines[[ii]]$x, conLines[[ii]]$y)
            temp_pts3 <- temp_pts[seq(nrow(temp_pts),1,by=-1),]
            temp_pts2 <- con_xy[[mergekey[mkIndx,1]]]
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts3) + nrow(temp_pts2),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts3[,1], temp_pts2[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts3[,2], temp_pts2[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx,] <- c(counter_xy, ii, 2, 1, mergekey[mkIndx, 5:7]) # overwrite the row
            counter_xy <- counter_xy + 1
            
          }else if(kcase == 3){
            mkIndx <- which(mergekey[,5] == i) # row in mergekey
            # append segment ii onto the existing
            temp_pts <- con_xy[[mergekey[mkIndx,1]]]
            temp_pts2 <- cbind(conLines[[ii]]$x, conLines[[ii]]$y)
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts) + nrow(temp_pts2),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts[,1], temp_pts2[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts[,2], temp_pts2[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx,] <- c(counter_xy, mergekey[mkIndx, 2:4], ii, 1, 2) # overwrite the row
            counter_xy <- counter_xy + 1
          }else if(kcase == 4){
            mkIndx <- which(mergekey[,2] == ii) # row in mergekey
            # reverse order of segment i points and append in front of existing
            temp_pts <- cbind(conLines[[i]]$x, conLines[[i]]$y)
            temp_pts3 <- temp_pts[seq(nrow(temp_pts),1,by=-1),]
            temp_pts2 <- con_xy[[mergekey[mkIndx,1]]]
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts3) + nrow(temp_pts2),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts3[,1], temp_pts2[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts3[,2], temp_pts2[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx,] <- c(counter_xy, i, 2, 1, mergekey[mkIndx, 5:7]) # overwrite the row
            counter_xy <- counter_xy + 1
          }else if(kcase == 5){
            mkIndx <- which(mergekey[,5] == ii)
            # append segment i at end of existing
            temp_pts <- con_xy[[mergekey[mkIndx,1]]]
            temp_pts2 <- cbind(conLines[[i]]$x, conLines[[i]]$y)
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts) + nrow(temp_pts2),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts[,1], temp_pts2[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts[,2], temp_pts2[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx,] <- c(counter_xy, mergekey[mkIndx, 2:4], i, 1, 2) # overwrite the row
            counter_xy <- counter_xy + 1
          }else if(kcase == 6){
            mkIndx1 <- which(mergekey[,2] == i)
            mkIndx2 <- which(mergekey[,5] == ii)
            # need to join two existing and then flag both rows for removal
            # append first to second
            temp_pts <- con_xy[[mergekey[mkIndx2,1]]]
            temp_pts2 <- con_xy[[mergekey[mkIndx1,1]]]
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts) + nrow(temp_pts2),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts[,1], temp_pts2[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts[,2], temp_pts2[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx1,1]
            rfcount2 <- rfcount2 + 1
            removeflag2[rfcount2] <- mergekey[mkIndx2,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx1,] <- c(counter_xy, mergekey[mkIndx2, 2:4], mergekey[mkIndx1,5:7]) # overwrite the row
            mergekey[mkIndx2,] <- c(0,0,0,0,0,0,0) # overwrite the row
            counter_xy <- counter_xy + 1
          }else if(kcase == 7){
            mkIndx1 <- which(mergekey[,2] == ii)
            mkIndx2 <- which(mergekey[,5] == i)
            # need to join two existing and then flag both rows for removal
            # append second to first
            temp_pts <- con_xy[[mergekey[mkIndx2,1]]]
            temp_pts2 <- con_xy[[mergekey[mkIndx1,1]]]
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts) + nrow(temp_pts2),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts[,1], temp_pts2[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts[,2], temp_pts2[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx1,1]
            rfcount2 <- rfcount2 + 1
            removeflag2[rfcount2] <- mergekey[mkIndx2,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx2,] <- c(counter_xy, mergekey[mkIndx2, 2:4], mergekey[mkIndx1,5:7]) # overwrite the row
            mergekey[mkIndx1,] <- c(0,0,0,0,0,0,0) # overwrite the row
            counter_xy <- counter_xy + 1
          }else if(kcase == 8){
            mkIndx1 <- which(mergekey[,2] == ii)
            mkIndx2 <- which(mergekey[,2] == i)
            # need to join two existing and then flag both rows for removal
            # reverse and append second to first
            temp_pts <- con_xy[[mergekey[mkIndx1,1]]]
            temp_pts3 <- temp_pts[seq(nrow(temp_pts),1,by=-1),]
            temp_pts2 <- con_xy[[mergekey[mkIndx2,1]]]
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts3) + nrow(temp_pts2),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts3[,1], temp_pts2[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts3[,2], temp_pts2[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx1,1]
            rfcount2 <- rfcount2 + 1
            removeflag2[rfcount2] <- mergekey[mkIndx2,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx2,] <- c(counter_xy, mergekey[mkIndx1, seq(7,5,by=-1)], mergekey[mkIndx2,5:7]) # overwrite the row
            mergekey[mkIndx1,] <- c(0,0,0,0,0,0,0) # overwrite the row
            counter_xy <- counter_xy + 1
          }else if(kcase == 9){
            mkIndx1 <- which(mergekey[,5] == ii)
            mkIndx2 <- which(mergekey[,5] == i)
            # need to join two existing and then flag both rows for removal
            # reverse and append second after first
            temp_pts <- con_xy[[mergekey[mkIndx1,1]]]
            temp_pts3 <- temp_pts[seq(nrow(temp_pts),1,by=-1),]
            temp_pts2 <- con_xy[[mergekey[mkIndx2,1]]]
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts2) + nrow(temp_pts3),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts2[,1], temp_pts3[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts2[,2], temp_pts3[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx1,1]
            rfcount2 <- rfcount2 + 1
            removeflag2[rfcount2] <- mergekey[mkIndx2,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx2,] <- c(counter_xy, mergekey[mkIndx2, 2:4], mergekey[mkIndx1,seq(4,2,by=-1)]) # overwrite the row
            mergekey[mkIndx1,] <- c(0,0,0,0,0,0,0) # overwrite the row
            counter_xy <- counter_xy + 1
          }
        }else if(mergetype[i,ii] == 2){ # tip to tail ----------------------------
          if(kcase == 1){
            temp_pts <- cbind(conLines[[i]]$x, conLines[[i]]$y)
            temp_pts2 <- cbind(conLines[[ii]]$x, conLines[[ii]]$y)
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts) + nrow(temp_pts2),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts2[,1], temp_pts[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts2[,2], temp_pts[,2])
            mergekey[mergekeyIndx,] <- c(counter_xy, ii, 1, 2, i, 1, 2)
            mergekeyIndx <- mergekeyIndx + 1
            counter_xy <- counter_xy + 1
          }else if(kcase == 2){
            # need to merge with existing merger and then flag the existing one for deletion; overwrite the entry in mergekey
            mkIndx <- which(mergekey[,2] == i) # row in mergekey
            # append ii onto the existing list with i at start
            temp_pts <- cbind(conLines[[ii]]$x, conLines[[ii]]$y)
            temp_pts2 <- con_xy[[mergekey[mkIndx,1]]]
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts) + nrow(temp_pts2),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts[,1], temp_pts2[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts[,2], temp_pts2[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx,] <- c(counter_xy, ii, 1, 2, mergekey[mkIndx, 5:7]) # overwrite the row
            counter_xy <- counter_xy + 1
            
          }else if(kcase == 3){
            mkIndx <- which(mergekey[,5] == i) # row in mergekey
            # reverse and append segment ii onto the existing
            temp_pts <- con_xy[[mergekey[mkIndx,1]]]
            temp_pts2 <- cbind(conLines[[ii]]$x, conLines[[ii]]$y)
            temp_pts3 <- temp_pts2[seq(nrow(temp_pts2),1,by=-1),]
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts) + nrow(temp_pts2),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts[,1], temp_pts3[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts[,2], temp_pts3[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx,] <- c(counter_xy, mergekey[mkIndx, 2:4], ii, 2, 1) # overwrite the row
            counter_xy <- counter_xy + 1
          }else if(kcase == 4){
            mkIndx <- which(mergekey[,2] == ii) # row in mergekey
            # reverse order of segment i points and append in front of existing
            temp_pts <- cbind(conLines[[i]]$x, conLines[[i]]$y)
            temp_pts3 <- temp_pts[seq(nrow(temp_pts),1,by=-1),]
            temp_pts2 <- con_xy[[mergekey[mkIndx,1]]]
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts3) + nrow(temp_pts2),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts3[,1], temp_pts2[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts3[,2], temp_pts2[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx,] <- c(counter_xy, i, 2, 1, mergekey[mkIndx, 5:7]) # overwrite the row
            counter_xy <- counter_xy + 1
          }else if(kcase == 5){
            mkIndx <- which(mergekey[,5] == ii)
            # append segment i at end of existing
            temp_pts <- con_xy[[mergekey[mkIndx,1]]]
            temp_pts2 <- cbind(conLines[[i]]$x, conLines[[i]]$y)
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts) + nrow(temp_pts2),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts[,1], temp_pts2[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts[,2], temp_pts2[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx,] <- c(counter_xy, mergekey[mkIndx, 2:4], i, 1, 2) # overwrite the row
            counter_xy <- counter_xy + 1
          }else if(kcase == 6){
            mkIndx1 <- which(mergekey[,2] == i)
            mkIndx2 <- which(mergekey[,5] == ii)
            # need to join two existing and then flag both rows for removal
            # append second to first
            temp_pts <- con_xy[[mergekey[mkIndx2,1]]]
            temp_pts2 <- con_xy[[mergekey[mkIndx1,1]]]
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts) + nrow(temp_pts2),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts[,1], temp_pts2[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts[,2], temp_pts2[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx1,1]
            rfcount2 <- rfcount2 + 1
            removeflag2[rfcount2] <- mergekey[mkIndx2,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx1,] <- c(counter_xy, mergekey[mkIndx2, 2:4], mergekey[mkIndx1,5:7]) # overwrite the row
            mergekey[mkIndx2,] <- c(0,0,0,0,0,0,0) # overwrite the row
            counter_xy <- counter_xy + 1
          }else if(kcase == 7){
            mkIndx1 <- which(mergekey[,2] == ii)
            mkIndx2 <- which(mergekey[,5] == i)
            # need to join two existing and then flag both rows for removal
            # append second to first
            temp_pts <- con_xy[[mergekey[mkIndx2,1]]]
            temp_pts2 <- con_xy[[mergekey[mkIndx1,1]]]
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts) + nrow(temp_pts2),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts[,1], temp_pts2[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts[,2], temp_pts2[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx1,1]
            rfcount2 <- rfcount2 + 1
            removeflag2[rfcount2] <- mergekey[mkIndx2,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx2,] <- c(counter_xy, mergekey[mkIndx2, 2:4], mergekey[mkIndx1,5:7]) # overwrite the row
            mergekey[mkIndx1,] <- c(0,0,0,0,0,0,0) # overwrite the row
            counter_xy <- counter_xy + 1
          }else if(kcase == 8){
            mkIndx1 <- which(mergekey[,2] == ii)
            mkIndx2 <- which(mergekey[,2] == i)
            # need to join two existing and then flag both rows for removal
            # reverse and append second to first
            temp_pts <- con_xy[[mergekey[mkIndx1,1]]]
            temp_pts3 <- temp_pts[seq(nrow(temp_pts),1,by=-1),]
            temp_pts2 <- con_xy[[mergekey[mkIndx2,1]]]
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts3) + nrow(temp_pts2),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts3[,1], temp_pts2[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts3[,2], temp_pts2[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx1,1]
            rfcount2 <- rfcount2 + 1
            removeflag2[rfcount2] <- mergekey[mkIndx2,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx2,] <- c(counter_xy, mergekey[mkIndx1, seq(7,5,by=-1)], mergekey[mkIndx2,5:7]) # overwrite the row
            mergekey[mkIndx1,] <- c(0,0,0,0,0,0,0) # overwrite the row
            counter_xy <- counter_xy + 1
          }else if(kcase == 9){
            mkIndx1 <- which(mergekey[,5] == ii)
            mkIndx2 <- which(mergekey[,5] == i)
            # need to join two existing and then flag both rows for removal
            # reverse and append second after first
            temp_pts <- con_xy[[mergekey[mkIndx1,1]]]
            temp_pts3 <- temp_pts[seq(nrow(temp_pts),1,by=-1),]
            temp_pts2 <- con_xy[[mergekey[mkIndx2,1]]]
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts2) + nrow(temp_pts3),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts2[,1], temp_pts3[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts2[,2], temp_pts3[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx1,1]
            rfcount2 <- rfcount2 + 1
            removeflag2[rfcount2] <- mergekey[mkIndx2,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx2,] <- c(counter_xy, mergekey[mkIndx2, 2:4], mergekey[mkIndx1,seq(4,2,by=-1)]) # overwrite the row
            mergekey[mkIndx1,] <- c(0,0,0,0,0,0,0) # overwrite the row
            counter_xy <- counter_xy + 1
          }
        }else if(mergetype[i,ii] == 3){ # tail to tip ----------------------------
          
          
          if(kcase == 1){
            temp_pts <- cbind(conLines[[i]]$x, conLines[[i]]$y)
            temp_pts2 <- cbind(conLines[[ii]]$x, conLines[[ii]]$y)
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts) + nrow(temp_pts2),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts[,1], temp_pts2[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts[,2], temp_pts2[,2])
            mergekey[mergekeyIndx,] <- c(counter_xy, i, 1, 2, ii, 1, 2)
            mergekeyIndx <- mergekeyIndx + 1
            counter_xy <- counter_xy + 1
          }else if(kcase == 2){
            # need to merge with existing merger and then flag the existing one for deletion; overwrite the entry in mergekey
            mkIndx <- which(mergekey[,2] == i) # row in mergekey
            # append ii onto the existing list that has i at start
            temp_pts <- cbind(conLines[[ii]]$x, conLines[[ii]]$y)
            temp_pts2 <- con_xy[[mergekey[mkIndx,1]]]
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts) + nrow(temp_pts2),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts[,1], temp_pts2[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts[,2], temp_pts2[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx,] <- c(counter_xy, ii, 1, 2, mergekey[mkIndx, 5:7]) # overwrite the row
            counter_xy <- counter_xy + 1
            
          }else if(kcase == 3){
            mkIndx <- which(mergekey[,5] == i) # row in mergekey
            # append segment ii onto the existing
            temp_pts <- con_xy[[mergekey[mkIndx,1]]]
            temp_pts2 <- cbind(conLines[[ii]]$x, conLines[[ii]]$y)
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts) + nrow(temp_pts2),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts[,1], temp_pts2[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts[,2], temp_pts2[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx,] <- c(counter_xy, mergekey[mkIndx, 2:4], ii, 1, 2) # overwrite the row
            counter_xy <- counter_xy + 1
          }else if(kcase == 4){
            mkIndx <- which(mergekey[,2] == ii) # row in mergekey
            # append segment i points onto existing
            temp_pts <- cbind(conLines[[i]]$x, conLines[[i]]$y)
            temp_pts2 <- con_xy[[mergekey[mkIndx,1]]]
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts) + nrow(temp_pts2),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts[,1], temp_pts2[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts[,2], temp_pts2[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx,] <- c(counter_xy, i, 1, 2, mergekey[mkIndx, 5:7]) # overwrite the row
            counter_xy <- counter_xy + 1
          }else if(kcase == 5){
            mkIndx <- which(mergekey[,5] == ii)
            # reverse and append segment i at end of existing
            temp_pts <- con_xy[[mergekey[mkIndx,1]]]
            temp_pts2 <- cbind(conLines[[i]]$x, conLines[[i]]$y)
            temp_pts3 <- temp_pts2[seq(nrow(temp_pts2),1,by=-1),]
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts) + nrow(temp_pts3),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts[,1], temp_pts3[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts[,2], temp_pts3[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx,] <- c(counter_xy, mergekey[mkIndx, 2:4], i, 2, 1) # overwrite the row
            counter_xy <- counter_xy + 1
          }else if(kcase == 6){
            mkIndx1 <- which(mergekey[,2] == i)
            mkIndx2 <- which(mergekey[,5] == ii)
            # need to join two existing and then flag both rows for removal
            # append second to first
            temp_pts <- con_xy[[mergekey[mkIndx2,1]]]
            temp_pts2 <- con_xy[[mergekey[mkIndx1,1]]]
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts) + nrow(temp_pts2),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts[,1], temp_pts2[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts[,2], temp_pts2[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx1,1]
            rfcount2 <- rfcount2 + 1
            removeflag2[rfcount2] <- mergekey[mkIndx2,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx1,] <- c(counter_xy, mergekey[mkIndx2, 2:4], mergekey[mkIndx1,5:7]) # overwrite the row
            mergekey[mkIndx2,] <- c(0,0,0,0,0,0,0) # overwrite the row
            counter_xy <- counter_xy + 1
          }else if(kcase == 7){
            mkIndx1 <- which(mergekey[,2] == ii)
            mkIndx2 <- which(mergekey[,5] == i)
            # need to join two existing and then flag both rows for removal
            # append second after first
            temp_pts <- con_xy[[mergekey[mkIndx2,1]]]
            temp_pts2 <- con_xy[[mergekey[mkIndx1,1]]]
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts) + nrow(temp_pts2),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts[,1], temp_pts2[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts[,2], temp_pts2[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx1,1]
            rfcount2 <- rfcount2 + 1
            removeflag2[rfcount2] <- mergekey[mkIndx2,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx2,] <- c(counter_xy, mergekey[mkIndx2, 2:4], mergekey[mkIndx1,5:7]) # overwrite the row
            mergekey[mkIndx1,] <- c(0,0,0,0,0,0,0) # overwrite the row
            counter_xy <- counter_xy + 1
          }else if(kcase == 8){
            mkIndx1 <- which(mergekey[,2] == ii)
            mkIndx2 <- which(mergekey[,2] == i)
            # need to join two existing and then flag both rows for removal
            # reverse and append second to first
            temp_pts <- con_xy[[mergekey[mkIndx1,1]]]
            temp_pts3 <- temp_pts[seq(nrow(temp_pts),1,by=-1),]
            temp_pts2 <- con_xy[[mergekey[mkIndx2,1]]]
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts3) + nrow(temp_pts2),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts3[,1], temp_pts2[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts3[,2], temp_pts2[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx1,1]
            rfcount2 <- rfcount2 + 1
            removeflag2[rfcount2] <- mergekey[mkIndx2,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx2,] <- c(counter_xy, mergekey[mkIndx1, seq(7,5,by=-1)], mergekey[mkIndx2,5:7]) # overwrite the row
            mergekey[mkIndx1,] <- c(0,0,0,0,0,0,0) # overwrite the row
            counter_xy <- counter_xy + 1
          }else if(kcase == 9){
            mkIndx1 <- which(mergekey[,5] == ii)
            mkIndx2 <- which(mergekey[,5] == i)
            # need to join two existing and then flag both rows for removal
            # reverse and append second after first
            temp_pts <- con_xy[[mergekey[mkIndx1,1]]]
            temp_pts3 <- temp_pts[seq(nrow(temp_pts),1,by=-1),]
            temp_pts2 <- con_xy[[mergekey[mkIndx2,1]]]
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts2) + nrow(temp_pts3),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts2[,1], temp_pts3[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts2[,2], temp_pts3[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx1,1]
            rfcount2 <- rfcount2 + 1
            removeflag2[rfcount2] <- mergekey[mkIndx2,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx2,] <- c(counter_xy, mergekey[mkIndx2, 2:4], mergekey[mkIndx1,seq(4,2,by=-1)]) # overwrite the row
            mergekey[mkIndx1,] <- c(0,0,0,0,0,0,0) # overwrite the row
            counter_xy <- counter_xy + 1
          }
          
          
        }else if(mergetype[i,ii] == 4){ # tail to tail ---------------------------
          
          
          if(kcase == 1){
            temp_pts <- cbind(conLines[[i]]$x, conLines[[i]]$y)
            temp_pts2 <- cbind(conLines[[ii]]$x, conLines[[ii]]$y)
            temp_pts3 <- temp_pts2[seq(nrow(temp_pts2),1,by=-1),]
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts) + nrow(temp_pts2),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts[,1], temp_pts3[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts[,2], temp_pts3[,2])
            mergekey[mergekeyIndx,] <- c(counter_xy, i, 1, 2, ii, 2, 1)
            mergekeyIndx <- mergekeyIndx + 1
            counter_xy <- counter_xy + 1
            
          }else if(kcase == 2){
            # need to merge with existing merger and then flag the existing one for deletion; overwrite the entry in mergekey
            mkIndx <- which(mergekey[,2] == i) # row in mergekey
            # append ii onto the existing list that has i at start
            temp_pts <- cbind(conLines[[ii]]$x, conLines[[ii]]$y)
            temp_pts2 <- con_xy[[mergekey[mkIndx,1]]]
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts) + nrow(temp_pts2),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts[,1], temp_pts2[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts[,2], temp_pts2[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx,] <- c(counter_xy, ii, 1, 2, mergekey[mkIndx, 5:7]) # overwrite the row
            counter_xy <- counter_xy + 1
            
          }else if(kcase == 3){
            mkIndx <- which(mergekey[,5] == i) # row in mergekey
            # reverse and append segment ii after the existing
            temp_pts <- con_xy[[mergekey[mkIndx,1]]]
            temp_pts2 <- cbind(conLines[[ii]]$x, conLines[[ii]]$y)
            temp_pts3 <- temp_pts2[seq(nrow(temp_pts2),1,by=-1),]
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts) + nrow(temp_pts3),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts[,1], temp_pts3[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts[,2], temp_pts3[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx,] <- c(counter_xy, mergekey[mkIndx, 2:4], ii, 2, 1) # overwrite the row
            counter_xy <- counter_xy + 1
          }else if(kcase == 4){
            mkIndx <- which(mergekey[,2] == ii) # row in mergekey
            # append segment i points onto existing
            temp_pts <- cbind(conLines[[i]]$x, conLines[[i]]$y)
            temp_pts2 <- con_xy[[mergekey[mkIndx,1]]]
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts) + nrow(temp_pts2),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts[,1], temp_pts2[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts[,2], temp_pts2[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx,] <- c(counter_xy, i, 1, 2, mergekey[mkIndx, 5:7]) # overwrite the row
            counter_xy <- counter_xy + 1
          }else if(kcase == 5){
            mkIndx <- which(mergekey[,5] == ii)
            # reverse and append segment i at end of existing
            temp_pts <- con_xy[[mergekey[mkIndx,1]]]
            temp_pts2 <- cbind(conLines[[i]]$x, conLines[[i]]$y)
            temp_pts3 <- temp_pts2[seq(nrow(temp_pts2),1,by=-1),]
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts) + nrow(temp_pts3),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts[,1], temp_pts3[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts[,2], temp_pts3[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx,] <- c(counter_xy, mergekey[mkIndx, 2:4], i, 2, 1) # overwrite the row
            counter_xy <- counter_xy + 1
          }else if(kcase == 6){
            mkIndx1 <- which(mergekey[,2] == i)
            mkIndx2 <- which(mergekey[,5] == ii)
            # need to join two existing and then flag both rows for removal
            # append second to first
            temp_pts <- con_xy[[mergekey[mkIndx2,1]]]
            temp_pts2 <- con_xy[[mergekey[mkIndx1,1]]]
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts) + nrow(temp_pts2),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts[,1], temp_pts2[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts[,2], temp_pts2[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx1,1]
            rfcount2 <- rfcount2 + 1
            removeflag2[rfcount2] <- mergekey[mkIndx2,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx1,] <- c(counter_xy, mergekey[mkIndx2, 2:4], mergekey[mkIndx1,5:7]) # overwrite the row
            mergekey[mkIndx2,] <- c(0,0,0,0,0,0,0) # overwrite the row
            counter_xy <- counter_xy + 1
          }else if(kcase == 7){
            mkIndx1 <- which(mergekey[,2] == ii)
            mkIndx2 <- which(mergekey[,5] == i)
            # need to join two existing and then flag both rows for removal
            # append second after first
            temp_pts <- con_xy[[mergekey[mkIndx2,1]]]
            temp_pts2 <- con_xy[[mergekey[mkIndx1,1]]]
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts) + nrow(temp_pts2),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts[,1], temp_pts2[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts[,2], temp_pts2[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx1,1]
            rfcount2 <- rfcount2 + 1
            removeflag2[rfcount2] <- mergekey[mkIndx2,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx2,] <- c(counter_xy, mergekey[mkIndx2, 2:4], mergekey[mkIndx1,5:7]) # overwrite the row
            mergekey[mkIndx1,] <- c(0,0,0,0,0,0,0) # overwrite the row
            counter_xy <- counter_xy + 1
          }else if(kcase == 8){
            mkIndx1 <- which(mergekey[,2] == ii)
            mkIndx2 <- which(mergekey[,2] == i)
            # need to join two existing and then flag both rows for removal
            # reverse and append second to first
            temp_pts <- con_xy[[mergekey[mkIndx1,1]]]
            temp_pts3 <- temp_pts[seq(nrow(temp_pts),1,by=-1),]
            temp_pts2 <- con_xy[[mergekey[mkIndx2,1]]]
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts3) + nrow(temp_pts2),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts3[,1], temp_pts2[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts3[,2], temp_pts2[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx1,1]
            rfcount2 <- rfcount2 + 1
            removeflag2[rfcount2] <- mergekey[mkIndx2,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx2,] <- c(counter_xy, mergekey[mkIndx1, seq(7,5,by=-1)], mergekey[mkIndx2,5:7]) # overwrite the row
            mergekey[mkIndx1,] <- c(0,0,0,0,0,0,0) # overwrite the row
            counter_xy <- counter_xy + 1
          }else if(kcase == 9){
            mkIndx1 <- which(mergekey[,5] == ii)
            mkIndx2 <- which(mergekey[,5] == i)
            # need to join two existing and then flag both rows for removal
            # reverse and append second after first
            temp_pts <- con_xy[[mergekey[mkIndx1,1]]]
            temp_pts3 <- temp_pts[seq(nrow(temp_pts),1,by=-1),]
            temp_pts2 <- con_xy[[mergekey[mkIndx2,1]]]
            con_xy[[counter_xy]] <- matrix(0,nrow(temp_pts2) + nrow(temp_pts3),2)
            con_xy[[counter_xy]][,1] <- c(temp_pts2[,1], temp_pts3[,1])
            con_xy[[counter_xy]][,2] <- c(temp_pts2[,2], temp_pts3[,2])
            removeflag2[rfcount2] <- mergekey[mkIndx1,1]
            rfcount2 <- rfcount2 + 1
            removeflag2[rfcount2] <- mergekey[mkIndx2,1]
            rfcount2 <- rfcount2 + 1
            mergekey[mkIndx2,] <- c(counter_xy, mergekey[mkIndx2, 2:4], mergekey[mkIndx1,seq(4,2,by=-1)]) # overwrite the row
            mergekey[mkIndx1,] <- c(0,0,0,0,0,0,0) # overwrite the row
            counter_xy <- counter_xy + 1
          }
        }
      }
    }
  }
  
  # remove duplicates in removeflag before proceeding
  
  # if(rfcount > 1){
  #   con_xy <- con_xy[- removeflag[1,1:(rfcount - 1)]]
  # }
  
  ### Instead, just ignore the rows of conLines that are listed in the removeflag vector
  # print(length(con_xy))
  ### Add this: Now verify that the chosen combined segments go through a stagnation point
  for(i in 1:length(con_xy)){
    mindist <- 1
    
    # print(i)
    
    for(ii in 1:nrow(con_xy[[i]])){
      dist <- sqrt((con_xy[[i]][ii,1] - xy_stg[,1])^2 + (con_xy[[i]][ii,2] - xy_stg[,2])^2)
      if(min(dist) < mindist){
        mindist <- min(dist)
      }
      
    }
    
    if(mindist < toldist){
      # keep it
    }else{
      removeflag2[rfcount2] <- i
      rfcount2 <- rfcount2 + 1
    }
  }
  
  # for(i in 1:length(con_xy)){
  #   lines(con_xy[[i]][,1], con_xy[[i]][,2], col="black", lwd=2)#, col=colours[i])
  # }
  
  ### Now search for connections to the existing portions of the capture zone boundaries
  toldistedge <- 0.01
  con_xy_ext <- list()
  segments <- matrix(0, length(con_xy) * 2, nrow(wells) * 2) # one row for each end of the segment in con_xy_ext (the first connections with index [1])
  segCounter <- matrix(0, length(con_xy) * 2, 1)
  
  for(i in 1:length(con_xy)){
    if(length(which(removeflag2 == i)) == 0){ # add this check
      temp_pts <- cbind(con_xy[[i]][,1], con_xy[[i]][,2])
      
      # lines(con_xy[[i]][,1], con_xy[[i]][,2], col="black", lwd = 2)
      
      for(ii in 1:2){ # extend each end of the segment
        if(ii == 1){
          x1 <- temp_pts[1,1]
          y1 <- temp_pts[1,2]
          x2 <- temp_pts[2,1]
          y2 <- temp_pts[2,2]
          slope1 <- (y1 - y2)/(x1 - x2)
        }else{
          x1 <- temp_pts[nrow(temp_pts),1]
          y1 <- temp_pts[nrow(temp_pts),2]
          x2 <- temp_pts[(nrow(temp_pts) - 1),1]
          y2 <- temp_pts[(nrow(temp_pts) - 1),2]
          
          # print(paste(x1,y1,x2,y2))
          
          breakIndx <- 0
          
          for(iii in 2:(nrow(temp_pts) - 1)){
            i1 <- iii + 1
            i2 <- iii
            theta1 <- atan2((temp_pts[i1,2] - temp_pts[(i1-1),2]), (temp_pts[i1,1] - temp_pts[(i1-1),1]))
            theta2 <- atan2((temp_pts[i2,2] - temp_pts[(i2-1),2]), (temp_pts[i2,1] - temp_pts[(i2-1),1]))
            
            if(abs(theta1 - theta2) > threshAng && abs(theta1 - theta2) < threshAngPara){
              ### identify breakpoint to use instead of endpoint
              breakIndx <- i2 ## still need to test with Nagheli et al. (2020) case 2 option 4
            }
          }
          
          ### note: need some code above to ensure monotonicity of the curves (would avoid the need for this)!
          ### ensure no drastic breaks in slope
          
          if(breakIndx > 0){ ## still need to test with Nagheli et al. (2020) case 2 option 4
            x1 <- temp_pts[breakIndx,1]
            y1 <- temp_pts[breakIndx,2]
            x2 <- temp_pts[(breakIndx - 1),1]
            y2 <- temp_pts[(breakIndx - 1),2]
          }
          
          slope1 <- (y1-y2)/(x1-x2)
          # print(paste(x1,y1,x2,y2))
        }
        
        # if(x1 < toldistedge || y1 < toldistedge || x1 > (1 - toldistedge) || y1 > (1 - toldistedge)){ # fix
        if(x1 < (minx + toldistedge) || y1 < (miny + toldistedge) || x1 > (maxx - toldistedge) || y1 > (maxy - toldistedge)){ # already at boundary
          if(ii == 1){
            segments[(i * 2 - 1), segCounter[(i * 2 - 1),1]] <- 0
            segCounter[(i * 2 - 1),1] <- 0
          }else{
            segments[(i * 2), segCounter[(i * 2),1]] <- 0
            segCounter[(i * 2),1] <- 0
          }
          ### done
        }else{
          iii <- 1
          done <- FALSE
          while(iii <= nrow(wells) * 2 && done == FALSE) { ## would want to set iii max with a reasonable calculation based on the number of divisions
            
            ## linear extrap across gap
            distx <- x1 - xbreaks
            disty <- y1 - ybreaks
            
            indx <- 1
            for(i3 in 1:length(distx)){
              if(min(abs(distx)) == abs(distx[i3])){
                indx <- i3
              }
            }
            
            # gapx <- floor(xbreaks[indx] * n) # xNext must be > or < this value, depending on direction
            gapx <- floor((xbreaks[indx] - minx)/xstep) # xNext must be > or < this value, depending on direction
            xdir <- 0
            
            # if(x1 - xbreaks[indx] > 0){
            if(x1 < x2){
              xNexttrial <- (gapx - 1) * xstep + minx # going left
              if(xNexttrial == x1){
                xNexttrial <- (gapx - 2) * xstep + minx
              }
              xdir <- -1
            }else{
              xNexttrial <- (gapx + 1) * xstep + minx # going right
              if(xNexttrial == x1){
                xNexttrial <- (gapx + 2) * xstep + minx
              }
              xdir <- 1
            }
            
            ### check if mismatch caused by clustered xbreaks values
            if(xdir > 0 && floor((xNexttrial - minx)/xstep) < gapx){ # mismatch
              # choose a higher xbreaks value
              indx <- indx + 1 ## choose a larger xbreaks value
              if(indx > length(xbreaks)){
                indx <- length(xbreaks) # print("mismatch")
              }
              
              gapx <- floor((xbreaks[indx] - minx)/xstep) # xNext must be > or < this value, depending on direction
              xNexttrial <- (gapx + 1) * xstep + minx # going right
            }else if(xdir < 0 && floor((xNexttrial - minx)/xstep) > gapx){
              # choose a lower xbreaks value
              indx <- indx - 1 ## choose a larger xbreaks value
              if(indx <= 0){
                indx <- 1 # print("mismatch")
              }
              
              gapx <- floor((xbreaks[indx] - minx)/xstep) # xNext must be > or < this value, depending on direction
              xNexttrial <- (gapx - 1) * xstep + minx # going left
            }
            
            indx <- 1
            for(i3 in 1:length(disty)){
              if(min(abs(disty)) == abs(disty[i3])){
                indx <- i3
              }
            }
            
            # gapy <- floor(ybreaks[indx] * n) #yNext must be > or < this value, depending on direction
            gapy <- floor((ybreaks[indx] - miny)/ystep) #yNext must be > or < this value, depending on direction
            ydir <- 0
            
            if(y2 > y1){
              yNexttrial <- (gapy - 1) * ystep + miny # going down
              if(yNexttrial == y1){
                yNexttrial <- (gapy - 2) * ystep + miny
              }
              ydir <- -1
            }else{
              yNexttrial <- (gapy + 1) * ystep + miny # going up
              if(yNexttrial == y1){ # yNexttrial and yNext cannot be equal to y1
                yNexttrial <- (gapy + 2) * ystep + miny
              }
              ydir <- 1
            }
            
            ## if two or more ybreaks are very close together < 7 cells apart, according to the spacing,
            ## and if you are moving toward the confluence of the clustered breaks,
            ## then choose the outer break
            ## i.e., if slope1 and y1 and y2 relationships are out of sync with the ydir as estimated from (y1 - ybreaks[indx] > 0),
            ## then switch breaks
            if(ydir < 0 && floor((yNexttrial - miny)/ystep) > gapy){
            #   # choose a higher ybreak value
              indx <- indx - 1
              if(indx <= 0){
                indx <= 1 # print("mismatch")
              }

              gapy <- floor((ybreaks[indx] - miny)/ystep) #yNext must be > or < this value, depending on direction
              yNexttrial <- (gapy - 1) * ystep + miny # going down

            }else if(ydir > 0 && floor((yNexttrial - miny)/ystep) < gapy){
            # # choose a lower ybreak value
              indx <- indx + 1
              if(indx > length(ybreaks)){
                indx <- length(ybreaks) # print("mismatch")
              }
               
              gapy <- floor((ybreaks[indx] - miny)/ystep) #yNext must be > or < this value, depending on direction
              yNexttrial <- (gapy + 1) * ystep + miny # going up
            }
            
            #### if slope1 is positive and either (xgap is to the right and ygap is above) or (xgap is to the left and ygap is below) -- OK
            checkcombodir <- FALSE
            if(is.na(slope1)){
              checkcombodir <- FALSE
            }else if(is.infinite(slope1)){
              print("error - slope1 = 0/0") # ?????????????????????????? Try to prevent this
            }else if(slope1 > 0 && ((xdir > 0 && ydir > 0) || (xdir < 0 && ydir < 0))){
              checkcombodir <- TRUE
            }else if(slope1 < 0 && ((xdir < 0 && ydir > 0) || (xdir > 0 && ydir < 0))){
              ### if slope1 is negative and either(xgap is to the left and ygap is above) or (xgap is to the right and ygap is below) -- OK
              checkcombodir <- TRUE
            }
            
            checkvertdir <- FALSE
            ### if ygap is above and y1 > y2 upward
            ### if ygap is below and y1 < y2 downward
            # if((ydir > 0 && y1 > y2) ||(ydir < 0 && y1 < y2)){
            if(ydir > 0 || ydir < 0){
              checkvertdir <- TRUE
            }
            
            checkhorizdir <- FALSE
            ### if xgap is to the left and x1 < x2 left
            ### if xgap is to the right and x1 > x2 right
            # if((xdir < 0 && x1 < x2) ||(xdir > 0 && x1 > x2)){
            if(xdir < 0 || xdir > 0){
              checkhorizdir <- TRUE
            }
            
            if(min(abs(distx)) < toldist && min(abs(disty)) < toldist && checkcombodir == TRUE){ # add: and the slope indicates a transition across both gaps
              ### junction between horizontal and vertical bands
              # ** REVISE and use xNexttrial and yNexttrial, then test with Nagheli et al 2020 case 2 option 4
              indx <- 1
              for(i3 in 1:length(distx)){
                if(min(abs(distx)) == abs(distx[i3])){
                  indx <- i3
                }
              }
              
              # gapx <- floor(xbreaks[indx] * n) # xNext must be > or < this value, depending on direction
              gapx <- floor((xbreaks[indx] - minx)/xstep) # xNext must be > or < this value, depending on direction
              xdir <- 0
              
              if(x1 - xbreaks[indx] > 0){
                # xNexttrial <- (gapx - 1) / n # going left
                xNexttrial <- (gapx - 1) * xstep + minx # going left
                xdir <- -1
              }else{
                # xNexttrial <- (gapx + 1) / n # going right
                xNexttrial <- (gapx + 1) * xstep + minx # going right
                xdir <- 1
              }
              
              
              indx <- 1
              for(i3 in 1:length(disty)){
                if(min(abs(disty)) == abs(disty[i3])){
                  indx <- i3
                }
              }
              
              # gapy <- floor(ybreaks[indx] * n) #yNext must be > or < this value, depending on direction
              gapy <- floor((ybreaks[indx] - miny)/ystep) #yNext must be > or < this value, depending on direction
              ydir <- 0
              
              if(y1 - ybreaks[indx] > 0){
                # yNexttrial <- (gapy - 1) / n # going down
                yNexttrial <- (gapy - 1) * ystep + miny # going down
                ydir <- -1
              }else{
                # yNexttrial <- (gapy + 1) / n # going up
                yNexttrial <- (gapy + 1) * ystep + miny # going up
                ydir <- 1
              }
              
              yNext <- slope1 * (xNexttrial - x1) + y1
              xNext <- ((yNext - y1) / slope1) + x1
              
              y_OK <- FALSE
              
              # if(floor(yNext * n) != gapy){
              if(floor((yNext - miny)/ystep) != gapy){
                y_OK <- TRUE
              }else{
                i2 <- 1
                while(y_OK == FALSE && i2 < 10){
                  if(xdir < 0){
                    xNexttrial <- xNexttrial - 1/n # going left
                  }else{
                    xNexttrial <- xNexttrial + 1/n # going right
                  }
                  
                  yNext <- slope1 * (xNexttrial - x1) + y1
                  xNext <- ((yNext - y1) / slope1) + x1
                  
                  # if(floor(yNext * n) != gapy){
                  if(floor((yNext - miny)/ystep) != gapy){
                    y_OK <- TRUE
                  }
                  
                  i2 <- i2 + 1
                }
              }
              
            }else if(min(abs(distx)) < toldist && checkhorizdir == TRUE){ # consider vertical
              ### transition to next vertical band
              xNext <- xNexttrial
              yNext <- slope1 * (xNext - x1) + y1
            }else if(min(abs(disty)) < toldist && checkvertdir == TRUE){
              ### transition to next horizontal band
              yNext <- yNexttrial
              xNext <- ((yNext - y1) / slope1) + x1
            }
            
            # psiVal <- psi[floor(xNext * n), floor(yNext * n)]
            # print(paste(xNext, floor((xNext - minx)/xstep), yNext, floor((yNext - miny)/ystep)))
            psiVal <- psi[floor((xNext - minx)/xstep), floor((yNext - miny)/ystep)]
            
            ### fix the sequences!
            
            if(yNext > ybreaks[1]){
              # ysequence <- seq(floor(yNext * n),n,by=1) # adjust based on points
              ysequence <- seq(floor((yNext - miny)/ystep),n,by=1) # adjust based on points
            }else if(yNext > ybreaks[length(ybreaks)]){
              # ystart <- floor(ybreaks[min(which(ybreaks < yNext))] * n) + 2 ## this version for i == 2 (works for i == 1 as well)
              # yend <- floor(ybreaks[max(which(ybreaks > yNext))] * n) - 1  ## this version for i == 2 (works for i == 1 as well)
              if(ydir > 0){
                # ystart <- floor((ybreaks[min(which(ybreaks < yNext))] - miny)/ystep) + 2 ## this version for i == 2 (works for i == 1 as well)
                ystart <- floor((yNext - miny)/ystep)
                yend <- floor((ybreaks[max(which(ybreaks > yNext))] - miny)/ystep) - 1  ## this version for i == 2 (works for i == 1 as well)
                
              }else{
                ystart <- floor((ybreaks[min(which(ybreaks < yNext))] - miny)/ystep) + 2 ## this version for i == 2 (works for i == 1 as well)
                yend <- floor((yNext - miny)/ystep)
              }
              
              if(yend > ystart){
                ysequence <- seq(ystart,yend,by=1)
              }else{
                ysequence <- yend
              }
            }else if(yNext > miny){ # don't use zero here
              ystart <- 1
              # yend <- floor(yNext * n)
              yend <- floor((yNext - miny)/ystep)
              ysequence <- seq(ystart,yend,by=1)
            }else{ # add this
              ystart <- 1
              yend <- 1
              ysequence <- 1
            }
            
            if(xNext > xbreaks[length(xbreaks)]){
              # xsequence <- seq(floor(xNext * n),n,by=1)
              xsequence <- seq(floor((xNext - miny)/ystep),n,by=1)
            }else if(xNext > xbreaks[1]){
              if(x1 > xNext){ # moving left
                # xstart <- floor(xbreaks[max(which(xbreaks < xNext))] * n) + 2
                # xend <- floor(xNext * n)
                xstart <- floor((xbreaks[max(which(xbreaks < xNext))] - minx)/xstep) + 2
                xend <- floor((xNext - minx)/xstep)
              }else{
                # xstart <- floor(xNext * n)
                # xend <- floor(xbreaks[min(which(xbreaks > xNext))] * n) - 1
                xstart <- floor((xNext - minx)/xstep)
                xend <- floor((xbreaks[min(which(xbreaks > xNext))] - minx)/xstep) - 1
              }
              
              if(xend > xstart){
                xsequence <- seq(xstart,xend,by=1)
              }else{
                xsequence <- xend
              }
            }else if(xNext > minx){ # don't use zero here
              # xsequence <- seq(1,floor(xNext * n),by=1)
              xsequence <- seq(1,floor((xNext - minx)/xstep),by=1)
            }else{
              xsequence <- 1
            }
            
            if(length(xsequence) > 1 && length(ysequence) > 1){ # need a matrix for psi to calc contour lines
              tempLines <- contourLines(zD_xy[xsequence,1], zD_xy[ysequence,2], psi[xsequence,ysequence], levels=psiVal) # changed the zD_xy index for x
              # lines(tempLines[[1]]$x, tempLines[[1]]$y, col="black", lwd = 2, xlab="X", ylab="Y", xlim=c(0,1), ylim=c(0,1))
              
              con_xy_ext_len <- length(con_xy_ext)
              # con_xy_ext[[(con_xy_ext_len + 1):(con_xy_ext_len + length(tempLines))]] <- tempLines # fails for Whitehorse case
              con_xy_ext[[con_xy_ext_len + 1]] <- tempLines[[1]]
              
              if(ii == 1){
                segments[(i * 2 - 1), segCounter[(i * 2 - 1),1]] <- con_xy_ext_len + 1
                segCounter[(i * 2 - 1),1] <- segCounter[(i * 2 - 1),1] + 1
              }else{
                segments[(i * 2), (segCounter[(i * 2),1] + 1)] <- con_xy_ext_len + 1
                segCounter[(i * 2),1] <- segCounter[(i * 2),1] + 1
              }
              
              # recalcuate x1 and y1 (assume list of segments is in order)
              dist1 <- sqrt((tempLines[[1]]$x[1] - x1)^2 + (tempLines[[1]]$y[1] - y1)^2)
              dist2 <- sqrt((tempLines[[length(tempLines)]]$x[length(tempLines[[length(tempLines)]]$x)] - x1)^2 + (tempLines[[length(tempLines)]]$y[length(tempLines[[length(tempLines)]]$y)] - y1)^2)
              
              if(dist2 < dist1){ # check this; incorrect point selected
                x1 <- tempLines[[1]]$x[1]
                y1 <- tempLines[[1]]$y[1]
                x2 <- tempLines[[1]]$x[2]
                y2 <- tempLines[[1]]$y[2]
                slope1 <- (y1 - y2)/(x1 - x2)
              }else{
                # x1 <- tempLines[[length(tempLines)]]$x[length(tempLines[[length(tempLines)]]$x)]
                # y1 <- tempLines[[length(tempLines)]]$y[length(tempLines[[length(tempLines)]]$y)]
                # x2 <- tempLines[[length(tempLines)]]$x[(length(tempLines[[length(tempLines)]]$x) - 1)]
                # y2 <- tempLines[[length(tempLines)]]$y[(length(tempLines[[length(tempLines)]]$y) - 1)]
                x1 <- tempLines[[1]]$x[length(tempLines[[1]]$x)]
                y1 <- tempLines[[1]]$y[length(tempLines[[1]]$y)]
                x2 <- tempLines[[1]]$x[(length(tempLines[[1]]$x) - 1)]
                y2 <- tempLines[[1]]$y[(length(tempLines[[1]]$y) - 1)]
                slope1 <- (y1 - y2)/(x1 - x2)
              }
              
              if(x1 < (minx + toldistedge) || y1 < (miny + toldistedge) || x1 > (maxx - toldistedge) || y1 > (maxy - toldistedge)){ # fixed
                done <- TRUE
              }
            }else{
              #if either xsequence or ysequence is only one element in length
              # print("advance to xNext,yNext")
              # x2 <- x1
              # y2 <- ?
              # xnew <- x1 + xstep * xdir
              # ynew <- slope1 * (xnew - x1) + y1
              # x1 <- xnew
              # y1 <- ynew
              x1 <- xNext
              y1 <- yNext
              slope1 <- (y1 - y2)/(x1 - x2)
              # need a check on the toldistedge?
            }
            
            iii <- iii + 1
            # print(iii)
          }
          
          
          ### identify psi at the new location and contour that value in the new vert or horiz band
          
        }
      }
    }
  }
  
  # lines(c(0,1),c(0,0), lwd = 3, col="black")
  # dev.off()
  # 
  # for(i in 1:length(con_xy_ext)){
  #   points(con_xy_ext[[i]]$x,con_xy_ext[[i]]$y, col="blue")
  # }
  # 
  ### splice together the sets of points making up the capture zone boundary curves
  
  capzones <- list()
  counter <- 1
  
  for(i in 1:length(con_xy)){
    if(length(which(removeflag2 == i)) == 0){ # add this check
      # endpts <- matrix(0, 2 + (segCounter[(i * 2 - 1)] - 1) * 2 + (segCounter[(i * 2)] - 1) * 2, 5) # format: x,y,con_xy(1) or con_xy_ext(2), index (i if con_xy or index in con_xy_ext), order(point at start - 1, point at end of list - 2)
      endpts <- matrix(0, 2 + (segCounter[(i * 2 - 1)]) * 2 + (segCounter[(i * 2)]) * 2, 5) # format: x,y,con_xy(1) or con_xy_ext(2), index (i if con_xy or index in con_xy_ext), order(point at start - 1, point at end of list - 2)
      
      endpts[1,] <- c(con_xy[[i]][1,1], con_xy[[i]][1,2], 1, i, 1)
      endpts[2,] <- c(con_xy[[i]][length(con_xy[[i]][,1]),1], con_xy[[i]][length(con_xy[[i]][,2]),2], 1, i, 2) # added comma before 2 in y coord length selection
      endpts_count <- 3
      
      if(i > 1  && segCounter[(i * 2 - 1)] > 0){ ## add this check
        # for(ii in 1:(segCounter[(i * 2 - 1)] - 1)){ # problem here if only one item in con_xy
        for(ii in 1:(segCounter[(i * 2 - 1)])){
          # endpts[endpts_count,] <- c(con_xy_ext[[segments[(i * 2 - 1),ii]]][[1]]$x[1], con_xy_ext[[segments[(i * 2 - 1),ii]]][[1]]$y[1], 2, segments[(i*2 - 1),ii], 1)
          endpts[endpts_count,] <- c(con_xy_ext[[segments[(i * 2 - 1),ii]]]$x[1], con_xy_ext[[segments[(i * 2 - 1),ii]]]$y[1], 2, segments[(i*2 - 1),ii], 1)
          endpts_count <- endpts_count + 1
          # endpts[endpts_count,] <- c(con_xy_ext[[segments[(i * 2 - 1),ii]]][[1]]$x[length(con_xy_ext[[segments[(i * 2 - 1),ii]]][[1]]$x)], con_xy_ext[[segments[(i * 2 - 1),ii]]][[1]]$y[length(con_xy_ext[[segments[(i * 2 - 1),ii]]][[1]]$y)], 2, segments[(i*2 - 1),ii], 2)
          endpts[endpts_count,] <- c(con_xy_ext[[segments[(i * 2 - 1),ii]]]$x[length(con_xy_ext[[segments[(i * 2 - 1),ii]]]$x)], con_xy_ext[[segments[(i * 2 - 1),ii]]]$y[length(con_xy_ext[[segments[(i * 2 - 1),ii]]]$y)], 2, segments[(i*2 - 1),ii], 2)
          endpts_count <- endpts_count + 1
        }
      }
      
      # for(ii in 1:(segCounter[(i * 2)] - 1)){ ### add the points before the start of the segment through the stagnation point
      if(segCounter[(i*2)] > 0){
        for(ii in 1:(segCounter[(i * 2)])){
          endpts[endpts_count,] <- c(con_xy_ext[[segments[(i * 2),ii]]]$x[1], con_xy_ext[[segments[(i * 2),ii]]]$y[1], 2, segments[(i*2),ii], 1)
          endpts_count <- endpts_count + 1
          endpts[endpts_count,] <- c(con_xy_ext[[segments[(i * 2),ii]]]$x[length(con_xy_ext[[segments[(i * 2),ii]]]$x)], con_xy_ext[[segments[(i * 2),ii]]]$y[length(con_xy_ext[[segments[(i * 2),ii]]]$y)], 2, segments[(i*2),ii], 2)
          endpts_count <- endpts_count + 1
        }
      }
      
      ### inspiration from:
      ### https://stackoverflow.com/questions/67497664/how-do-i-sort-points-in-clockwise-order-in-r-with-respect-to-the-center
      ### calc centre of set of endpoints; pick them in clockwise order
      centre_xy <- c(sum(endpts[,1])/length(endpts), sum(endpts[,2])/length(endpts))
      cen_angle2 <- atan2(centre_xy[2] - endpts[,2], centre_xy[1] - endpts[,1])
      ### now pick the largest angle and proceed toward the smallest, sort endpts descending
      ### https://www.statology.org/sort-matrix-in-r/
      ord1 <- order(cen_angle2, decreasing = TRUE)
      sorted_endpts <- endpts[ord1, ]
      
      flagOrd <- matrix(0, nrow(sorted_endpts), 1)
      
      for(ii in 1:nrow(sorted_endpts)){
        ## first verify that the other endpoint has not been processed
        if (flagOrd[ii,1] == 0){
        
          if(sorted_endpts[ii,3] == 1){ # use con_xy
            if(ii == 1){ # first row
              if(sorted_endpts[ii,5] == 1){ # in order
                pts <- cbind(con_xy[[sorted_endpts[ii,4]]][,1], con_xy[[sorted_endpts[ii,4]]][,2])
              }else{ # reverse order
                seqIndx <- seq(length(con_xy[[sorted_endpts[ii,4]]][,1]),1,by=-1)
                pts <- cbind(con_xy[[sorted_endpts[ii,4]]][seqIndx,1], con_xy[[sorted_endpts[ii,4]]][seqIndx,2])
              }
            }else{
              if(sorted_endpts[ii,5] == 1){
                pts <- rbind(pts, cbind(con_xy[[sorted_endpts[ii,4]]][,1], con_xy[[sorted_endpts[ii,4]]][,2]))
              }else{ # reverse order
                seqIndx <- seq(length(con_xy[[sorted_endpts[ii,4]]][,1]),1,by=-1)
                pts <- rbind(pts, cbind(con_xy[[sorted_endpts[ii,4]]][seqIndx,1], con_xy[[sorted_endpts[ii,4]]][seqIndx,2]))
              }
            }
          }else{ # use con_xy_ext
            if(ii == 1){ # first row
              if(sorted_endpts[ii,5] == 1){
                pts <- cbind(con_xy_ext[[sorted_endpts[ii,4]]]$x, con_xy_ext[[sorted_endpts[ii,4]]]$y)
              }else{
                seqIndx <- seq(length(con_xy_ext[[sorted_endpts[ii,4]]]$x),1,by=-1)
                pts <- cbind(con_xy_ext[[sorted_endpts[ii,4]]]$x[seqIndx], con_xy_ext[[sorted_endpts[ii,4]]]$y[seqIndx])
              }
            }else{
              if(sorted_endpts[ii,5] == 1){
                pts <- rbind(pts, cbind(con_xy_ext[[sorted_endpts[ii,4]]]$x, con_xy_ext[[sorted_endpts[ii,4]]]$y))
              }else{
                seqIndx <- seq(length(con_xy_ext[[sorted_endpts[ii,4]]]$x),1,by=-1)
                pts <- rbind(pts, cbind(con_xy_ext[[sorted_endpts[ii,4]]]$x[seqIndx], con_xy_ext[[sorted_endpts[ii,4]]]$y[seqIndx]))
              }
            }
          }
          
          ### flag rows that are endpoints of the same set of points that was just processed
          fIndx <- which(sorted_endpts[,4] == sorted_endpts[ii,4])
          fIndx2 <- which(sorted_endpts[,3] == sorted_endpts[ii,3])
          fIndx3 <- c(fIndx,fIndx2)
          # print(fIndx3)
          for(iii in 1:length(fIndx3)){
            if(length(which(fIndx3 == fIndx3[iii])) > 1){
              flagOrd[fIndx3[iii],1] <- 1
            }
          }
          # print(flagOrd)
        }
      }
      
      capzones[[counter]] <- pts
      counter <- counter + 1
    }
  }
  
  returnList <- list()
  returnList[[1]] <- psi3
  returnList[[2]] <- capzones
  
  return(returnList)
}