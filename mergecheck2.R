### mergecheck.R
### Andrew J. Wiebe, 12 May 2023
#
# Objective: Given two sets of points that denote segments, determine whether it is reasonable to merge the two 
# based on the shortest distance between endpoints and the angles of lines through the first or last few points of each
#
# Tip = start of series of (x,y) coordinates; Tail = end of series
#
# Consider adding code: the line segment between the two points to be joined should also have an acceptable angle
#
# Parameters:
#    temp_pts - points for segment 1
#    temp_pts2 - points for segment 2
#    toldist - distance tolerance for a merger
#    indxskip - number of points to skip when drawing a line through the tip or tail of the segment
#    threshAng - angle difference threshold in radians that is a minimum value for an acceptable merger
#    flagexttail - flag indicating whether the tail of the first segment is to be extended (TRUE) or is already merged (FALSE)
#
# Return value: 0 - merge not recommended
#               1 - recommended merge tip to tip
#               2 - recommended merge tip to tail
#               3 - recommended merge tail to tip
#               4 - recommended merge tail to tail
#
#

mergecheck2 <- function(temp_pts, temp_pts2, toldist, indxskip, threshAng, flagexttail){

  mergetype <- 0
  rangeMin <- 10 * pi # dummy value
  rangeMax <- -10 * pi # dummy value
  
  x1 <- temp_pts[1,1]
  y1 <- temp_pts[1,2]
  x2 <- temp_pts[nrow(temp_pts),1]
  y2 <- temp_pts[nrow(temp_pts),2]
  x3 <- temp_pts2[1,1]
  y3 <- temp_pts2[1,2]
  x4 <- temp_pts2[nrow(temp_pts2),1]
  y4 <- temp_pts2[nrow(temp_pts2),2]
  
  ### compare segment endpoints
  dist1 <- sqrt((x3 - x1)^2 + (y3 - y1)^2) # tip to tip
  dist2 <- sqrt((x4 - x1)^2 + (y4 - y1)^2) # tip to tail
  dist3 <- sqrt((x3 - x2)^2 + (y3 - y2)^2) # tail to tip
  dist4 <- sqrt((x4 - x2)^2 + (y4 - y2)^2) # tail to tail
  
  if(dist1 < toldist && dist1 <= dist2 && flagexttail == FALSE){ # merge tip to tip
    i1 <- indxskip + 1 # 27 Jan 2023
    i2 <- indxskip + 1 # 27 Jan 2023
    theta1 <- atan2((temp_pts[i1,2] - temp_pts[(i1-indxskip),2]),(temp_pts[i1,1] - temp_pts[(i1-indxskip),1])) # 12 May 2023 - calc with respect to "origin" between trial points, i.e., non-origin point minus origin point
    theta2 <- atan2((temp_pts2[i2,2] - temp_pts2[(i2-indxskip),2]),(temp_pts2[i2,1] - temp_pts2[(i2-indxskip),1])) # 12 May 2023
    
    # ensure angles are between 0 and 2 * pi
    if(theta1 < 0){
      theta1 <- 2 * pi + theta1 # 8 May 2023
    }
    
    if(theta2 < 0){
      theta2 <- 2 * pi + theta2 # 8 May 2023
    }
    
    # calculate expected angle at 180 degrees away from theta1, compare with theta2
    if(theta1 > (pi + threshAng)){
      expectedAng <- theta1 - pi
      rangeMin <- expectedAng - threshAng
      rangeMax <- expectedAng + threshAng
    }else if(theta1 < (pi - threshAng)){
      expectedAng <- theta1 + pi
      rangeMin <- expectedAng - threshAng
      rangeMax <- expectedAng + threshAng
    }else{ # expected angle is around zero degrees and two cases must be considered (Quad 1 and Quad 4)
      if(theta2 <= 2 * threshAng){ # Quad 1 (+x,+y)
        expectedAng <- theta1 - pi
        rangeMin <- expectedAng - threshAng # could be negative
        rangeMax <- expectedAng + threshAng
      }else if(theta2 >= (2 * pi - 2 * threshAng)){
        expectedAng <- theta1 + pi
        rangeMin <- expectedAng - threshAng
        rangeMax <- expectedAng + threshAng
      }
    }
    
    if(theta2 >= rangeMin && theta2 <= rangeMax){
      mergetype <- 1
      # print(paste("tip-tip", x1,y1,x2,y2,x3,y3,x4,y4, theta1, theta2))
    }
    
  }else if(dist2 < toldist && dist2 <= dist1 && flagexttail == FALSE){ # merge tip to tail # 3 May 2023
    i1 <- indxskip + 1 # 27 Jan 2023
    i2 <- nrow(temp_pts2)
    theta1 <- atan2((temp_pts[i1,2] - temp_pts[(i1-indxskip),2]),(temp_pts[i1,1] - temp_pts[(i1-indxskip),1])) # 12 May 2023 - subtract the "origin"
    theta2 <- atan2((temp_pts2[(i2-indxskip),2] - temp_pts2[i2,2]),(temp_pts2[(i2-indxskip),1] - temp_pts2[i2,1])) # 12 May 2023
    
    # ensure angles are between 0 and 2 * pi
    if(theta1 < 0){
      theta1 <- 2 * pi + theta1 # 8 May 2023
    }
    
    if(theta2 < 0){
      theta2 <- 2 * pi + theta2 # 8 May 2023
    }
    
    # calculate expected angle at 180 degrees away from theta1, compare with theta2
    if(theta1 > (pi + threshAng)){
      expectedAng <- theta1 - pi
      rangeMin <- expectedAng - threshAng
      rangeMax <- expectedAng + threshAng
    }else if(theta1 < (pi - threshAng)){
      expectedAng <- theta1 + pi
      rangeMin <- expectedAng - threshAng
      rangeMax <- expectedAng + threshAng
    }else{ # expected angle is around zero degrees and two cases must be considered (Quad 1 and Quad 4)
      if(theta2 <= 2 * threshAng){ # Quad 1 (+x,+y)
        expectedAng <- theta1 - pi
        rangeMin <- expectedAng - threshAng # could be negative
        rangeMax <- expectedAng + threshAng
      }else if(theta2 >= (2 * pi - 2 * threshAng)){
        expectedAng <- theta1 + pi
        rangeMin <- expectedAng - threshAng
        rangeMax <- expectedAng + threshAng
      }
    }
    
    if(theta2 >= rangeMin && theta2 <= rangeMax){
      mergetype <- 2
      # print(paste("tip-tail",x1,y1,x2,y2,x3,y3,x4,y4, theta1, theta2))
    }
    
  }else if(dist3 < toldist && dist3 <= dist4 && flagexttail == TRUE){ # merge tail to tip
    i1 <- nrow(temp_pts)
    i2 <- indxskip + 1 # 27 Jan 2023
    theta1 <- atan2((temp_pts[(i1-indxskip),2] - temp_pts[i1,2]),(temp_pts[(i1-indxskip),1] - temp_pts[i1,1])) # 3 May 2023
    theta2 <- atan2((temp_pts2[i2,2] - temp_pts2[(i2-indxskip),2]), (temp_pts2[i2,1] - temp_pts2[(i2-indxskip),1])) # 27 Jan 2023
    
    # ensure angles are between 0 and 2 * pi
    if(theta1 < 0){
      theta1 <- 2 * pi + theta1 # 8 May 2023
    }
    
    if(theta2 < 0){
      theta2 <- 2 * pi + theta2 # 8 May 2023
    }
    
    # calculate expected angle at 180 degrees away from theta1, compare with theta2
    if(theta1 > (pi + threshAng)){
      expectedAng <- theta1 - pi
      rangeMin <- expectedAng - threshAng
      rangeMax <- expectedAng + threshAng
    }else if(theta1 < (pi - threshAng)){
      expectedAng <- theta1 + pi
      rangeMin <- expectedAng - threshAng
      rangeMax <- expectedAng + threshAng
    }else{ # expected angle is around zero degrees and two cases must be considered (Quad 1 and Quad 4)
      if(theta2 <= 2 * threshAng){ # Quad 1 (+x,+y)
        expectedAng <- theta1 - pi
        rangeMin <- expectedAng - threshAng # could be negative
        rangeMax <- expectedAng + threshAng
      }else if(theta2 >= (2 * pi - 2 * threshAng)){
        expectedAng <- theta1 + pi
        rangeMin <- expectedAng - threshAng
        rangeMax <- expectedAng + threshAng
      }
    }
    
    if(theta2 >= rangeMin && theta2 <= rangeMax){
      mergetype <- 3
      # print(paste("tail-tip", x1,y1,x2,y2,x3,y3,x4,y4, theta1, theta2))
    }
    
  }else if(dist4 < toldist && dist4 <= dist3 && flagexttail == TRUE){ # merge tail to tail
    i1 <- nrow(temp_pts)
    i2 <- nrow(temp_pts2)
    theta1 <- atan2((temp_pts[(i1-indxskip),2] - temp_pts[i1,2]), (temp_pts[(i1-indxskip),1] - temp_pts[i1,1])) # 3 Apr 2023
    theta2 <- atan2((temp_pts2[(i2-indxskip),2] - temp_pts2[i2,2]), (temp_pts2[(i2-indxskip),1] - temp_pts2[i2,1])) # 27 Jan 2023
    
    # ensure angles are between 0 and 2 * pi
    if(theta1 < 0){
      theta1 <- 2 * pi + theta1 # 8 May 2023
    }
    
    if(theta2 < 0){
      theta2 <- 2 * pi + theta2 # 8 May 2023
    }
    
    # calculate expected angle at 180 degrees away from theta1, compare with theta2
    if(theta1 > (pi + threshAng)){
      expectedAng <- theta1 - pi
      rangeMin <- expectedAng - threshAng
      rangeMax <- expectedAng + threshAng
    }else if(theta1 < (pi - threshAng)){
      expectedAng <- theta1 + pi
      rangeMin <- expectedAng - threshAng
      rangeMax <- expectedAng + threshAng
    }else{ # expected angle is around zero degrees and two cases must be considered (Quad 1 and Quad 4)
      if(theta2 <= 2 * threshAng){ # Quad 1 (+x,+y)
        expectedAng <- theta1 - pi
        rangeMin <- expectedAng - threshAng # could be negative
        rangeMax <- expectedAng + threshAng
      }else if(theta2 >= (2 * pi - 2 * threshAng)){
        expectedAng <- theta1 + pi
        rangeMin <- expectedAng - threshAng
        rangeMax <- expectedAng + threshAng
      }
    }
    
    if(theta2 >= rangeMin && theta2 <= rangeMax){
      mergetype <- 4
      # print(paste("tail-tail", x1,y1,x2,y2,x3,y3,x4,y4, theta1, theta2))
    }
    # }
  }
  
  return(mergetype)
}