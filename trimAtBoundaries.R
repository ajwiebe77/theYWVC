### holzbecher_v3.R
### Andrew J. Wiebe, 9 Dec 2022
#
# Parameters:
#     list1 - a list of sets of points
#     case - the case from Nagheli et al. (2020)
#     option - the option from Nagheli et al. (2020) demonstration examples, needed for triangular aquifers
#     alpha - the angle of the triangular aquifer
#     xy - boolean value denoting whether the x and y vectors in the list are labelled or not
# 
# Return values:
#    list2 <- a list of the sets of points that are within the aquifer


trimAtBoundaries <- function(list1, case, option, alpha, xy){
  list2 <- list1 
  
  if(case == 1){
    ptA <- c(0, 1)
    ptB <- c(0, 0)
    ptC <- c(1, 0)
    ptD <- c(1, 1)
    boundary <- rbind(ptA, ptB, ptC, ptD)
    
    for(i in 1:length(list1)){
      if(xy == TRUE){
        for(ii in 1:length(list1[[i]]$x)){ # for each point
          if(list1[[i]]$x[ii] > 0){
            # fails if x > 1
            # fails if y < 0
            # fails if x < 0
            # fails if y > 1
            if(list1[[i]]$x[ii] > 1.0 || list1[[i]]$y[ii] < 0 || list1[[i]]$x[ii] < 0 || list1[[i]]$y[ii] > 1){
              list2[[i]]$x[ii] <- NaN
              list2[[i]]$y[ii] <- NaN
            }
          }
        }
      }
    }
    
  }else if(case == 2 && (option == 1 || option == 2 || option == 3 || option == 4)){
    ptA <- c(cos(alpha), sin(alpha))
    ptB <- c(0, 0)
    ptC <- c(1, 0) # mirror image across y axis
    boundary <- rbind(ptA, ptB, ptC)
    outerarc <- cbind(cos(seq(0,alpha,alpha/10)), sin(seq(0,alpha,alpha/10)))
    
    for(i in 1:length(list1)){
      if(xy == TRUE){
        for(ii in 1:length(list1[[i]]$x)){ # for each point
          if(list1[[i]]$x[ii] > 0){
            # fails if x > 1
            # fails if y < 0
            # fails if y > outerarc
            if(list1[[i]]$x[ii] > 1.0 || list1[[i]]$y[ii] < 0 || list1[[i]]$y[ii] > sin(acos(list1[[i]]$x[ii]))){
              list2[[i]]$x[ii] <- NaN
              list2[[i]]$y[ii] <- NaN
            }
          }else if(list1[[i]]$x[ii] < 0){ # x < 0
            # fails if x < cos(alpha)
            # fails if y < x * tan(alpha)
            # fails if y > outerarc
            if(list1[[i]]$x[ii] < cos(alpha) || list1[[i]]$y[ii] < list1[[i]]$x[ii] * tan(alpha) || list1[[i]]$y[ii] > sin(acos(list1[[i]]$x[ii]))){
              list2[[i]]$x[ii] <- NaN
              list2[[i]]$y[ii] <- NaN
            }
          }else if(list1[[i]]$x[ii] == 0){
            # fails if y > 1
            # fails if y < 0
            if(list1[[i]]$y[ii] > 1 || list1[[i]]$y[ii] < 0){
              list2[[i]]$x[ii] <- NaN
              list2[[i]]$y[ii] <- NaN
            }
          }
        }
      }else{
        for(ii in 1:length(list1[[i]][,1])){ # for each point
          if(list1[[i]][ii,1] > 0){
            # fails if x > 1
            # fails if y < 0
            # fails if y > outerarc
            if(list1[[i]][ii,1] > 1.0 || list1[[i]][ii,2] < 0 || list1[[i]][ii,2] > sin(acos(list1[[i]][ii,1]))){
              list2[[i]][ii,1] <- NaN
              list2[[i]][ii,2] <- NaN
            }
          }else if(list1[[i]][ii,1] < 0){ # x < 0
            # fails if x < cos(alpha)
            # fails if y < x * tan(alpha)
            # fails if y > outerarc
            if(list1[[i]][ii,1] < cos(alpha) || list1[[i]][ii,2] < list1[[i]][ii,2] * tan(alpha) || list1[[i]][ii,2] > sin(acos(list1[[i]][ii,1]))){
              list2[[i]][ii,1] <- NaN
              list2[[i]][ii,2] <- NaN
            }
          }else if(list1[[i]][ii,1] == 0){
            # fails if y > 1
            # fails if y < 0
            if(list1[[i]][ii,2] > 1 || list1[[i]][ii,2] < 0){
              list2[[i]][ii,1] <- NaN
              list2[[i]][ii,2] <- NaN
            }
          }
        }
      }
    }
  }
  
  
  return(list2)
}