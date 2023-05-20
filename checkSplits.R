### checkSplits.R
### Andrew J. Wiebe, 14 May 2023
#
# Objective: Given two values and a two-column vector, check whether the two values are in the same row in the vector and return TRUE or FALSE.
#            Also deal with the case where one or both values are not present in the vector's columns.
#            Assume that each value is only in each column zero or one times.
# 
# Parameters:
#    vect - a two-column vector
#    val1 - first value to check
#    val2 - second value to check
#
#
# Return value:
#    TRUE - val1 and val2 are in the same row somewhere in vect
#    FALSE - val1 and val2 are not in the same row in vect
#

checkSplits <- function(vect, val1, val2){
  
  sameRow <- FALSE
  
  if(length(which(vect[,1] == val1)) >= 1 && length(which(vect[,2] == val2)) >= 1){ # val1 is in column 1 of vect and val2 is in column2 of vect
    if(which(vect[,1] == val1) == which(vect[,2] == val2)){
      sameRow <- TRUE
    }
  }
  
  if(length(which(vect[,2] == val1)) >= 1 && length(which(vect[,1] == val2)) >= 1){ # val1 is in column 1 of vect and val2 is in column2 of vect
    if(which(vect[,2] == val1) == which(vect[,1] == val2)){
      sameRow <- TRUE
    }
  }
  
  return(sameRow)
}
