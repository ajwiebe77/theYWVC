#### dispot_str1_v1.R
#
# Andrew J. Wiebe, 12 May 2022; modifies dispot_str2_v3.m, 9 May 2022
#
# Calculate the discharge potential and equipotential lines for Equation 4 in Nagheli et al (2020)
# for complex numbers representing a given zD location.
# This is for case 1 in Nagheli et al. (2020), a semi-infinite aquifer with one constant head boundary and 
# the other three sides open.
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
  # Based on the Octave code version

dispot_str1_v1 <- function(wells, consts, x){
  zDx <- x[,1]
  zDy <- x[,2]
  phi0D <- x[,3]
  
  r <- consts[1]
  beta <- consts[2]
  K <- consts[3]
  b <- consts[4]
  q0 <- consts[5]
  q0D <- q0 / (K * b)
  
  Qw <- wells[,3]
  QwD <- Qw / (K * b * r)
  
  zD <- cbind(complex(real = zDx, imaginary = zDy))
  zwD <- cbind(complex(real = wells[,1], imaginary = wells[,2]))
  
  comOmegaDval <- matrix(0.0,nrow(x),1)
  
  for (n in 1:nrow(x)){
    sum1 <- 0
    for (w in 1:nrow(wells)){
      sum1 <- sum1 + (QwD[w]/(2*pi)) * log((zD[n] - zwD[w])/(zD[n] - Conj(zwD[w])))
    }
  
    comOmegaDval[n] <- phi0D[n] - q0D * zD[n] * exp(complex(real = 0, imaginary = - beta)) + sum1
  }

  
  return(comOmegaDval)
}
