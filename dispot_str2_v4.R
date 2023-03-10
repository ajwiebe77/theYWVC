#### dispot_str2_v4.m
#
# Andrew J. Wiebe, 13 May 2022; modifies dispot_str2_v3_imagX.m, 9 May 2022; modifies dispot_str2_v2.m, 28 Apr 2022 (which seems to be a draft that looks more like case 3b than case 2); modifies dispot_str2.m, 10 Mar 2022
#
# Solve equation 16 in Nagheli et al (2020) for complex numbers representing the
# dimensionless discharge potential and the stream function. This is for case 2 in the
# Nagheli et al. (2020) paper, where the wedge-shaped aquifer is bounded by 
# head boundaries on two sides and is open (semi-infinite) on the third side (outer arc).
#
# Parameters:
  #    wells(:,1) = real coordinate for each well location (real(zwD))
  #    wells(:,2) = complex coordinate for each well location (imag(zwD))
  #    wells(:,3) = pumping rate for each well (>0 for injection, < 0 for extraction)
  #    consts(1,1) - r_wedge, i.e., the length of the wedge in the actual system
  #    consts(2,1) - alpha, i.e., the angle of the triangular aquifer
  #    consts(3,1) - beta, i.e., the angle between the ambient aquifer flow direction and the x-axis
  #    consts(4,1) - K, i.e., the hydraulic conductivity of the aquifer
  #    consts(5,1) - b, i.e., the thickness of the aquifer
  #    consts(6,1) - q0, i.e., regional uniform flow per unit width
  #    x(:,1) = zDx, list of dimensionless x coordinates at which to evaluate the equation
  #    x(:,2) = zDy, list of dimensionless y coordinates at which to evaluate the equation
  #    x(:,3) = phiD0 - the initial dimensionless discharge potential with no pumping, for each (real(zD)=zDx,imag(zD)=zDy) pair
  
  # 
  # Returns a complex number corresponding to the zD location. The real part of the
  #    number is the dimensionless discharge potential, and the imaginary part is the stream function.
  #
  # Coded for R on 7 June 2022
  # Based on Octave code
  
dispot_str2_v4 <- function(wells, consts, x){
  zDx <- x[,1]
  zDy <- x[,2]
  phi0D <- x[,3]
  
  r_wedge <- consts[1]
  alpha <- consts[2]
  beta <- consts[3]
  K <- consts[4]
  b <- consts[5]
  q0 <- consts[6]
  q0D <- q0 / (K * b)
  t <- pi / alpha
  
  Qw <- wells[,3]
  QwD <- Qw / (K * b * r_wedge)
  
  zD <- cbind(complex(real = zDx, imaginary = zDy))
  zwD <- cbind(complex(real = wells[,1], imaginary = wells[,2]))
  
  comOmegaDval <- matrix(0.0,nrow(x),1)
  
  for (n in 1:nrow(x)){
    sum1 <- 0
    for (w in 1:nrow(wells)){
      mag_zwDw <- abs(zwD[w])
      # sum1 <- sum1 + (QwD[w]/(2*pi)) * log((power(zD[n] / mag_zwDw,pi/alpha) - power(zwD[w] / mag_zwDw,pi/alpha))/(power(zD[n] / mag_zwDw, pi/alpha) - power(Conj(zwD[w]) / mag_zwDw,pi/alpha)))
      sum1 <- sum1 + (QwD[w]/(2*pi)) * log(((zD[n] / mag_zwDw)^t - (zwD[w] / mag_zwDw)^t)/((zD[n] / mag_zwDw)^t - (Conj(zwD[w]) / mag_zwDw)^t))
    }
  
    comOmegaDval[n] <- phi0D[n] - q0D * zD[n] * exp(complex(real = 0, imaginary = - beta)) + sum1
  }
  
  return(comOmegaDval)
}