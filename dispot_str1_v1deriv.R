## ---------------------------------------------------------------------
## dispot_str1_v1deriv.m
#
# Andrew J. Wiebe, 12 May 2022; modifies dispot_str2_v3.m, 9 May 2022
#
# Return the derivative of equation 4 in Nagheli et al (2020) for complex numbers representing the
# dimensionless discharge potential and the stream function (i.e., return value of equation 7 for a given zD location.
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
#    zD          = a pair of values representing the real and imaginary components of a complex number representing a point at which to evaluate the gradient of the complex potential
# 
# Coded for R on 1 June 2022

dispot_str1_v1deriv <- function(wells, consts, zD){
  if(length(consts) == 5){
    r <- consts[1]
    beta <- consts[2]
    K <- consts[3]
    b <- consts[4]
    q0 <- consts[5]
  }else{
    r <- consts[1]
    beta <- consts[3]
    K <- consts[4]
    b <- consts[5]
    q0 <- consts[6]
  }
  
  q0D <- q0 / (K * b)
  
  Qw <- wells[,3]
  QwD <- Qw / (K * b * r)
  
  #zwD <- complex(real = wells[,1], imaginary = wells[,2]) # this is not a column vector as expected
  zwD <- cbind(complex(real = wells[,1], imaginary = wells[,2]))
  
  sum1 <- 0
  
  for (w in 1:nrow(wells)){
    sum1 <- sum1 + (QwD[w]/(2*pi)) * (1 /(zD - zwD[w]) - 1/(zD - Conj(zwD[w])))
  }
  
  dOdz <- - q0D * exp(complex(real = 0, imaginary = - beta)) + sum1
}