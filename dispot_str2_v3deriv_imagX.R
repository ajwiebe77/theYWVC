#### dispot_str2_v3deriv_imagX.R
#
# Andrew J. Wiebe, 9 May 2022
#
# Return the derivative of equation 19 in Nagheli et al (2020) for complex numbers representing the
# dimensionless discharge potential and the stream function (i.e., return value of equation 19 for a given zD location.
# This is for case 2 in Nagheli et al. (2020), where the wedge-shaped aquifer is bounded by 
# head boundaries on two sides and is open (semi-infinite) on the third side (outer arc).
#
# Parameters:
#    wells(:,1) = real coordinate for each well location (real(zwD))
#    wells(:,2) = complex coordinate for each well location (imag(zwD))
#    wells(:,3) = pumping rate for each well (>0 for injection, < 0 for extraction)
#    consts(1,1) = r_wedge, i.e., the length of the wedge in the actual system
#    consts(2,1) = alpha, i.e., the angle of the triangular aquifer
#    consts(3,1) = beta, i.e., the angle between the ambient aquifer flow direction and the x-axis
#    consts(4,1) = K, i.e., the hydraulic conductivity of the aquifer
#    consts(5,1) = b, i.e., the thickness of the aquifer
#    consts(6,1) = q0, i.e., regional uniform flow per unit width
#    zD          = a complex number representing the real and imaginary coordinates at which to evaluate the equation
# 
#    Coded for R on 7 June 2022
#    Based on Octave code

dispot_str2_v3deriv_imagX <- function(wells, consts, zD){
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
  
  zwD <- cbind(complex(real = wells[,1], imaginary = wells[,2]))
  
  s <- 0
  
  for (w in 1:nrow(wells)){
    # s <- s + (QwD[w]/(2*pi)) * ((pi * power(zD,t - 1) / (alpha *(power(zD,t) - power(zwD[w],t)))) - (pi * power(zD,t - 1) / (alpha * (power(zD,t) - power(Conj(zwD[w]),t)))))
    s <- s + (QwD[w]/(2*pi)) * ((pi * zD^(t - 1) / (alpha *((zD^t) - (zwD[w])^t))) - (pi * zD^(t - 1) / (alpha * ((zD^t) - (Conj(zwD[w]))^t))))
  }
  
  dOdz <- - q0D * exp(complex(real = 0, imaginary = - beta)) + s
  return(dOdz)
}