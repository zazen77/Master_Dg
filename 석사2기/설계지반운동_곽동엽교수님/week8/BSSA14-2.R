BSSA14 <- function(M,Vs30,Rjb,ftype,region,dz1){
# BSSA14 GMPE
# Mag : Moment Magnitude
# VS30 : VS30 (m/s)
# Rjb : Joyner-Boore Distance (km)
# RS : Reverse-Fault tag
# NS : Normal-Fault tag
# SS : Strike Slip-Fault tag
# region (0:CA,TW,NZ,global, 1:CHtur, 3:italyjapan)
# dz1 : differential from the average basin-depth predictor
# prd : target period for IM (0 = PGA, -1 = PGV)

# rm(list=ls(all=TRUE)) 
  
# Example input parameters
# M <- 6
# Vs30 <- 760
# Rjb <- 1
# ftype = 2
# region = 0
# dz1 = 0

  # Read Coefficients 
  coef <- read.csv('./GMM/BSSA14/coeffs_BSSA14.csv')
  
  period <- coef$Period
  # For linear site response -- ln(Flin)  
  clin <- coef$c
  Vc <- coef$Vc
  Vref <- 760
  
  # For nonlinear site amplification  --ln(Fnl)
  f1 = coef$f1
  f3 = coef$f3
  f4 = coef$f4
  f5 = coef$f5

  # For distance terms
  c1 = coef$c1
  c2 = coef$c2
  c3 = coef$c3
  h = coef$h
  Mref = coef$Mref
  Rref = coef$Rref

  # For magnitude scaling and focal mechanism terms
  e0 = coef$e0
  e1 = coef$e1
  e2 = coef$e2
  e3 = coef$e3
  e4 = coef$e4
  e5 = coef$e5
  e6 = coef$e6
  Mh = coef$Mh

  # For aleatory uncertainties
  R1 = coef$R1
  R2 = coef$R2
  DfR = coef$DfR
  DfV = coef$DfV
  V1 = coef$V1
  V2 = coef$V2
  phi1 = coef$phi1
  phi2 = coef$phi2
  tau1 = coef$tau1
  tau2 = coef$tau2

  # region Dc3 , period region code (0:CA,TW,NZ,global, 1:CHtur, 3:italyjapan)
  regionCHtur = coef$Dc3.ChinaTurkey
  regionJPit = coef$Dc3.ItalyJapan

  # rbasin depth, period, f6, f7, f7/f6
  f6 = coef$f6
  f7 = coef$f7
  f7f6 = f7 / f6
  
  if (ftype == 1) {
    NS = 1
    RS = 0
    SS = 0
    U = 0
  } else if (ftype == 2) {
    NS = 0
    RS = 1
    SS = 0
    U = 0
  } else if (ftype == 3) {
    NS = 0
    RS = 0
    SS = 1
    U = 0
  } else {
    NS = 0
    RS = 0
    SS = 0
    U = 1
  }
    
  ##### Calculate pga4nl ######  validated
  if(region == 1){
    Dc3 = regionCHtur
  } else if (region == 3){
    Dc3 = regionJPit
  } else {
    Dc3 = array(0,dim=length(regionCHtur))
  }

  R = sqrt(Rjb^2 + (h[1]^2))
  Fd1 = (c1[1] + c2[1] * (M - Mref[1])) * log(R / Rref[1]) + (c3[1] + Dc3[1]) * (R - Rref[1])
  if(M <= Mh[1]){
    Fm1 = e0[1] * U + e1[1] * SS + e2[1] * NS + e3[1] * RS + e4[1] * (M - Mh[1]) + e5[1] * ((M - Mh[1])^2)
  } else {
    Fm1 = e0[1] * U + e1[1] * SS + e2[1] * NS + e3[1] * RS + e6[1] * (M - Mh[1])
  }

  pga4nl = exp(Fm1 + Fd1)

  ###### Calculate SA for various Ts, PGA, and PGV ######
  
  # Magnitude Term - validated
  l = length(coef$Mh)
  Fm <-c()
  for(i in 1:l){
    if (M <= Mh[i]) {
      Fm[i] = e0[i] * U + e1[i] * SS + e2[i] * NS + e3[i] * RS + e4[i] * (M - Mh[i]) + e5[i] * (M - Mh[i])^2
    } else {
      Fm[i] = e0[i] * U + e1[i] * SS + e2[i] * NS + e3[i] * RS + e6[i] * (M - Mh[i])
    }
  }
  
  # Distance Term - validated
  R = sqrt(Rjb ^ 2 + (h ^ 2)) 
  if (region == 0) { #global, ca, tw, nz
      Fp = (c1 + c2 * (M - Mref)) * log(R / Rref) + (c3 + Dc3) * (R - Rref)
  } else if (region == 3 | region == 5) { #ch, tur
    Fp = (c1 + c2 * (M - Mref)) * log(R / Rref) + (c3 + Dc3chtur) * (R - Rref)
  } else if (region == 1 | region == 4) { #jp, italy
    Fp = (c1 + c2 * (M - Mref)) * log(R / Rref) + (c3 + Dc3jpit) * (R - Rref)
  } else {
    Fp = (c1 + c2 * (M - Mref)) * log(R / Rref) + (c3 + Dc3) * (R - Rref)
  }

    # Basin depth term - validated
    # unknown basin depth input will be treated as 0, i.e. centered model  
    if(dz1 == 0 && period < 0.65){
      Fbd1 = 0   
    } else if(period >= 0.65 && dz1 <= f7f6){
      Fbd1 = f6 * dz1
    } else if(period >= 0.65 && dz1 > f7f6){
      Fbd1 = f7
    } else {
      Fbd1 =  0
    }
  
    # Linear Site Response - validated
  Flin <-c()
  for(i in 1:l){
    if(Vs30 <= Vc[i]){
      Flin[i] = clin[i] * log(Vs30 / Vref) #if Vs30 = Vref = 760, Flin = 0 (in ln units)    
    } else {
      Flin[i] = clin[i] * log(Vc[i] / Vref)
    }
  }

    # Nonlinear Site Term - validated
    if(Vs30 >= Vref){   #if Vs30 = Vref = 760, Fnl = 0 (in ln units)
      Fnl = 0
    } else {
      f2 = f4 * (exp(f5 * (min(Vs30, Vref) - 360)) - exp(f5 * (Vref - 360)))
      Fnl = f1 + f2 * log((pga4nl + f3) / f3)
    }

    # Aleatory Uncertainty
#     if(M <= 4.5){
#       tauM = tau1    
#     } else if(M > 4.5 && M<5.5){
#       tauM = tau1+(tau2-tau1)*(M-4.5)
#     } else {
#       tauM = tau2
#     }
# 
#     if(M <= 4.5){
#       phiM = phi1
#     } else if(M > 4.5 && M<5.5){
#       phiM= phi1+(phi2-phi1)*(M-4.5)
#     } else {
#       phiM = phi2
#     }
# 
#     if(Rjb <= R1){
#       phiMR =  phiM
#     } else if(Rjb > R1 && Rjb <= R2){
#       phiMR = phiM + DfR*(log(Rjb/R1)/log(R2/R1))
#     } else {
#       phiMR = phiM + DfR
#     }
#   
#     if(Vs30 >= V2){
#       phiMRV = phiMR
#     } else if(Vs30 >= V1 && Vs30 < V2){
#       phiMRV = phiMR - DfV*(log(V2/Vs30)/log(V2/V1))
#     } else {
#       phiMRV = phiMR - DfV
#     }

    #### Results of the model - validated
    # Median SA
    lnY = Fp + Fm + Fbd1 + Fnl + Flin # for cases except Vs30 = Vref = 760 m/s
    expY = exp(lnY)
    out = cbind(coef$Period, round(lnY,digits=4), round(expY,digits=4)) # in g exp
    colnames(out) = c('T (sec)','SA (ln g)','SA (exp g)')
    
    # Site amplification
    # siteterm = Fnl + Flin + Fbd1 # for cases except Vs30 = Vref = 760 m/s
  
    # Aleatory uncertainties
    # tau = tauM # inter-event uncertainty
    # phi = phiMRV # intra-event uncertainty
    # sigma = sqrt(tau^2+phi^2) # combined uncertainty
    return(as.data.frame(out))
    
  }

  # idx <- which(prd == period)
  # lnYi <- lnY[idx]


JEJB <- BSSA14(4.1,689.8,71.2,0,0,0)
#YOCB <- BSSA14(4.1,537,87.8,0,0,0)
#PHA2 <- BSSA14(4.1,750.2,47.7,0,0,0)
#USN2 <- BSSA14(4.1,731.2,86.6,0,0,0)
#SMP36 <- BSSA14(4.1,642.4,38.6,0,0,0)
#SMP106 <- BSSA14(4.1,399,46.2,0,0,0)
#SMP261 <- BSSA14(4.1,607.7,76.2,0,0,0)
#JEJB2 <- BSSA14(2.5,689.8,26.4,0,0,0)
#USN22 <- BSSA14(2.5,731.2,8.8,0,0,0)
SMP73 <- BSSA14(2.5,400,11.4,0,0,0)
#SMP185 <- BSSA14(2.5,731.2,8.8,0,0,0)

#BSSA14(3.4,760,61.65672,0,0,0)

