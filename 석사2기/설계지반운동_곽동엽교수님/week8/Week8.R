rm(list=ls())

library(pracma)
library(dplyr)

TE <- function(x) { (beta*exp(-beta*(x-Mmin))) / (1 - exp(-beta*(Mmax-Mmin))) }
b <- 0.9
beta <- b*log(10)
Mmin <- 5
Mmax <- 7

P_M1 <- integrate(TE, 5, 5.5)$value
f_M1 <- P_M1/0.5

P_M2 <- integrate(TE, 5.5, 6)$value
fM2 <- P_M2/0.5

P_M3 <- integrate(TE, 6, 6.5)$value
fM3 <- P_M3/0.5

P_M4 <- integrate(TE, 6.5, 7)$value
fM4 <- P_M4/0.5




# 4. ----------------------------------------------------------------------

BSSA14 <- function(M, Vs30, Rjb, fault.type, region, z1, data.type) {
  
  library(data.table)
  # BSSA14 GMM
  
  # Recommendation
  # SS & RS : M = 3-8.5
  # NS : M = 3-7
  # Distance, Rjb : 0-400 km
  # Time-averaged Vs30 : 150-1500 m/s
  # Basin effect, z1 : 0-3 km
  # Main shock and afeter shock
  
  # M <- 6
  # Vs30 <- 760
  # Rjb <- 20
  # fault.type <- 0
  # region <- 0
  # z1 <- 0
  # data.type <- 0
  
  # Read Coefficients 
  # coef <- read.csv('./GMM/BSSA14/coeffs_BSSA14.csv')
  # coef <- read.csv('./GMM/BSSA14/05_eqs-070113eqs184m_suppl_es1_030915.csv')
  coef <- fread('BSSA14_Coefficients_071314_Revisedf4_071514.csv')
  period <- coef$period
  
  # For linear site response : ln(Flin)  
  clin <- coef$clin
  Vc <- coef$Vc
  Vref <- 760
  
  # For nonlinear site amplification : ln(Fnl)
  f1 <- coef$f1
  f3 <- coef$f3
  f4 <- coef$f4
  f5 <- coef$f5
  
  # For path terms
  c1 <- coef$c1
  c2 <- coef$c2
  c3 <- coef$c3
  h <- coef$h
  Mref <- coef$Mref
  Rref <- coef$Rref
  
  # For magnitude scaling and focal mechanism terms
  e0 <- coef$e0
  e1 <- coef$e1
  e2 <- coef$e2
  e3 <- coef$e3
  e4 <- coef$e4
  e5 <- coef$e5
  e6 <- coef$e6
  Mh <- coef$Mh
  
  # For aleatory uncertainties
  R1 <- coef$R1
  R2 <- coef$R2
  DfR <- coef$DfR
  DfV <- coef$DfV
  V1 <- coef$V1
  V2 <- coef$V2
  phi1 <- coef$phi1
  phi2 <- coef$phi2
  tau1 <- coef$tau1
  tau2 <- coef$tau2
  
  # region Dc3 , period region code (0:CA,TW,NZ,global, 1:CHtur, 3:italyjapan)
  regionChTu <- coef$Dc3.ChinaTurkey
  regionItJp <- coef$Dc3.ItalyJapan
  
  # rbasin depth, period, f6, f7, f7/f6
  f6 <- coef$f6
  f7 <- coef$f7
  f7f6 <- f7 / f6
  
  
  # Fault Type
  # 0:Unspecifeied, 1:NS, 2:RS, 3:SS
  if (fault.type == 1) {
    NS <- 1
    RS <- 0
    SS <- 0
    U <- 0
  } else if (fault.type == 2) {
    NS <- 0
    RS <- 1
    SS <- 0
    U <- 0
  } else if (fault.type == 3) {
    NS <- 0
    RS <- 0
    SS <- 1
    U <- 0
  } else {
    NS <- 0
    RS <- 0
    SS <- 0
    U <- 1
  }
  
  
  # Fe:Souece(e:Event), Fp:Path, Fs: Site effects
  # Epsilon:fractional number of standard deviations
  
  # Region
  # region Dc3, 0:CA,TW,NZ,global, 1:China Turkey, 3:Italy, Japan)
  if(region == 1){
    Dc3 <- regionChTu
  } else if (region == 3){
    Dc3 <- regionItJp
  } else {
    Dc3 = array(0,dim=length(regionChTu))
  }
  
  
  # Magnitude Term (Event)
  l = length(coef$Mh)
  Fe <- c()
  for(i in 1:l){
    if (M <= Mh[i]) {
      Fe[i] <- e0[i] * U + e1[i] * SS + e2[i] * NS + e3[i] * RS + e4[i] * (M - Mh[i]) + e5[i] * (M - Mh[i])^2
    } else {
      Fe[i] <- e0[i] * U + e1[i] * SS + e2[i] * NS + e3[i] * RS + e6[i] * (M - Mh[i])
    }
  }
  
  # Path Term (Distance)
  R <- sqrt(Rjb ^ 2 + (h ^ 2)) 
  Fp <- (c1 + c2 * (M - Mref)) * log(R / Rref) + (c3 + Dc3) * (R - Rref)
  
  # Site Term : Fs = Flin + Fnl + Fdz1
  
  # Linear Site Amplification
  Flin <-c()
  for(i in 1:l){
    if(Vs30 <= Vc[i]){
      Flin[i] <- clin[i] * log(Vs30 / Vref) 
    } else {
      Flin[i] <- clin[i] * log(Vc[i] / Vref)
    }
  }
  
  
  # PGAr for non-linear : for rock, Vs = 760 = Vref. Fs=0. So, lnY = Fe + Fp
  R.a <- sqrt(Rjb^2 + (h[1]^2))
  Fp.a <- (c1[2] + c2[2] * (M - Mref[2])) * log(R.a / Rref[2]) + (c3[2] + Dc3[2]) * (R.a - Rref[2])
  if(M <= Mh[2]){
    Fe.a <- e0[2] * U + e1[2] * SS + e2[2] * NS + e3[2] * RS + e4[2] * (M - Mh[2]) + e5[2] * ((M - Mh[2])^2)
  } else {
    Fe.a <- e0[2] * U + e1[2] * SS + e2[2] * NS + e3[2] * RS + e6[2] * (M - Mh[2])
  }
  pgar <- exp(Fe.a + Fp.a)
  
  # Nonlinear Site Amplification
  if(Vs30 >= Vref){   
    Fnl <- f1
  } else {
    f2 <- f4 * (exp(f5 * (min(Vs30, Vref) - 360)) - exp(f5 * (Vref - 360)))
    Fnl <- f1 + f2 * log((pgar + f3) / f3)
  }
  
  
  # Basin depth on GM amplitude
  
  # # basin depth z1, 0:Unknown 
  if(z1 == 0){
    dz1 <- 0
  } else {
    z1 <- z1
    
    # data.type, 0:Unlnown, 1:California, 2:Japan
    if(data.type == 1){
      mz1 <- (-7.15/4)*log((Vs30^4+570.94^4)/(1360^4+570.94^4))-log(1000)
      dz1 <- z1-mz1
    } else if (data.type == 1) {
      mz1 <- (-5.23/2)*log((Vs30^2+412.39^2)/(1360^2+412.39^2))-log(1000)
      dz1 <- z1-mz1
    } else {
      dz1 <- 0
    }
    
  }
  
  # unknown basin depth input will be treated as 0. i.e. Fdz1=0 for "centered" model
  if(dz1 == 0 && period < 0.65){
    Fdz1 <- 0   
  } else if(period >= 0.65 && dz1 <= f7f6){
    Fdz1 <- f6 * dz1
  } else if(period >= 0.65 && dz1 > f7f6){
    Fdz1 <- f7
  } else {
    Fdz1 <- 0
  }
  
  
  # Aleatory Uncertainty Function : sigma = sqrt(tau^2+phi^2)
  
  if(M <= 4.5){
    tau_M <- tau1
  } else if(M > 4.5 && M<5.5){
    tau_M <- tau1+(tau2-tau1)*(M-4.5)
  } else {
    tau_M <- tau2
  }
  
  if(M <= 4.5){
    phi_M <- phi1
  } else if(M > 4.5 && M<5.5){
    phi_M <- phi1+(phi2-phi1)*(M-4.5)
  } else {
    phi_M<- phi2
  }
  
  if(Rjb <= R1){
    phi_MR <- phi_M
  } else if(Rjb > R1 && Rjb <= R2){
    phi_MR<- phi_M + DfR*(log(Rjb/R1)/log(R2/R1))
  } else {
    phi_MR <- phi_M + DfR
  }
  
  if(Vs30 >= V2){
    phi_MRV = phi_MR
  } else if(Vs30 >= V1 && Vs30 < V2){
    phi_MRV <- phi_MR - DfV*(log(V2/Vs30)/log(V2/V1))
  } else {
    phi_MRV <- phi_MR - DfV
  }
  
  # Result
  Fs <- Flin + Fnl + Fdz1
  lnY <- Fp + Fe + Fs
  Y <- exp(lnY)
  
  Tau <- tau_M                # Between-event Uncertainty
  Phi <- phi_MRV              # Within-event Uncertainty
  sigma <- sqrt(Tau^2+Phi^2)  # Total Uncertainty
  
  result <- cbind("Period (s)" = coef$period,
                  "SA (ln g)" = round(lnY, 4),
                  "SA (g)" = round(Y, 4),
                  "STD (ln)" = round(sigma, 4))
  
  out <- cbind("SA (g)" = round(Y, 4),
               "STD (ln)" = round(sigma, 4))
  
  return(as.data.frame(result))
  # return(as.data.frame(out))
  
}

M.med <- c(5.25, 5.75, 6.25, 6.75)
R.jb <- c(17.205, 10.770, 10.770, 17.205, 26.000, 35.440,
          16.008, 10.308, 10.308, 16.008, 24.622, 34.004,
          12.806, 10.000, 10.000, 13.454, 21.471, 30.676,
          10.000, 10.000, 10.000, 10.050, 14.866, 15.620)

IMs <- list()

for(i in 1:length(M.med)) {
  
  for(j in 1:length(R.jb)){
    IMs[[j]] <- BSSA14(M.med[i], 760, R.jb[j], 3, 0, 0, 0)
  }
  
  return(IMs)
}


pga <- c()
std <- c()
for(i in 1:length(IMs)){
  
  pga[i] <- IMs[[i]]$`SA (ln g)`[2] 
  std[i] <- IMs[[i]]$`STD (ln)`[2]
  
}

out <- cbind("PGA (ln)" = pga,
             "Std (ln)" = std)

# 5. ----------------------------------------------------------------------

z <- 0.1    # PGA.test, g unit

# e: Epsilon = (ln(z) - PGA) / Std
e <- (log(z)- pga)/std
e

# Standard Normal dist.
# P(pga>z|M,Rjb)
P_pga <- 1-pnorm(e,mean = 0,sd = 1)

# P(M), P(h)
P_M <- c(P_M1, P_M2, P_M3, P_M4)
P_h <- rep(1/6, 6)

# Evaluate P(M)*P(h)*P(pga>z|M,Rjb)
P_M.h <- cbind(rep(P_M,6),P_h)

P.all <- data.frame(P_M.h, P_pga)
colnames(P.all) <- c("P(M)", "P(h)", "P(PGA>z|M,Rjb)")

P <- P_M.h[,1]*P_M.h[,2]*P_pga
sum(P)

lambda <- 0.03
v <- sum(P)*lambda

fx.poisson <- function(v,dt) {1-exp(-v*dt)}
dt <- 1   # 1 year
P.fin <- fx.poisson(v,dt)
P.fin
  
# 6. ----------------------------------------------------------------------

# Repeat above for value z : 0.001-1
z <- seq(0.001,1,0.001)
e <- list()
P_pga <- list()
P <- data.frame()
sig.P <- data.frame()

for (i in 1:length(z)) {
  
  e[[i]] <- (log(z[i]) - pga)/std
  P_pga[[i]] <- 1-pnorm(e[[i]],mean = 0,sd = 1)
  
  # P(M), P(h)
  P_M <- c(P_M1, P_M2, P_M3, P_M4)
  P_h <- rep(1/6, 6)
  
  # Evaluate P(M)*P(h)*P(pga>z|M,Rjb)
  P_M.h <- cbind(rep(P_M,6),P_h)
  
  lambda <- 0.03
  
  v <- sum(P_M.h[,1]*P_M.h[,2]*P_pga[[i]])*lambda
  sigma.P <- v/lambda
  
  fx.poisson <- function(v,dt) {1-exp(-v*dt)}
  dt <- 1   # 1 year
  P.poi <- fx.poisson(v,dt)
  P <- rbind(P,P.poi)
  sig.P <- rbind(sig.P,sigma.P)
}

P <- cbind(z,P, sig.P)
colnames(P) <- c('PGA','P','sigma.P')

# P$PGA <- round(P$PGA, 3)
# P$PGA[which(P$P == round(P.50, 3))]


# PGA with 10% prob. of exceedence in 50 years
dt <- 50
p <- 0.1
v.50 <- -(log(1-p))/dt
P.50 <- fx.poisson(v.50,1)



library(ggplot2)
breaks <- 10^(-10:10)
minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))

p <- ggplot(P,aes(PGA,P)) +
  theme_bw() +
  geom_line() +
  geom_segment(aes(x=0, xend=0.197, y=P.50, yend=P.50),color="red", size=1) +
  geom_segment(aes(x=0.197, xend=0.197, y=P.50, yend=0),color="red",arrow=arrow(), size=1) +
  labs(x = "PGA (g)", y='APE') +
  scale_y_log10(breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(breaks = breaks, minor_breaks = minor_breaks) +
  theme(axis.title = element_text(size=15),
        title = element_text(size=18),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15))
p

jpeg('6.jpg',width=3000,height=1800,res=350)
print(p)
dev.off()





# 7. ----------------------------------------------------------------------

# Repeat above for value z : 0.001-1
z <- seq(0.001,1,0.001)
e <- list()
P_pga <- list()
P <- data.frame()
sig.P <- data.frame()

for (i in 1:length(z)) {
  
  e[[i]] <- (log(z[i]) - pga)/std
  P_pga[[i]] <- 1-pnorm(e[[i]],mean = 0,sd = 1)
  
  # P(M), P(h)
  P_M <- c(P_M1, P_M2, P_M3, P_M4)
  P_h <- rep(1/6, 6)
  
  # Evaluate P(M)*P(h)*P(pga>z|M,Rjb)
  P_M.h <- cbind(rep(P_M,6),P_h)
  
  lambda <- 0.03
  
  v <- sum(P_M.h[,1]*P_M.h[,2]*P_pga[[i]])*lambda
  sigma.P <- v/lambda
  
  fx.poisson <- function(v,dt) {1-exp(-v*dt)}
  dt <- 1   # 1 year
  P.poi <- fx.poisson(v,dt)
  P <- rbind(P,P.poi)
  sig.P <- rbind(sig.P,sigma.P)
}

P <- cbind(z,P, sig.P)
colnames(P) <- c('PGA','P','sigma.P')


IMs <- BSSA14(6.75, 760, 10, 3, 0, 0, 0)
pga.bssa <- exp(IMs$`SA (ln g)`[2])
std.bssa <- exp(IMs$`SA (ln g)`[2]+IMs$`STD (ln)`[2])

pga.bssa <- round(pga.bssa, 3)
std.bssa <- round(std.bssa, 3)

P.med <- P$P[which(P$PGA == pga.bssa)]

a <- pga.bssa+std.bssa
P.std <- P$P[which(P$PGA == pga.bssa+std.bssa)]
P.std <- 6.385451e-06

# Return Period
dt <- 1
v.med <- -(log(1-P.med))/dt
v.std <- -(log(1-P.std))/dt

# Check for lambda
P.sig.med <- P$sigma.P[which(P$P == P.med)]
P.sig.std <- P$sigma.P[which(P$P == P.std)]
lambda.med <- v.med/P.sig.med
lambda.std <- v.std/P.sig.std


# Plot
p <- ggplot(P,aes(PGA,P)) +
  theme_bw() +
  geom_line() +
  # geom_segment(aes(x=0, xend=0.195, y=P.50, yend=P.50),color="red", size=1) +
  geom_segment(aes(x=pga.bssa, xend=0, y=P.med, yend=P.med, color="Median"),arrow=arrow(length = unit(0.5, 'cm')), size=1) +
  geom_segment(aes(x=pga.bssa, xend=pga.bssa, y=0, yend=P.med, color="Median"), size=1) +
  geom_segment(aes(x=pga.bssa+std.bssa, xend=0, y=P.std, yend=P.std, color="Median + Stand Deviation"),arrow=arrow(length = unit(0.5, 'cm')), size=1) +
  geom_segment(aes(x=pga.bssa+std.bssa, xend=pga.bssa+std.bssa, y=0, yend=P.std, color="Median + Stand Deviation"), size=1) +
  scale_color_manual("Type",values = c("Median"="blue","Median + Stand Deviation"="red")) +
  labs(x = "PGA (g)", y='APE') +
  scale_y_log10(breaks = breaks, minor_breaks = minor_breaks) +
  scale_x_log10(breaks = breaks, minor_breaks = minor_breaks) +
  theme(axis.title = element_text(size=15),
        title = element_text(size=18),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        legend.title = element_text(size = 10),
        legend.position = c(0.855,0.89),
        legend.text=element_text(size=10),
        legend.direction = 'vertical',
        legend.box.background = element_rect(color = "black",fill="white"),legend.box.margin = margin(0.1,0.1,0.1,0.1,"cm"))
p

jpeg('7.jpg',width=3000,height=1800,res=350)
print(p)
dev.off()





















