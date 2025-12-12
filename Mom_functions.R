####################################
#---Functions to compute moments---#
####################################

fnMom1 <- function(x) {                # first central moment
  mom1 <- mean(x-mean(x))
  return(mom1)
}

fnMom2 <- function(x) {                # second central moment
  mom2 <- mean((x-mean(x))^2)
  return(mom2)
}

fnACOV1 <- function(x) {               # autocovariance functions
  x.ts <- as.ts(x)
  acv <- acf(x.ts,lag.max=1,type="covariance",plot = FALSE)
  return(acv$acf[2])
}

fnACOV2 <- function(x) {               # autocovariance functions
  x.ts <- as.ts(x)
  acv <- acf(x.ts,lag.max=2,type="covariance",plot = FALSE)
  return(acv$acf[3])
}

fnACOV3 <- function(x) {               # autocovariance functions
  x.ts <- as.ts(x)
  acv <- acf(x.ts,lag.max=3,type="covariance",plot = FALSE)
  return(acv$acf[4])
}

fcACOR <- function(x) {									# autocorrelation function
  cor(x[-1],x[1:(length(x)-1)])
}

fcAC21 <- function(x) {									# E(Y_t{^2},Y_t+1)
  cor(x[-1],x[1:(length(x)-1)]^2)
}

fcAC12 <- function(x) {									# E(Y_t,Y_t+1{^2})
  cor(x[-1]^2,x[1:(length(x)-1)])
}

fcAC2 <- function(x) {									# E(Y_t{^2},Y_t+1{^2})
  cor(x[-1]^2,x[1:(length(x)-1)]^2)
}

fcAbsAC <- function(x) {									# abs autocorrelation function
  cor(abs(x[-1]),abs(x[1:(length(x)-1)]))
}

fcAC22 <- function(x) {									# E(Y_t{^2},Y_t+2{^2})
  cor(x[-(1:2)]^2,x[1:(length(x)-2)]^2)
}

fcAC32 <- function(x) {									# E(Y_t{^2},Y_t+3{^2})
  cor(x[-(1:3)]^2,x[1:(length(x)-3)]^2)
}

fcAC42 <- function(x) {									# E(Y_t{^2},Y_t+4{^2})
  cor(x[-(1:4)]^2,x[1:(length(x)-4)]^2)
}

fcAC52 <- function(x) {									# E(Y_t{^2},Y_t+5{^2})
  cor(x[-(1:5)]^2,x[1:(length(x)-5)]^2)
}

fcAbsAC2 <- function(x) {									# E|Y_t,Y_t+2|
  cor(abs(x[-(1:2)]),abs(x[1:(length(x)-2)]))
}

fcAbsAC3 <- function(x) {									# E|Y_t,Y_t+3|
  cor(abs(x[-(1:3)]),abs(x[1:(length(x)-3)]))
}

fcAbsAC4 <- function(x) {									# E|Y_t,Y_t+4|
  cor(abs(x[-(1:4)]),abs(x[1:(length(x)-4)]))
}

fcAbsAC5 <- function(x) {									# E|Y_t,Y_t+5|
  cor(abs(x[-(1:5)]),abs(x[1:(length(x)-5)]))
}

fcAC62 <- function(x) {									# E(Y_t{^2},Y_t+6{^2})
  cor(x[-(1:6)]^2,x[1:(length(x)-6)]^2)
}

fcAC72 <- function(x) {									# E(Y_t{^2},Y_t+7{^2})
  cor(x[-(1:7)]^2,x[1:(length(x)-7)]^2)
}

fcAC82 <- function(x) {									# E(Y_t{^2},Y_t+8{^2})
  cor(x[-(1:8)]^2,x[1:(length(x)-8)]^2)
}

fcAC92 <- function(x) {									# E(Y_t{^2},Y_t+9{^2})
  cor(x[-(1:9)]^2,x[1:(length(x)-9)]^2)
}

fcAC102 <- function(x) {									# E(Y_t{^2},Y_t+10{^2})
  cor(x[-(1:10)]^2,x[1:(length(x)-10)]^2)
}

fcAbsAC6 <- function(x) {									# E|Y_t,Y_t+6|
  cor(abs(x[-(1:6)]),abs(x[1:(length(x)-6)]))
}

fcAbsAC7 <- function(x) {									# E|Y_t,Y_t+7|
  cor(abs(x[-(1:7)]),abs(x[1:(length(x)-7)]))
}

fcAbsAC8 <- function(x) {									# E|Y_t,Y_t+8|
  cor(abs(x[-(1:8)]),abs(x[1:(length(x)-8)]))
}

fcAbsAC9 <- function(x) {									# E|Y_t,Y_t+9|
  cor(abs(x[-(1:9)]),abs(x[1:(length(x)-9)]))
}

fcAbsAC10 <- function(x) {									# E|Y_t,Y_t+10|
  cor(abs(x[-(1:10)]),abs(x[1:(length(x)-10)]))
}

fcAC152 <- function(x) {						 			# E(Y_t{^2},Y_t+15{^2})
  cor(x[-(1:15)]^2,x[1:(length(x)-15)]^2)
}

fcAbsAC15 <- function(x) {									# E|Y_t,Y_t+15|
  cor(abs(x[-(1:15)]),abs(x[1:(length(x)-15)]))
}

fcAC202 <- function(x) {						 			# E(Y_t{^2},Y_t+20{^2})
  cor(x[-(1:20)]^2,x[1:(length(x)-20)]^2)
}

fcAbsAC20 <- function(x) {									# E|Y_t,Y_t+20|
  cor(abs(x[-(1:20)]),abs(x[1:(length(x)-20)]))
}

fcAC252 <- function(x) {						 			# E(Y_t{^2},Y_t+25{^2})
  cor(x[-(1:25)]^2,x[1:(length(x)-25)]^2)
}

fcAbsAC25 <- function(x) {									# E|Y_t,Y_t+25|
  cor(abs(x[-(1:25)]),abs(x[1:(length(x)-25)]))
}

##############################################
# Remove transient and subsample simulations #
##############################################

fnF <- function(z,t) {z[-(1:t)]}
fnSub <- function(y,s) {y[1:s]}

