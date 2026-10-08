rm(list = ls())

wd <- "~/SSSUP Dropbox/Mario Martinoli/SNPMD/Submissions/Journals/JEBO/Revisions/Robustness_checks"
setwd(wd)

iS <- 10000
iN <- 100#0
iM <- 1000
iK <- 6
iP <- 100
iLim <- iP
vA <- (-iLim:iLim) / (1 + iLim)
iR <- length(vA)
dTheta <- 0.5

## Load files
# load(paste(wd,"/output/Est_MA_SNPMD_N100_M100_P201_het.RData",sep=""))
load(paste(wd,"/output/Est_MA_SNPMD_N100_M100_P201_het_cond.RData",sep=""))
mDist1 <- mDist
mfp1 <- mfp
vEst1 <- vEst
rm(mDist,mfp,vEst)
# load(paste(wd,"/output/Est_MA_SNPMD_N100_M1000_P201_het.RData",sep=""))
load(paste(wd,"/output/Est_MA_SNPMD_N100_M1000_P201_het_cond.RData",sep=""))
mDist2 <- mDist
mfp2 <- mfp
vEst2 <- vEst
rm(mDist,mfp,vEst)
# load(paste(wd,"/output/Est_MA_SNPMD_N1000_M1000_P201_het.RData",sep=""))
load(paste(wd,"/output/Est_MA_SNPMD_N1000_M1000_P201_het_cond.RData",sep=""))
mDist3 <- mDist
mfp3 <- mfp
vEst3 <- vEst
rm(mDist,mfp,vEst)

## Split the 10,000 replications
vReference <- 1:5000
vEvaluation <- 5001:10000

## Independent approximation of f(theta)
vF1.reference <- colMeans(mDist1[vReference, , drop = FALSE])
vF2.reference <- colMeans(mDist2[vReference, , drop = FALSE])
vF3.reference <- colMeans(mDist3[vReference, , drop = FALSE])

## Sampling errors: rows = replications; columns = theta values
mError1 <- sweep(
  mDist1[vEvaluation, , drop = FALSE],
  MARGIN = 2,
  STATS = vF1.reference,
  FUN = "-"
)

mError2 <- sweep(
  mDist2[vEvaluation, , drop = FALSE],
  MARGIN = 2,
  STATS = vF2.reference,
  FUN = "-"
)

mError3 <- sweep(
  mDist3[vEvaluation, , drop = FALSE],
  MARGIN = 2,
  STATS = vF3.reference,
  FUN = "-"
)

## Mean errors
vMeanError1 <- colMeans(mError1)
vMeanError2 <- colMeans(mError2)
vMeanError3 <- colMeans(mError3)

## Include uncertainty from estimation of the reference function
vSE1 <- sqrt(
  apply(mDist1[vEvaluation, , drop = FALSE], 2, var) / length(vEvaluation) +
    apply(mDist1[vReference, , drop = FALSE], 2, var) / length(vReference)
)

vSE2 <- sqrt(
  apply(mDist2[vEvaluation, , drop = FALSE], 2, var) / length(vEvaluation) +
    apply(mDist2[vReference, , drop = FALSE], 2, var) / length(vReference)
)

vSE3 <- sqrt(
  apply(mDist3[vEvaluation, , drop = FALSE], 2, var) / length(vEvaluation) +
    apply(mDist3[vReference, , drop = FALSE], 2, var) / length(vReference)
)

## Simultaneous Bonferroni confidence bands
dAlpha <- 0.05
dCritical <- qnorm(1 - dAlpha/2) #/ (2 * iR))
vLow1 <- vMeanError1 - dCritical * vSE1
vUp1 <- vMeanError1 + dCritical * vSE1
vLow2 <- vMeanError2 - dCritical * vSE1
vUp2 <- vMeanError2 + dCritical * vSE1
vLow3 <- vMeanError3 - dCritical * vSE1
vUp3 <- vMeanError3 + dCritical * vSE1

## Estimate Sigma
mSigma1 <- cov(mError1)
mCorrelation1 <- cov2cor(mSigma1)
mSigma2 <- cov(mError2)
mCorrelation2 <- cov2cor(mSigma2)
mSigma2 <- cov(mError2)
mCorrelation2 <- cov2cor(mSigma2)
mSigma3 <- cov(mError3)
mCorrelation3 <- cov2cor(mSigma3)

vEigenvalues1 <- eigen(
  (mSigma1 + t(mSigma1)) / 2,
  symmetric = TRUE,
  only.values = TRUE
)$values

vEigenvalues2 <- eigen(
  (mSigma2 + t(mSigma2)) / 2,
  symmetric = TRUE,
  only.values = TRUE
)$values

vEigenvalues3 <- eigen(
  (mSigma3 + t(mSigma3)) / 2,
  symmetric = TRUE,
  only.values = TRUE
)$values

dLambdaMin1 <- min(vEigenvalues1)
dLambdaMax1 <- max(vEigenvalues1)
dTrace1 <- sum(vEigenvalues1)
dConditionNumber1 <- dLambdaMax1 / dLambdaMin1

dLambdaMin2 <- min(vEigenvalues2)
dLambdaMax2 <- max(vEigenvalues2)
dTrace2 <- sum(vEigenvalues2)
dConditionNumber2 <- dLambdaMax2 / dLambdaMin2

dLambdaMin3 <- min(vEigenvalues3)
dLambdaMax3 <- max(vEigenvalues3)
dTrace3 <- sum(vEigenvalues3)
dConditionNumber3 <- dLambdaMax3 / dLambdaMin3

mOffDiagonal1 <- row(mCorrelation1) != col(mCorrelation1)
dMeanAbsoluteCorrelation1 <- mean(
  abs(mCorrelation1[mOffDiagonal1])
)

mOffDiagonal2 <- row(mCorrelation2) != col(mCorrelation2)
dMeanAbsoluteCorrelation2 <- mean(
  abs(mCorrelation2[mOffDiagonal2])
)

mOffDiagonal3 <- row(mCorrelation3) != col(mCorrelation3)
dMeanAbsoluteCorrelation3 <- mean(
  abs(mCorrelation3[mOffDiagonal3])
)

c(
  lambda_min = dLambdaMin1,
  lambda_max = dLambdaMax1,
  trace = dTrace1,
  condition_number = dConditionNumber1,
  mean_absolute_correlation = dMeanAbsoluteCorrelation1
)

c(
  lambda_min = dLambdaMin2,
  lambda_max = dLambdaMax2,
  trace = dTrace2,
  condition_number = dConditionNumber2,
  mean_absolute_correlation = dMeanAbsoluteCorrelation2
)

c(
  lambda_min = dLambdaMin3,
  lambda_max = dLambdaMax3,
  trace = dTrace3,
  condition_number = dConditionNumber3,
  mean_absolute_correlation = dMeanAbsoluteCorrelation3
)

# pdf("MC_exp_hetero_std.pdf",width=10,height=10)
pdf("MC_exp_hetero_cond.pdf",width=10,height=10)
par(mfrow = c(3, 3), 
    mar = c(0, 1, 0, 1),  # Reduced margins: bottom=5, left=4, top=2, right=1
    oma = c(4, 4, 1, 1),  # Outer margin at bottom for legend
    xpd = NA)

## 1. Mean-zero condition
ylim <- range(vLow1, vUp1)

plot(
  vA, vMeanError1,
  type = "l", lwd = 2,
  ylim = ylim,
  xlab = "",
  ylab = expression(bar(epsilon)(theta)),
  xaxt = "n",
  col = adjustcolor("gray", alpha.f = 0.75)
)
polygon(
  c(vA, rev(vA)),
  c(vLow1, rev(vUp1)),
  col = adjustcolor("gray", alpha.f = 0.25),
  border = NA
)
lines(vA, vMeanError1, lwd = 2)
segments(-1,0,1,0,lwd=2,lty=2)

ylim <- range(vLow2, vUp2)

plot(
  vA, vMeanError2,
  type = "l", lwd = 2,
  ylim = ylim,
  xlab = "",
  ylab = "",
  xaxt = "n",
  col = adjustcolor("gray", alpha.f = 0.75),
)
polygon(
  c(vA, rev(vA)),
  c(vLow2, rev(vUp2)),
  col = adjustcolor("gray", alpha.f = 0.25),
  border = NA
)
lines(vA, vMeanError2, lwd = 2)
# abline(h = 0, lty = 2)
segments(-1,0,1,0,lwd=2,lty=2)

ylim <- range(vLow3, vUp3)

plot(
  vA, vMeanError3,
  type = "l", lwd = 2,
  ylim = ylim,
  xlab = "",
  ylab = "",
  xaxt = "n",
  col = adjustcolor("gray", alpha.f = 0.75),
)
polygon(
  c(vA, rev(vA)),
  c(vLow3, rev(vUp3)),
  col = adjustcolor("gray", alpha.f = 0.25),
  border = NA
)
lines(vA, vMeanError3, lwd = 2)
segments(-1,0,1,0,lwd=2,lty=2)

## 2. Heteroskedasticity
plot(
  vA, diag(mSigma1),
  type = "l", lwd = 2,
  ylab = expression(hat(Sigma)),
  xaxt = "n"
)

plot(
  vA, diag(mSigma2),
  type = "l", lwd = 2,
  xlab = "",
  ylab = "",
  xaxt = "n"
)

plot(
  vA, diag(mSigma3),
  type = "l", lwd = 2,
  xlab = "",
  ylab = "",
  xaxt = "n",
)

## 3. Correlation across parameter values
image(
  vA, vA, mCorrelation1,
  zlim = c(-1, 1),
  col = hcl.colors(100, "Blue-Red 3"),
  xlab = "",
  ylab = expression(theta[j])
)
contour(vA, vA, mCorrelation1, add = TRUE, drawlabels = FALSE)

image(
  vA, vA, mCorrelation2,
  zlim = c(-1, 1),
  col = hcl.colors(100, "Blue-Red 3"),
  xlab = expression(theta[i]),
  ylab = ""
)
contour(vA, vA, mCorrelation2, add = TRUE, drawlabels = FALSE)

image(
  vA, vA, mCorrelation3,
  zlim = c(-1, 1),
  col = hcl.colors(100, "Blue-Red 3"),
  xlab = "",
  ylab = ""
)
contour(vA, vA, mCorrelation3, add = TRUE, drawlabels = FALSE)

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, oma = c(0, 0, 0, 0), xpd = FALSE)
dev.off()
