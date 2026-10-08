rm(list=ls())
library(xtable)

wd <- "~"
setwd(wd)

## Parameters
iS <- 10000
iN <- 1000
iM <- 1000
vK <- 2:6
## Fifty logarithmically spaced actual grid sizes in [21, 2001]
vGridSize <- 2L * round(
  (exp(seq(log(21), log(2001), length.out = 50L)) - 1) / 2
) + 1L

stopifnot(
  length(vGridSize) == 50L,
  length(unique(vGridSize)) == 50L,
  min(vGridSize) == 21L,
  max(vGridSize) == 2001L,
  all(vGridSize %% 2L == 1L)
)

## Half-width needed by the existing symmetric-grid construction
vHalfGrid <- (vGridSize - 1L) %/% 2L
lA <- list()
for (p in 1:length(vGridSize)) {
  lA[[p]] <- (-vGridSize[p]:vGridSize[p]) / (1 + vGridSize[p])
}
vP <- as.numeric(unlist(lapply(lA,length)))
dTheta <- 0.5

fnMCSummary <- function(aEst, dTheta, vGridSize, vK) {
  vGridSize <- length(vGridSize)
  iK <- length(vK)
  
  mMSE <- matrix(
    NA_real_, vGridSize, iK,
    dimnames = list(
      paste0("P=", seq_len(vGridSize)),
      paste0("K=", vK)
    )
  )
  mMCSE <- mMSE
  
  for (p in seq_len(vGridSize)) {
    for (k in seq_len(iK)) {
      vSquaredError <- (aEst[, k, p] - dTheta)^2
      
      mMSE[p, k] <- mean(vSquaredError)
      mMCSE[p, k] <- sd(vSquaredError) / sqrt(length(vSquaredError))
    }
  }
  
  list(MSE = mMSE, MCSE = mMCSE)
}

fnBenchmarkMSE <- function(vEstimate, dTheta) {
  vSquaredError <- (as.numeric(vEstimate) - dTheta)^2
  
  c(
    MSE = mean(vSquaredError),
    MCSE = sd(vSquaredError) / sqrt(length(vSquaredError))
  )
}

## Read data II
vEst1.II <- read.csv(paste(wd,"Est_MA_II_N100_M100.csv",sep=""),header = T)[,-1]
vEst2.II <- read.csv(paste(wd,"Est_MA_II_N100_M1000.csv",sep=""),header = T)[,-1]
vEst3.II <- read.csv(paste(wd,"Est_MA_II_N1000_M1000.csv",sep=""),header = T)[,-1]

dBias1.II <- mean(vEst1.II)-dTheta
dVar1.II <- var(vEst1.II)
dMSE1.II <- dBias1.II^2+dVar1.II

dBias2.II <- mean(vEst2.II)-dTheta
dVar2.II <- var(vEst2.II)
dMSE2.II <- dBias2.II^2+dVar2.II

dBias3.II <- mean(vEst3.II)-dTheta
dVar3.II <- var(vEst3.II)
dMSE3.II <- dBias3.II^2+dVar3.II

## Read data MLE
vEst1.mle <- read.csv(paste(wd,"Est_MA_MLE_N100.csv",sep=""),header = T)[,-1]
vEst2.mle <- read.csv(paste(wd,"Est_MA_MLE_N1000.csv",sep=""),header = T)[,-1]

dBias1.mle <- mean(vEst1.mle)-dTheta
dVar1.mle <- var(vEst1.mle)
dMSE1.mle <- dBias1.mle^2+dVar1.mle

dBias2.mle <- mean(vEst2.mle)-dTheta
dVar2.mle <- var(vEst2.mle)
dMSE2.mle <- dBias2.mle^2+dVar2.mle

## Read data SMM
vEst1.SMM <- read.csv(paste(wd,"Est_MA_SMM_N100_M100.csv",sep=""),header = T)[,-1]
vEst2.SMM <- read.csv(paste(wd,"Est_MA_SMM_N100_M1000.csv",sep=""),header = T)[,-1]
vEst3.SMM <- read.csv(paste(wd,"Est_MA_SMM_N1000_M1000.csv",sep=""),header = T)[,-1]

dBias1.SMM <- mean(vEst1.SMM)-dTheta
dVar1.SMM <- var(vEst1.SMM)
dMSE1.SMM <- dBias1.SMM^2+dVar1.SMM

dBias2.SMM <- mean(vEst2.SMM)-dTheta
dVar2.SMM <- var(vEst2.SMM)
dMSE2.SMM <- dBias2.SMM^2+dVar2.SMM

dBias3.SMM <- mean(vEst3.SMM)-dTheta
dVar3.SMM <- var(vEst3.SMM)
dMSE3.SMM <- dBias3.SMM^2+dVar3.SMM

wd <- "~"
setwd(wd)

lRun1 <- readRDS("Est_MA_SNPMD_N100_M100_logP.rds")
lRun2 <- readRDS("Est_MA_SNPMD_N100_M1000_logP.rds")
lRun3 <- readRDS("Est_MA_SNPMD_N1000_M1000_logP.rds")

vGridSize <- lRun1$grid_size
vK <- lRun1$degree

stopifnot(
  identical(vGridSize, lRun2$grid_size),
  identical(vGridSize, lRun3$grid_size)
)

lStats1 <- fnMCSummary(lRun1$estimates, dTheta, vGridSize, vK)
lStats2 <- fnMCSummary(lRun2$estimates, dTheta, vGridSize, vK)
lStats3 <- fnMCSummary(lRun3$estimates, dTheta, vGridSize, vK)

pdf("MC_exp_MA_logP.pdf",width=10,height=10)
par(mfrow = c(3, 1), 
    mar = c(0, 0, 0, 0),  # Reduced margins: bottom=5, left=4, top=2, right=1
    oma = c(4, 4, 1, 1),  # Outer margin at bottom for legend
    xpd = NA)

plot(lStats1$MSE[,1],type="l",lwd=2,ylim=c(0,0.5),xlab="",ylab="",xaxt="n")
for (k in 1:5) {
  lines(lStats1$MSE[,k],lwd=2,lty=k)
}
segments(x0=1,y0=dMSE1.II,x1=50,y1=dMSE1.II,lwd=2,col="darkgray")
segments(x0=1,y0=dMSE1.mle,x1=50,y1=dMSE1.mle,lwd=2,lty=2,col="darkgray")
segments(x0=1,y0=dMSE1.SMM,x1=50,y1=dMSE1.SMM,lwd=2,lty=3,col="darkgray")

legend("topright",legend=c(paste0("K=",2:6),"Ind. Inf.","MLE","SMM"),
       lty = c(1:5, 1, 2, 3),col = c(rep("black",5), rep("darkgray",3)),lwd = 2,
       bty = "n", cex = 0.8)

plot(lStats2$MSE[,1],type="l",lwd=2,ylim=c(0,0.5),xlab="",ylab="",xaxt="n")
for (k in 1:5) {
  lines(lStats2$MSE[,k],lwd=2,lty=k)
}
segments(x0=1,y0=dMSE2.II,x1=50,y1=dMSE2.II,lwd=2,col="darkgray")
segments(x0=1,y0=dMSE1.mle,x1=50,y1=dMSE1.mle,lwd=2,lty=2,col="darkgray")
segments(x0=1,y0=dMSE2.SMM,x1=50,y1=dMSE2.SMM,lwd=2,lty=3,col="darkgray")

plot(lStats3$MSE[,1],type="l",lwd=2,ylim=c(0,0.5),xlab="",ylab="",xaxt="n")
axis(side = 1, at = which(vGridSize%in%vTicks),
     labels=vTicks)
for (k in 1:5) {
  lines(lStats3$MSE[,k],lwd=2,lty=k)
}
segments(x0=1,y0=dMSE3.II,x1=50,y1=dMSE3.II,lwd=2,col="darkgray")
segments(x0=1,y0=dMSE2.mle,x1=50,y1=dMSE2.mle,lwd=2,lty=2,col="darkgray")
segments(x0=1,y0=dMSE3.SMM,x1=50,y1=dMSE3.SMM,lwd=2,lty=3,col="darkgray")

# Reset layout to default
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, oma = c(0, 0, 0, 0), xpd = FALSE)
dev.off()

dfMSE1 <- as.data.frame(lStats1$MSE[c(1,5,10,15,20,25,30,35,40,45,50),1:5])
dfMSE2 <- as.data.frame(lStats2$MSE[c(1,5,10,15,20,25,30,35,40,45,50),1:5])
dfMSE3 <- as.data.frame(lStats3$MSE[c(1,5,10,15,20,25,30,35,40,45,50),1:5])

## ============================================================
## Multi-panel LaTeX MSE table
## ============================================================

## Indices already used to construct dfMSE1, dfMSE2 and dfMSE3
vTableIndex <- c(1,5,10,15,20,25,30,35,40,45,50)
vTableP <- vGridSize[vTableIndex]

## Add P labels and the appropriate benchmarks to each panel
fnCreateMSEPanel <- function(
    dfMSE,
    dMSE.SMM,
    dMSE.II,
    dMSE.MLE
) {
  stopifnot(
    nrow(dfMSE) == length(vTableP),
    ncol(dfMSE) >= 5L
  )
  
  dfPanel <- data.frame(
    P = vTableP,
    dfMSE[, 1:5, drop = FALSE],
    SMM = rep(dMSE.SMM, length(vTableP)),
    II = rep(dMSE.II, length(vTableP)),
    MLE = rep(dMSE.MLE, length(vTableP)),
    check.names = FALSE
  )
  
  names(dfPanel) <- c(
    "$P$",
    paste0("$K=", 2:6, "$"),
    "SMM",
    "II",
    "MLE"
  )
  
  dfPanel
}

## Panel A: N = 100, M = 100
dfPanelA <- fnCreateMSEPanel(
  dfMSE = dfMSE1,
  dMSE.SMM = dMSE1.SMM,
  dMSE.II = dMSE1.II,
  dMSE.MLE = dMSE1.mle
)

## Panel B: N = 100, M = 1000
## The MLE depends only on N, so the N = 100 MLE is used.
dfPanelB <- fnCreateMSEPanel(
  dfMSE = dfMSE2,
  dMSE.SMM = dMSE2.SMM,
  dMSE.II = dMSE2.II,
  dMSE.MLE = dMSE1.mle
)

## Panel C: N = 1000, M = 1000
dfPanelC <- fnCreateMSEPanel(
  dfMSE = dfMSE3,
  dMSE.SMM = dMSE3.SMM,
  dMSE.II = dMSE3.II,
  dMSE.MLE = dMSE2.mle
)

## Combine the three panels
dfTableMSE <- rbind(dfPanelA, dfPanelB, dfPanelC)
rownames(dfTableMSE) <- NULL

iRowsPerPanel <- nrow(dfPanelA)
iNumberColumns <- ncol(dfTableMSE)

## Insert panel headings into the LaTeX tabular
lPanelHeaders <- list(
  pos = list(
    0,
    iRowsPerPanel,
    2L * iRowsPerPanel
  ),
  command = c(
    paste0(
      "\\multicolumn{", iNumberColumns,
      "}{c}{$N=100$, $M=100$}",
      "\\\\ \n"
    ),
    paste0(
      "\\addlinespace\n",
      "\\multicolumn{", iNumberColumns,
      "}{c}{$N=100$, $M=1000$}",
      "\\\\ \n"
    ),
    paste0(
      "\\addlinespace\n",
      "\\multicolumn{", iNumberColumns,
      "}{c}{$N=1000$, $M=1000$}",
      "\\\\ \n"
    )
  )
)

## Create xtable object
lTableMSE <- xtable(
  dfTableMSE,
  caption = paste(
    "Monte Carlo mean squared errors for the SNPMD estimator of $\\theta$",
    "in the MA(1) experiment with a logarithmically-spaced grid.",
    "The table reports results for different representative grid sizes $P$,",
    "polynomial degree $P$, and sample sizes $N,M$. The performance of the",
    " SNPMD estimator is compared with II, KF-MLE, and SMM."
  ),
  label = "tab:MA_MSE_logP",
  align = c(
    "l",                     # row-name column used internally by xtable
    "r",                     # P
    rep("r", 8L)             # K=2,...,6, II, MLE and SMM
  ),
  digits = c(
    0,                       # row-name column
    0,                       # P
    rep(4, 8L)               # MSE values
  )
)

## Export table
pTableFile <- file.path(wd, "Table_MA_MSE_logP.tex")

print(
  lTableMSE,
  file = pTableFile,
  type = "latex",
  floating = TRUE,
  table.placement = "!htbp",
  caption.placement = "top",
  include.rownames = FALSE,
  sanitize.text.function = identity,
  sanitize.colnames.function = identity,
  add.to.row = lPanelHeaders,
  booktabs = TRUE,
  hline.after = c(-1, 0, nrow(dfTableMSE)),
  size = "scriptsize"
)

cat("LaTeX table saved to:\n", pTableFile, "\n")
