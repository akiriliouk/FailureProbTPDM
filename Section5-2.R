source("functions.R")

############################################################
####### Data pre-processing & plotting #######
############################################################

load("Wind.RData")
d <- nrow(coord)
n <- nrow(gust)

xlon <- sort(unique(coord$LON))
ylat <- sort(unique(coord$LAT))
ylatlim <- c(floor(min(ylat)), ceiling(max(ylat)))
xlonlim <- c(floor(min(xlon)), ceiling(max(xlon)))
indsea <- c(1,2,3,4,5, 6, 7, 9, 11, 14, 21, 22, 23, 24, 25, 26, 27)
indland <- c(1:d)[-indsea]
LonLat <- coord[,2:3]
d1 <- length(indsea)
d2 <- length(indland)

map("world", col= "white", fill = T, bg = "white", xlim = c(3.3,7.2), ylim = c(50.775, 53.45))
title("d = 35 stations", cex.main = 2)
points(LonLat[indsea,], pch = 19, col = "chocolate", cex = 1.5)
points(LonLat[indland,], pch = 19, col = "Forestgreen", cex = 1.5)

# data organized as sea stations (17) and then inland stations (18)
coord <- coord[c(indsea,indland),]

dataorig <- cbind(gust[,indsea+1], gust[,indland+1])
dataoland <- dataorig[,(d1+1):d]
dataosea <- dataorig[,1:d1]

# standardized to alpha = 2
data <- apply(dataorig, 2, function(i) (-log(Fgpd(i)))^(-1/2))
dataland <- data[,(d1+1):d]
datasea <- data[,1:d1]


########################################################
####### Chi and TPDM #######
############################################################

Sigma <- tailDepMatrix(data, alpha = 2, qu = 0.95, mest = FALSE)$Sigma
Sigmasea <- tailDepMatrix(datasea, alpha = 2, qu = 0.95, mest = FALSE)$Sigma
Sigmaland <- tailDepMatrix(dataland, alpha = 2, qu = 0.95, mest = FALSE)$Sigma

chis <- chiPairs(data, coord$LAT, coord$LON)
chisea <- chiPairs(data[,1:d1], coord$LAT[1:d1], coord$LON[1:d1])
chiland <- chiPairs(data[,(d1+1):d], coord$LAT[(d1+1):d], coord$LON[(d1+1):d])

pdf("WINDchi.pdf")
par(cex.lab=2,cex.axis=2,cex.main=2,mar=c(5,5.5,4,2))
plot(chis, pch = 19, xlab = 'distance (in km)', ylab = expression(paste(chi[jk])), main = expression(paste("Tail dependence coefficients ", chi[jk])), ylim = c(0.3, 1))
points(chisea, col = "brown", pch = 19)
points(chiland, col = "Forestgreen", pch = 19)
tmpsea <- loess.smooth(chisea[,1], chisea[,2])
tmpland <- loess.smooth(chiland[,1], chiland[,2])
lines(tmpsea, lwd = 3, col = "brown")
lines(tmpland, lwd = 3, col = "Forestgreen")
legend("topright", legend = c("coast", "inland"), col = c("brown", "Forestgreen"), lwd = c(3,3), cex = 2)
dev.off()

tpdms <- cbind(chis[,1], Sigma[lower.tri(Sigma)])
tpdmsea <- cbind(chisea[,1], Sigmasea[lower.tri(Sigmasea)])
tpdmland <- cbind(chiland[,1], Sigmaland[lower.tri(Sigmaland)])

pdf("WINDtpdm.pdf")
par(cex.lab=2,cex.axis=2,cex.main=2,mar=c(5,5.5,4,2))
plot(tpdms, pch = 19, xlab = 'distance (in km)', ylab = expression(paste(sigma[jk])), main = expression(paste("TPDM entries ", sigma[jk])), ylim = c(0.625, 1))
points(tpdmsea, col = "brown", pch = 19)
points(tpdmland, col = "ForestGreen", pch = 19)
tmpsea <- loess.smooth(tpdmsea[,1], tpdmsea[,2])
tmpland <- loess.smooth(tpdmland[,1], tpdmland[,2])
lines(tmpsea, lwd = 3, col = "brown")
lines(tmpland, lwd = 3, col = "Forestgreen")
legend("topright", legend = c("coast", "land"), col = c("brown", "Forestgreen"), lwd = c(3,3), cex = 2)
dev.off()


############################################################

x1 <- 1
x2 <- 1.2
x3 <- 1.3

x1sum <- 0.8
x2sum <- 0.9
x3sum <- 1

nsim <- 100
N <- 10^5

############################################################
####### Inland stations analysis #######
############################################################

Atildeland <- Ahatland <- Yoland <- Yotildeland <- vector('list', length = nsim+1)
Atildeland[[1]] <- tailDepMatrix(dataland, alpha = 2, qu = 0.95, mest = FALSE)$A
Ahatland[[1]] <- decompOne(Sigmaland)$A
Ystarland <- MaxLinearSim(Ahatland[[1]], N = N)
Ystartildeland <- MaxLinearSim(Atildeland[[1]], N = N)
Yoland[[1]] <- sapply(c(1:d2), function(i) FgpdInv(dataoland[,i], exp(-Ystarland[,i]^(-2))))
Yotildeland[[1]] <- sapply(c(1:d2), function(i) FgpdInv(dataoland[,i], exp(-Ystartildeland[,i]^(-2))))


for(I in 1:nsim){ #takes a couple of minutes
  print(I)
  set.seed(I)
  nsample <- sample(1:n,size=n,replace=T)
  newdatao <- dataoland[nsample,]
  newdata <- apply(newdatao, 2, function(i) (-log(Fgpd(i, qu = 0.95)))^(-1/2))
  tdm <- tailDepMatrix(newdata, alpha = 2, qu = 0.95, mest = FALSE)
  Ahatland[[I+1]] <- decompOne(tdm$Sigma)$A
  Atildeland[[I+1]] <- tdm$A
  Ystarland <- MaxLinearSim(Ahatland[[I+1]], N = N)
  Ystartildeland <- MaxLinearSim(Atildeland[[I+1]], N = N)
  Yoland[[I+1]] <- sapply(c(1:d2), function(i) FgpdInv(newdatao[,i], exp(-Ystarland[,i]^(-2))))
  Yotildeland[[I+1]] <- sapply(c(1:d2), function(i) FgpdInv(newdatao[,i], exp(-Ystartildeland[,i]^(-2))))
  }

### max-region
xq1land <- apply(dataoland, 2, function(i) (-log(Fgpd(i, y = x1, qu = 0.95)))^(-1/2))
xq2land <- apply(dataoland, 2, function(i) (-log(Fgpd(i, y = x2, qu = 0.95)))^(-1/2))
xq3land <- apply(dataoland, 2, function(i) (-log(Fgpd(i, y = x3, qu = 0.95)))^(-1/2))
empq1land <- length(which(apply(dataoland, 1, max) >= x1))/n # 45 points
empq2land <- length(which(apply(dataoland, 1, max) >= x2))/n # 5 points 
empq3land <- length(which(apply(dataoland, 1, max) >= x3))/n # 1 point 
probahatq1land <- sapply(Ahatland, function(i) MaxLinearProba(xq1land, i, 2, "max"))
probahatq2land <- sapply(Ahatland, function(i) MaxLinearProba(xq2land, i, 2, "max"))
probahatq3land <- sapply(Ahatland, function(i) MaxLinearProba(xq3land, i, 2, "max"))
probatilq1land <- sapply(Atildeland, function(i) MaxLinearProba(xq1land, i, 2, "max"))
probatilq2land <- sapply(Atildeland, function(i) MaxLinearProba(xq2land, i, 2, "max"))
probatilq3land <- sapply(Atildeland, function(i) MaxLinearProba(xq3land, i, 2, "max"))

### sum-region
empq1land2 <- length(which(apply(dataoland, 1, mean) >= x1sum))/n # 48 points
empq2land2 <- length(which(apply(dataoland, 1, mean) >= x2sum))/n # 15 points 
empq3land2 <- length(which(apply(dataoland, 1, mean) >= x3sum))/n # 4 points
probahatq1land2 <- sapply(1:(nsim+1), function(i) mean(rowMeans(Yoland[[i]]) >= x1sum))
probahatq2land2 <- sapply(1:(nsim+1), function(i) mean(rowMeans(Yoland[[i]]) >= x2sum))
probahatq3land2 <- sapply(1:(nsim+1), function(i) mean(rowMeans(Yoland[[i]]) >= x3sum))
probatilq1land2 <- sapply(1:(nsim+1), function(i) mean(rowMeans(Yotildeland[[i]]) >= x1sum))
probatilq2land2 <- sapply(1:(nsim+1), function(i) mean(rowMeans(Yotildeland[[i]]) >= x2sum))
probatilq3land2 <- sapply(1:(nsim+1), function(i) mean(rowMeans(Yotildeland[[i]]) >= x3sum))


namesbox <- c(expression(paste(hat(A), " (exact)")), expression(tilde(A)))

pdf("WINDmaxlandq1.pdf")
par(cex.lab=2,cex.axis=2,cex.main=1.5,mar=c(5,5.5,4,2),mgp=c(3, 2, 0))
boxplot(cbind(probahatq1land,probatilq1land), xlab = "", ylim = c(0.008,0.014),
        main = "Inland, max, 100 km/h", names = namesbox)
points(rep(empq1land,2), pch = 15, cex = 1.75, col = "red")
dev.off()

pdf("WINDmaxlandq2.pdf")
par(cex.lab=2,cex.axis=2,cex.main=1.5,mar=c(5,5.5,4,2),mgp=c(3, 2, 0))
boxplot(cbind(probahatq2land,probatilq2land), xlab = "", ylim = c(0.0007,0.0014),
        main = "Inland, max, 120 km/h", names = namesbox)
points(rep(empq2land,2), pch = 15, cex = 1.75, col = "red")
dev.off()

pdf("WINDmaxlandq3.pdf")
par(cex.lab=2,cex.axis=2,cex.main=1.5,mar=c(5,5.5,4,2),mgp=c(3, 2, 0))
boxplot(cbind(probahatq3land,probatilq3land), xlab = "", ylim = c(0.00019,0.0003),
        main = "Inland, max, 130 km/h", names = namesbox)
points(rep(empq3land,2), pch = 15, cex = 1.75, col = "red")
dev.off()

pdf("WINDsumlandq1.pdf")
par(cex.lab=2,cex.axis=2,cex.main=1.5,mar=c(5,5.5,4,2),mgp=c(3, 2, 0))
boxplot(cbind(probahatq1land2,probatilq1land2), xlab = "", ylim = c(0.008,0.02),
        main = "Inland, sum, 80 km/h", names = namesbox)
points(rep(empq1land2,2), pch = 15, cex = 1.75, col = "red")
dev.off()

pdf("WINDsumlandq2.pdf")
par(cex.lab=2,cex.axis=2,cex.main=1.5,mar=c(5,5.5,4,2),mgp=c(3, 2, 0))
boxplot(cbind(probahatq2land2,probatilq2land2), xlab = "", ylim = c(0.002,0.009),
        main = "Inland, sum, 90 km/h", names = namesbox)
points(rep(empq2land2,2), pch = 15, cex = 1.75, col = "red")
dev.off()

pdf("WINDsumlandq3.pdf")
par(cex.lab=2,cex.axis=2,cex.main=1.5,mar=c(5,5.5,4,2),mgp=c(3, 2, 0))
boxplot(cbind(probahatq3land2,probatilq3land2), xlab = "", ylim = c(0,0.0035),
        main = "Inland, sum, 100 km/h", names = namesbox)
points(rep(empq3land2,2), pch = 15, cex = 1.75, col = "red")
dev.off()

############################################################
####### Coastal stations analysis #######
############################################################

Atildesea <- Ahatsea <- Yosea <- Yotildesea <- vector('list', length = nsim+1)
Atildesea[[1]] <- tailDepMatrix(datasea, alpha = 2, qu = 0.95, mest = FALSE)$A
Ahatsea[[1]] <- decompOne(Sigmasea)$A
Ystarsea <- MaxLinearSim(Ahatsea[[1]], N = N)
Ystartildesea <- MaxLinearSim(Atildesea[[1]], N = N)
Yosea[[1]] <- sapply(c(1:d1), function(i) FgpdInv(dataosea[,i], exp(-Ystarsea[,i]^(-2))))
Yotildesea[[1]] <- sapply(c(1:d1), function(i) FgpdInv(dataosea[,i], exp(-Ystartildesea[,i]^(-2))))


for(I in 1:nsim){ #takes a couple of minutes
  print(I)
  set.seed(I)
  nsample <- sample(1:n,size=n,replace=T)
  newdatao <- dataosea[nsample,]
  newdata <- apply(newdatao, 2, function(i) (-log(Fgpd(i, qu = 0.95)))^(-1/2))
  tdm <- tailDepMatrix(newdata, alpha = 2, qu = 0.95, mest = FALSE)
  Ahatsea[[I+1]] <- decompOne(tdm$Sigma)$A
  Atildesea[[I+1]] <- tdm$A
  Ystarsea <- MaxLinearSim(Ahatsea[[I+1]], N = N)
  Ystartildesea <- MaxLinearSim(Atildesea[[I+1]], N = N)
  Yosea[[I+1]] <- sapply(c(1:d1), function(i) FgpdInv(newdatao[,i], exp(-Ystarsea[,i]^(-2))))
  Yotildesea[[I+1]] <- sapply(c(1:d1), function(i) FgpdInv(newdatao[,i], exp(-Ystartildesea[,i]^(-2))))
}

### max-region
xq1sea <- apply(dataosea, 2, function(i) (-log(Fgpd(i, y = x1, qu = 0.95)))^(-1/2))
xq2sea <- apply(dataosea, 2, function(i) (-log(Fgpd(i, y = x2, qu = 0.95)))^(-1/2))
xq3sea <- apply(dataosea, 2, function(i) (-log(Fgpd(i, y = x3, qu = 0.95)))^(-1/2))
empq1sea <- length(which(apply(dataosea, 1, max) >= x1))/n # 147 points
empq2sea <- length(which(apply(dataosea, 1, max) >= x2))/n # 25 points 
empq3sea <- length(which(apply(dataosea, 1, max) >= x3))/n # 12 points 
probahatq1sea <- sapply(Ahatsea, function(i) MaxLinearProba(xq1sea, i, 2, "max"))
probahatq2sea <- sapply(Ahatsea, function(i) MaxLinearProba(xq2sea, i, 2, "max"))
probahatq3sea <- sapply(Ahatsea, function(i) MaxLinearProba(xq3sea, i, 2, "max"))
probatilq1sea <- sapply(Atildesea, function(i) MaxLinearProba(xq1sea, i, 2, "max"))
probatilq2sea <- sapply(Atildesea, function(i) MaxLinearProba(xq2sea, i, 2, "max"))
probatilq3sea <- sapply(Atildesea, function(i) MaxLinearProba(xq3sea, i, 2, "max"))

### sum-region
empq1sea2 <- length(which(apply(dataosea, 1, mean) >= x1sum))/n # 148 points
empq2sea2 <- length(which(apply(dataosea, 1, mean) >= x2sum))/n # 55 points 
empq3sea2 <- length(which(apply(dataosea, 1, mean) >= x3sum))/n # 17 points
probahatq1sea2 <- sapply(1:(nsim+1), function(i) mean(rowMeans(Yosea[[i]]) >= x1sum))
probahatq2sea2 <- sapply(1:(nsim+1), function(i) mean(rowMeans(Yosea[[i]]) >= x2sum))
probahatq3sea2 <- sapply(1:(nsim+1), function(i) mean(rowMeans(Yosea[[i]]) >= x3sum))
probatilq1sea2 <- sapply(1:(nsim+1), function(i) mean(rowMeans(Yotildesea[[i]]) >= x1sum))
probatilq2sea2 <- sapply(1:(nsim+1), function(i) mean(rowMeans(Yotildesea[[i]]) >= x2sum))
probatilq3sea2 <- sapply(1:(nsim+1), function(i) mean(rowMeans(Yotildesea[[i]]) >= x3sum))


namesbox <- c(expression(paste(hat(A), " (exact)")), expression(tilde(A)))

pdf("WINDmaxseaq1.pdf")
par(cex.lab=2,cex.axis=2,cex.main=1.5,mar=c(5,5.5,4,2),mgp=c(3, 2, 0))
boxplot(cbind(probahatq1sea,probatilq1sea), xlab = "", ylim = c(0.032,0.054),
        main = "Coastal, max, 100 km/h", names = namesbox)
points(rep(empq1sea,2), pch = 15, cex = 1.75, col = "red")
dev.off()

pdf("WINDmaxseaq2.pdf")
par(cex.lab=2,cex.axis=2,cex.main=1.5,mar=c(5,5.5,4,2),mgp=c(3, 2, 0))
boxplot(cbind(probahatq2sea,probatilq2sea), xlab = "", ylim = c(0.0055,0.009),
        main = "Coastal, max, 120 km/h", names = namesbox)
points(rep(empq2sea,2), pch = 15, cex = 1.75, col = "red")
dev.off()

pdf("WINDmaxseaq3.pdf")
par(cex.lab=2,cex.axis=2,cex.main=1.5,mar=c(5,5.5,4,2),mgp=c(3, 2, 0))
boxplot(cbind(probahatq3sea,probatilq3sea), xlab = "", ylim = c(0.002,0.0036),
        main = "Coastal, max, 130 km/h", names = namesbox)
points(rep(empq3sea,2), pch = 15, cex = 1.75, col = "red")
dev.off()

pdf("WINDsumseaq1.pdf")
par(cex.lab=2,cex.axis=2,cex.main=1.5,mar=c(5,5.5,4,2),mgp=c(3, 2, 0))
boxplot(cbind(probahatq1sea2,probatilq1sea2), xlab = "", ylim = c(0.025,0.048),
        main = "Coastal, sum, 80 km/h", names = namesbox)
points(rep(empq1sea2,2), pch = 15, cex = 1.75, col = "red")
dev.off()

pdf("WINDsumseaq2.pdf")
par(cex.lab=2,cex.axis=2,cex.main=1.5,mar=c(5,5.5,4,2),mgp=c(3, 2, 0))
boxplot(cbind(probahatq2sea2,probatilq2sea2), xlab = "", ylim = c(0.01,0.02),
        main = "Coastal, sum, 90 km/h", names = namesbox)
points(rep(empq2sea2,2), pch = 15, cex = 1.75, col = "red")
dev.off()

pdf("WINDsumseaq3.pdf")
par(cex.lab=2,cex.axis=2,cex.main=1.5,mar=c(5,5.5,4,2),mgp=c(3, 2, 0))
boxplot(cbind(probahatq3sea2,probatilq3sea2), xlab = "", ylim = c(0.003,0.009),
        main = "Coastal, sum, 100 km/h", names = namesbox)
points(rep(empq3sea2,2), pch = 15, cex = 1.75, col = "red")
dev.off()


