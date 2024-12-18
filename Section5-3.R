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
data <- apply(dataorig, 2, function(i) (-log(Fgpd(i, qu = 0.95)))^(-1/2))
dataland <- data[,(d1+1):d]
datasea <- data[,1:d1]

# thresholds
x1 <- 1
x2 <- 1.2
x3 <- 1.3


########################################################
####### Chi and TPRM #######
############################################################

Sigma <- tailDepMatrix(data, alpha = 2, qu = 0.95, mest = FALSE)$Sigma
Sigmasea <- tailDepMatrix(datasea, alpha = 2, qu = 0.95, mest = FALSE)$Sigma
Sigmaland <- tailDepMatrix(dataland, alpha = 2, qu = 0.95, mest = FALSE)$Sigma

chis <- chiPairs(data, coord$LAT, coord$LON)
chisea <- chiPairs(data[,1:d1], coord$LAT[1:d1], coord$LON[1:d1])
chiland <- chiPairs(data[,(d1+1):d], coord$LAT[(d1+1):d], coord$LON[(d1+1):d])

par(cex.lab=2,cex.axis=2,cex.main=2,mar=c(5,5.5,4,2))
plot(chis, pch = 19, xlab = 'distance (in km)', ylab = expression(paste(chi[jk])), main = expression(paste("Tail dependence coefficients ", chi[jk])), ylim = c(0.4, 0.9))
points(chisea, col = "brown", pch = 19)
points(chiland, col = "Forestgreen", pch = 19)
tmpsea <- loess.smooth(chisea[,1], chisea[,2])
tmpland <- loess.smooth(chiland[,1], chiland[,2])
lines(tmpsea, lwd = 3, col = "brown")
lines(tmpland, lwd = 3, col = "Forestgreen")
legend("topright", legend = c("coast", "inland"), col = c("brown", "Forestgreen"), lwd = c(3,3), cex = 2)

tpdms <- cbind(chis[,1], Sigma[lower.tri(Sigma)])
tpdmsea <- cbind(chisea[,1], Sigmasea[lower.tri(Sigmasea)])
tpdmland <- cbind(chiland[,1], Sigmaland[lower.tri(Sigmaland)])

par(cex.lab=2,cex.axis=2,cex.main=2,mar=c(5,5.5,4,2))
plot(tpdms, pch = 19, xlab = 'distance (in km)', ylab = expression(paste(sigma[jk])), main = expression(paste("TPDM entries ", sigma[jk])), ylim = c(0.7, 1.05))
points(tpdmsea, col = "brown", pch = 19)
points(tpdmland, col = "ForestGreen", pch = 19)
tmpsea <- loess.smooth(tpdmsea[,1], tpdmsea[,2])
tmpland <- loess.smooth(tpdmland[,1], tpdmland[,2])
lines(tmpsea, lwd = 3, col = "brown")
lines(tmpland, lwd = 3, col = "Forestgreen")
legend("topright", legend = c("coast", "land"), col = c("brown", "Forestgreen"), lwd = c(3,3), cex = 2)


############################################################
####### Inland stations analysis #######
############################################################

nsim <- 100
Atildeland <- Ahatland <- vector('list', length = nsim+1)
Atildeland[[1]] <- tailDepMatrix(dataland, alpha = 2, qu = 0.95, mest = FALSE)$A
Ahatland[[1]] <- decompOne(Sigmaland)$A

for(I in 1:nsim){ #takes a minute
  set.seed(I)
  nsample <- sample(1:n,size=n,replace=T)
  newdatao <- dataoland[nsample,]
  newdata <- apply(newdatao, 2, function(i) (-log(Fgpd(i, qu = 0.95)))^(-1/2))
  tdm <- tailDepMatrix(newdata, alpha = 2, qu = 0.95, mest = FALSE)
  Ahatland[[I+1]] <- decompOne(tdm$Sigma)$A
  Atildeland[[I+1]] <- tdm$A
}


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

namesbox <- c(expression(paste(hat(A), " (exact)")), expression(tilde(A)))

par(cex.lab=2,cex.axis=2,cex.main=1.5,mar=c(5,5.5,4,2),mgp=c(3, 2, 0))
boxplot(cbind(probahatq1land,probatilq1land), xlab = "", ylim = c(0.009,0.0155),
        main = "Inland, 100 km/h threshold", names = namesbox)
points(rep(empq1land,2), pch = 15, cex = 1.75, col = "red")

par(cex.lab=2,cex.axis=2,cex.main=1.5,mar=c(5,5.5,4,2),mgp=c(3, 2, 0))
boxplot(cbind(probahatq2land,probatilq2land), xlab = "", ylim = c(0.0007,0.0014),
        main = "Inland, 120 km/h threshold", names = namesbox)
points(rep(empq2land,2), pch = 15, cex = 1.75, col = "red")

par(cex.lab=2,cex.axis=2,cex.main=1.5,mar=c(5,5.5,4,2),mgp=c(3, 2, 0))
boxplot(cbind(probahatq3land,probatilq3land), xlab = "", ylim = c(0.00019,0.00032),
        main = "Inland, 144 km/h threshold", names = namesbox)
points(rep(empq3land,2), pch = 15, cex = 1.75, col = "red")

############################################################
####### Coastal stations analysis #######
############################################################

nsim <- 100
Atildesea <- Ahatsea <- vector('list', length = nsim+1)
Atildesea[[1]] <- tailDepMatrix(datasea, alpha = 2, qu = 0.95, mest = FALSE)$A
Ahatsea[[1]] <- decompOne(Sigmasea)$A

for(I in 1:nsim){ #takes a minute
  set.seed(I)
  nsample <- sample(1:n,size=n,replace=T)
  newdatao <- dataosea[nsample,]
  newdata <- apply(newdatao, 2, function(i) (-log(Fgpd(i, qu = 0.95)))^(-1/2)) 
  tdm <- tailDepMatrix(newdata, alpha = 2, qu = 0.95, mest = FALSE)
  Ahatsea[[I+1]] <- decompOne(tdm$Sigma)$A
  Atildesea[[I+1]] <- tdm$A
}


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

namesbox <- c(expression(paste(hat(A), " (exact)")), expression(tilde(A)))

par(cex.lab=2,cex.axis=2,cex.main=1.5,mar=c(5,5.5,4,2),mgp=c(3, 2, 0))
boxplot(cbind(probahatq1sea,probatilq1sea), xlab = "", ylim = c(0.035,0.06),
        main = "Coastal, 100 km/h threshold", names = namesbox)
points(rep(empq1sea,2), pch = 15, cex = 1.75, col = "red")

par(cex.lab=2,cex.axis=2,cex.main=1.5,mar=c(5,5.5,4,2),mgp=c(3, 2, 0))
boxplot(cbind(probahatq2sea,probatilq2sea), xlab = "", ylim = c(0.006,0.0105),
        main = "Coastal, 120 km/h threshold", names = namesbox)
points(rep(empq2sea,2), pch = 15, cex = 1.75, col = "red")

par(cex.lab=2,cex.axis=2,cex.main=1.5,mar=c(5,5.5,4,2),mgp=c(3, 2, 0))
boxplot(cbind(probahatq3sea,probatilq3sea), xlab = "", ylim = c(0.0024,0.0042),
        main = "Coastal, 144 km/h threshold", names = namesbox)
points(rep(empq3sea,2), pch = 15, cex = 1.75, col = "red")



