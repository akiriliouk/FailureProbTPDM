source("functions.R") #

############################################################
####### Data pre-processing & plotting #######
############################################################

load("FIN.RData")
n <- nrow(dataorig)
d <- ncol(dataorig)
data <- apply(dataorig, 2, function(i) pmax(-i, rep(0, n))) 

alphaci <- matrix(0, ncol = 3, nrow = d)
for(i in 1:d){ #takes a minute
  alpha <- Hilleye(data[,i], smooth = FALSE, knseq = c(15:1000))  
  bootalpha <- tsboot(data[,i], statistic = HillBoot, R = 500, sim = "geom",
                      l = 200, kn = alpha$keye)
  alphaci[i,1] <- bootalpha$t0
  alphaci[i,2] <- quantile(bootalpha$t, 0.025)
  alphaci[i,3] <- quantile(bootalpha$t, 0.925)
}

kall <- Hilleye(c(data), c(50:5000), smooth = FALSE)$keye
(alpha <- 1/Hill(sort(c(data), decreasing = TRUE), kall))

par(cex.lab=2,cex.axis=2,cex.main=1.5,mar=c(5,5.5,4,2))
plot(alphaci[,1], pch = 20, col = "red", ylim = c(0,7), ylab = expression(paste(hat(alpha)[j])), 
     xlab = "j", main = expression(paste("Estimates ", hat(alpha)[j], " with 95 % bootstrap confidence intervals")))
points(alphaci[,2], pch = 20, col = "blue")
points(alphaci[,3], pch = 20, col = "blue")
for(i in 1:d){
  segments(i, alphaci[i,2], i, alphaci[i,3], col = 'grey')
}
abline(h = alpha) 

R <- apply(data, 1, function(i) sum(i^alpha)^(1/alpha))
qu <- seq(0.95,0.999, by = 0.0005)
r <- quantile(R, qu)
m <- vector(length = length(qu))
for(i in 1:length(qu)){
  indx <- which(R > r[i])
  nexc <- length(indx)
  m[i] <- nexc*((r[i]^alpha)/n)
}

par(cex.lab=2,cex.axis=2,cex.main=1.5,mar=c(5,5.5,4,2))
plot(qu, m, ylab = expression(paste(hat(m))), type = "l", lwd = 2,
     xlab = expression(paste("quantile of ", R[1], ", ..., ", R[n])), 
     main = expression(paste("Estimated mass ", hat(m), " as a function of ", r[0])))
abline(h = 30)
dev.off()

#####################################################################
################ Main results ########################
#####################################################################

# The commented code below produces the file SPres.RData. It takes about 24h to run.

#Atilde <- Ahat <- Ahat5 <- vector('list', length = nsim+1)
#alpha <- vector(length = nsim+1)
#kall <- Hilleye(c(data), c(50:10000), smooth = FALSE)$keye
#alpha[1] <- 1/Hill(sort(c(data), decreasing = TRUE), kall)
#Sigma <- tailDepMatrix(data, alpha = alpha[1], qu = 0.975, mest = TRUE)$Sigma
#Atilde[[1]] <- tailDepMatrix(data, alpha = alpha[1], qu = 0.975, mest = TRUE)$A
#Ahat[[1]] <- decompOne(Sigma)$A
#Ahat5[[1]] <- decompOne(Sigma, tolr =  5)$A

#I <- count <- 1
#while(count <= nsim){
#  set.seed(I)
#  nsample <- sample(1:n,size=n,replace=T)
#  newdata <- data[nsample,]
#  kall <- Hilleye(c(newdata), c(50:10000), smooth = FALSE)$keye
#  alphatemp <- 1/Hill(sort(c(newdata), decreasing = TRUE), kall)
#  tdm <- tailDepMatrix(newdata, alpha = alphatemp, qu = 0.975, mest = TRUE)
#  Ahattemp <- decompOne(tdm$Sigma, maxsim = 1000)$A #tries during 5 minutes max 
#  if(!is.null(Ahattemp)){
#    alpha[count+1] <- alphatemp
#    Ahat[[count+1]] <- Ahattemp
#    Ahat5[[count+1]] <- decompOne(tdm$Sigma, tolr = 5)$A
#    Atilde[[count+1]] <- tdm$A
#    count <- count + 1
#   save(Ahat, Ahat5, Atilde, alpha, file = "FINres.RData")
#  }
#  I <- I + 1
#}

load("FINres.RData")

nsim <- 100
q <- 0.005
v <- rep(c(0.02,0.05,0.03), 10)

minsum <- function(x, v = rep(1/length(x), length(x))){
  return(min(sum(x[1:10]*v[1:10]),sum(x[11:20]*v[11:20]),sum(x[21:30]*v[21:30])))
}
maxsum <- function(x, v = rep(1/length(x), length(x))){
  return(max(sum(x[1:10]*v[1:10]),sum(x[11:20]*v[11:20]),sum(x[21:30]*v[21:30])))
}

qsum <- quantile(rowSums(data)/d, 1-q, type = 4) 
ressumAtilde <- sapply(c(1:nsim), function(k) MaxLinearProba(qsum, Atilde[[k]]^(2/alpha[k]), alpha[k], "sum", rep(1/d, d)))
ressumA <- sapply(c(1:nsim), function(k) MaxLinearProba(qsum, Ahat[[k]]^(2/alpha[k]), alpha[k], "sum", rep(1/d,d)))
ressumA5 <- sapply(c(1:nsim), function(k) MaxLinearProba(qsum, Ahat5[[k]]^(2/alpha[k]), alpha[k], "sum", rep(1/d,d)))

qminsum <- quantile(apply(data, 1, function(i) minsum(i)), 1-q, type = 5)
resminsumAtilde <- sapply(c(1:nsim), function(k) sum(apply(Atilde[[k]]^(2/alpha[k]),2,function(i) minsum(i/qminsum)^alpha[k])))
resminsumA <- sapply(c(1:nsim), function(k) sum(apply(Ahat[[k]]^(2/alpha[k]),2,function(i) minsum(i/qminsum)^alpha[k])))
resminsumA5 <- sapply(c(1:nsim), function(k) sum(apply(Ahat5[[k]]^(2/alpha[k]),2,function(i) minsum(i/qminsum)^alpha[k])))

qmaxsum <- quantile(apply(data, 1, function(i) maxsum(i)), 1-q, type = 5)
resmaxsumAtilde <- sapply(c(1:nsim), function(k) sum(apply(Atilde[[k]]^(2/alpha[k]),2,function(i) maxsum(i/qmaxsum)^alpha[k])))
resmaxsumA <- sapply(c(1:nsim), function(k) sum(apply(Ahat[[k]]^(2/alpha[k]),2,function(i) maxsum(i/qmaxsum)^alpha[k])))
resmaxsumA5 <- sapply(c(1:nsim), function(k) sum(apply(Ahat5[[k]]^(2/alpha[k]),2,function(i) maxsum(i/qmaxsum)^alpha[k])))

namesbox <- c(expression(paste(hat(A), " (exact)")), expression(paste(hat(A), " (approx)")), expression(tilde(A)))

par(cex.lab=1.5,cex.axis=2,cex.main=2,mar=c(5,5.5,4,2),mgp=c(3, 2, 0))
boxplot(cbind(ressumA, ressumA5, ressumAtilde), names = namesbox, ylim = c(0.003,0.0071),
        main = expression(paste("Estimates ", hat(p)[sum], " (equal ", v, ")")))
points(rep(q,3), pch = 15, cex = 1.75, col = "red")

par(cex.lab=1.5,cex.axis=2,cex.main=2,mar=c(5,5.5,4,2),mgp=c(3, 2, 0))
boxplot(cbind(resminsumA, resminsumA5, resminsumAtilde), names = namesbox, ylim = c(0.003,0.0071),
        main = expression(paste("Estimates ", hat(p)[minsum], " (equal ", v, ")")))
points(rep(q,3), pch = 15, cex = 1.75, col = "red")

par(cex.lab=1.5,cex.axis=2,cex.main=2,mar=c(5,5.5,4,2),mgp=c(3, 2, 0))
boxplot(cbind(resmaxsumA, resmaxsumA5, resmaxsumAtilde), names = namesbox, ylim = c(0.003,0.0071),
        main = expression(paste("Estimates ", hat(p)[maxsum], " (equal ", v, ")")))
points(rep(q,3), pch = 15, cex = 1.75, col = "red")


qsumv <- quantile(apply(data, 1, function(i) sum(v*i)),1-q, type = 4) 
ressumAtildev <- sapply(c(1:nsim), function(k) MaxLinearProba(qsumv, Atilde[[k]]^(2/alpha[k]), alpha[k], "sum", v))
ressumAv <- sapply(c(1:nsim), function(k) MaxLinearProba(qsumv, Ahat[[k]]^(2/alpha[k]), alpha[k], "sum", v))
ressumA5v <- sapply(c(1:nsim), function(k) MaxLinearProba(qsumv, Ahat5[[k]]^(2/alpha[k]), alpha[k], "sum", v))

qminsumv <- quantile(apply(data, 1, function(i) minsum(i, v)), 1-q, type = 5) 
resminsumAtildev <- sapply(c(1:nsim), function(k) sum(apply(Atilde[[k]]^(2/alpha[k]),2,function(i) minsum(i/qminsumv, v)^alpha[k])))
resminsumAv <- sapply(c(1:nsim), function(k) sum(apply(Ahat[[k]]^(2/alpha[k]),2,function(i) minsum(i/qminsum,v)^alpha[k])))
resminsumA5v <- sapply(c(1:nsim), function(k) sum(apply(Ahat5[[k]]^(2/alpha[k]),2,function(i) minsum(i/qminsum,v)^alpha[k])))

qmaxsumv <- quantile(apply(data, 1, function(i) maxsum(i, v)), 1-q, type = 5) 
resmaxsumAtildev <- sapply(c(1:nsim), function(k) sum(apply(Atilde[[k]]^(2/alpha[k]),2,function(i) maxsum(i/qmaxsumv, v)^alpha[k])))
resmaxsumAv <- sapply(c(1:nsim), function(k) sum(apply(Ahat[[k]]^(2/alpha[k]),2,function(i) maxsum(i/qmaxsum,v)^alpha[k])))
resmaxsumA5v <- sapply(c(1:nsim), function(k) sum(apply(Ahat5[[k]]^(2/alpha[k]),2,function(i) maxsum(i/qmaxsum,v)^alpha[k])))


par(cex.lab=1.5,cex.axis=2,cex.main=2,mar=c(5,5.5,4,2),mgp=c(3, 2, 0))
boxplot(cbind(ressumAv, ressumA5v, ressumAtildev), names = namesbox, ylim = c(0.0027,0.0069),
        main = expression(paste("Estimates ", hat(p)[sum], " (unequal ", v, ")")))
points(rep(q,3), pch = 15, cex = 1.75, col = "red")

par(cex.lab=1.5,cex.axis=2,cex.main=2,mar=c(5,5.5,4,2),mgp=c(3, 2, 0))
boxplot(cbind(resminsumAv, resminsumA5v, resminsumAtildev), names = namesbox, ylim = c(0.0027,0.0069),
        main = expression(paste("Estimates ", hat(p)[minsum], " (unequal ", v, ")")))
points(rep(q,3), pch = 15, cex = 1.75, col = "red")

par(cex.lab=1.5,cex.axis=2,cex.main=2,mar=c(5,5.5,4,2),mgp=c(3, 2, 0))
boxplot(cbind(resmaxsumAv, resmaxsumA5v, resmaxsumAtildev), names = namesbox, ylim = c(0.0027,0.0069),
        main = expression(paste("Estimates ", hat(p)[maxsum], " (unequal ", v, ")")))
points(rep(q,3), pch = 15, cex = 1.75, col = "red")
