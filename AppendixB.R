source("functions.R")

############################################################
####### Data pre-processing #######
############################################################

load("Rain.RData")
n <- nrow(dataorig)
d <- ncol(dataorig)


alphaci <- matrix(0, ncol = 3, nrow = d)
for(i in 1:d){
  alpha <- Hilleye(dataorig[,i], smooth = FALSE, knseq = c(10:1000))  
  bootalpha <- tsboot(dataorig[,i], statistic = HillBoot, R = 500, sim = "geom",
                            l = 200, kn = alpha$keye)
  alphaci[i,1] <- bootalpha$t0
  alphaci[i,2] <- quantile(bootalpha$t, 0.025)
  alphaci[i,3] <- quantile(bootalpha$t, 0.925)
}

kall <- Hilleye(c(dataorig), c(10:10000), smooth = FALSE)$keye
(alpha <- 1/Hill(sort(c(dataorig), decreasing = TRUE), kall))


plot(alphaci[,1], pch = 20, col = "red", ylim = c(0,7), ylab = expression(paste(hat(alpha)[i])), 
     xlab = "i", main = expression(paste("Estimates ", hat(alpha)[i], " with 95 % bootstrap confidence intervals")))
points(alphaci[,2], pch = 20, col = "blue")
points(alphaci[,3], pch = 20, col = "blue")
for(i in 1:d){
  segments(i, alphaci[i,2], i, alphaci[i,3], col = 'grey')
}
abline(h = alpha) 


data <- apply(dataorig, 2, function(i) (-log(Fgpd(i, qu = 0.95)))^(-1/2))

#################################################################################
########### Results ############################################################
###############################################################################

#nsim <- 100
#Atilde <- Ahat <- vector('list', length = nsim+1)
#Sigma <- tailDepMatrix(data, alpha = 2, qu = 0.95, mest = FALSE)$Sigma
#Atilde[[1]] <- tailDepMatrix(data, alpha = 2, qu = 0.95, mest = FALSE)$A
#Ahat[[1]] <- decompOne(Sigma, tolr = 5)$A

#for(I in 1:nsim){
#  set.seed(I)
#  nsample <- sample(1:n,size=n,replace=T)
#  newdatao <- dataorig[nsample,]
#  newdata <- apply(newdatao, 2, function(i) (-log(Fgpd(i, qu = 0.95)))^(-1/2))
#  tdm <- tailDepMatrix(newdata, alpha = 2, qu = 0.95, mest = FALSE)
#  Ahat[[I+1]] <- decompOne(tdm$Sigma, tolr = 5)$A
#  Atilde[[I+1]] <- tdm$A
#}
#save(Ahat, Atilde, file = "Rainres.RData")

load("Rainres.RData")

x1 <- apply(dataorig, 2, function(i) (-log(Fgpd(i, y = 80, qu = 0.95)))^(-1/2))
x2 <- apply(dataorig, 2, function(i) (-log(Fgpd(i, y = 110, qu = 0.95)))^(-1/2))
x3 <- apply(dataorig, 2, function(i) (-log(Fgpd(i, y = 140, qu = 0.95)))^(-1/2))
emp1 <- length(which(apply(data, 1, function(i) any(i >= x1))))/n
emp2 <- length(which(apply(data, 1, function(i) any(i >= x2))))/n
emp3 <- length(which(apply(data, 1, function(i) any(i >= x3))))/n

probahat1 <- sapply(Ahat, function(i) MaxLinearProba(x1, i, 2, "max"))
probahat2 <- sapply(Ahat, function(i) MaxLinearProba(x2, i, 2, "max"))
probahat3 <- sapply(Ahat, function(i) MaxLinearProba(x3, i, 2, "max"))

probatil1 <- sapply(Atilde, function(i) MaxLinearProba(x1, i, 2, "max"))
probatil2 <- sapply(Atilde, function(i) MaxLinearProba(x2, i, 2, "max"))
probatil3 <- sapply(Atilde, function(i) MaxLinearProba(x3, i, 2, "max"))

namesbox <- c(expression(paste(hat(A), " (approx)")), expression(tilde(A)))

par(cex.lab=2,cex.axis=2,cex.main=1.5,mar=c(5,5.5,4,2),mgp=c(3, 2, 0))
boxplot(cbind(probahat1,probatil1), xlab = "",
        main = expression(paste(hat(p)[max], ", 80 mm")), names = namesbox)
points(rep(emp1,2), pch = 15, cex = 1.75, col = "red")

par(cex.lab=2,cex.axis=2,cex.main=1.5,mar=c(5,5.5,4,2),mgp=c(3, 2, 0))
boxplot(cbind(probahat2,probatil2), xlab = "",
        main = expression(paste(hat(p)[max], ", 110 mm")), names = namesbox)
points(rep(emp2,2), pch = 15, cex = 1.75, col = "red")

par(cex.lab=2,cex.axis=2,cex.main=1.5,mar=c(5,5.5,4,2),mgp=c(3, 2, 0))
boxplot(cbind(probahat3,probatil3), xlab = "",
        main = expression(paste(hat(p)[max], ", 140 mm")), names = namesbox)
points(rep(emp3,2), pch = 15, cex = 1.75, col = "red")


x1 <- apply(dataorig, 2, function(i) (-log(Fgpd(i, y = 5, qu = 0.95)))^(-1/2))
x2 <- apply(dataorig, 2, function(i) (-log(Fgpd(i, y = 10, qu = 0.95)))^(-1/2))
x3 <- apply(dataorig, 2, function(i) (-log(Fgpd(i, y = 30, qu = 0.95)))^(-1/2))
emp1 <- length(which(apply(data, 1, function(i) all(i >= x1))))/n
emp2 <- length(which(apply(data, 1, function(i) all(i >= x2))))/n
emp3 <- length(which(apply(data, 1, function(i) all(i >= x3))))/n

probahat1 <- sapply(Ahat, function(i) MaxLinearProba(x1, i, 2, "min"))
probahat2 <- sapply(Ahat, function(i) MaxLinearProba(x2, i, 2, "min"))
probahat3 <- sapply(Ahat, function(i) MaxLinearProba(x3, i, 2, "min"))

probatil1 <- sapply(Atilde, function(i) MaxLinearProba(x1, i, 2, "min"))
probatil2 <- sapply(Atilde, function(i) MaxLinearProba(x2, i, 2, "min"))
probatil3 <- sapply(Atilde, function(i) MaxLinearProba(x3, i, 2, "min"))

namesbox <- c(expression(paste(hat(A), " (approx)")), expression(tilde(A)))

par(cex.lab=2,cex.axis=2,cex.main=1.5,mar=c(5,5.5,4,2),mgp=c(3, 2, 0))
boxplot(cbind(probahat1,probatil1), xlab = "", ylim = c(0,0.12),
        main = expression(paste(hat(p)[min], ", 10 mm")), names = namesbox)
points(rep(emp1,2), pch = 15, cex = 1.75, col = "red")

par(cex.lab=2,cex.axis=2,cex.main=1.5,mar=c(5,5.5,4,2),mgp=c(3, 2, 0))
boxplot(cbind(probahat2,probatil2), xlab = "", ylim = c(0,0.06),
        main = expression(paste(hat(p)[min], ", 20 mm")), names = namesbox)
points(rep(emp2,2), pch = 15, cex = 1.75, col = "red")

par(cex.lab=2,cex.axis=2,cex.main=1.5,mar=c(5,5.5,4,2),mgp=c(3, 2, 0))
boxplot(cbind(probahat3,probatil3), xlab = "", ylim = c(0,0.006),
        main = expression(paste(hat(p)[min], ", 30 mm")), names = namesbox)
points(rep(emp3,2), pch = 15, cex = 1.75, col = "red")



