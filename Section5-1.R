source("functions.R")

A1 <- rbind(c(1.00, 0.50, 0.00, 0.25, 1.75, 0.50, 0.75, 1.00, 1.00,  0.25,  1.75,  0.25,  1.50,  0.50,  0.25,  1.50),
            c(2.00, 0.00, 1.50, 1.00, 1.00, 0.25, 1.00, 1.00, 1.75,  0.25,  0.00,  2.00,  0.25,  0.50,  0.25,  0.75),
            c(1.75, 1.25, 0.75, 0.25, 0.50, 2.00, 1.75, 0.25, 0.75,  0.25,  1.25,  0.25,  1.00,  0.75,  1.00, 0.25),
            c(1.25, 0.25, 2.00, 0.25, 1.25, 2.00, 0.50, 0.25, 0.50,  0.50,  0.00,  0.75,  0.50,  0.25,  1.75,  0.50),
            c(1.75, 0.50, 0.75, 1.25, 0.25, 0.50, 1.75, 0.00, 2.00,  1.00,  1.50,  0.50,  0.00,  0.50,  1.25,  1.25))
A2 <- A1[,1:8]
A3 <- A1[,1:4]
SigmaA1 <- A1 %*% t(A1)
SigmaA2 <- A2 %*% t(A2)
SigmaA3 <- A3 %*% t(A3)

f1 <- function(A, alpha, cst){ #sum
  d <- nrow(A)
  b <- rep(1/d,d)
  temp <- apply(A, 2, function(i) (sum(b*i)^alpha))
  return(sum(temp)*(cst^(-alpha)))
}
f2 <- function(A, alpha, cst){ #product
  d <- nrow(A)
  temp <- apply(A, 2, function(i) prod(i^(alpha/d)))
  return(sum(temp)*(cst^(-alpha/d)))
}
f3 <- function(A, alpha, cst = 5){ #min
  temp <- apply(A, 2, function(i) (min(i/cst)^alpha))
  return(sum(temp))
}
f4 <- function(A, alpha, cst = 5){ #max
  temp <- apply(A, 2, function(i) (max(i/cst)^alpha))
  return(sum(temp))
}

simResults <- function(alpha, Sigma, A){
  d <- nrow(Sigma)
  perms <- permutations(d)
  rep <- factorial(d)
  Alist <- vector('list', length = rep)
  diff <- diffdiag <- vector(length = rep)
  nu <- matrix(0, ncol = 4, nrow = rep)
  for(i in 1:rep){
    Alist[[i]] <- decomp(Sigma, order = perms[i,])
    temp <- Alist[[i]] %*% t(Alist[[i]])
    diff[i] <- sqrt(sum((temp[lower.tri(temp)]-Sigma[lower.tri(Sigma)])^2))
    print(c(i,diff[i]))
    diffdiag[i] <- sqrt(sum((diag(Sigma) - diag(temp))^2))
    nu[i,] <- c(f1(Alist[[i]]^(2/alpha), alpha, uniroot(function(x) f1(A^(2/alpha), alpha, x) - 0.1, c(1,10^10))$root),
                 f2(Alist[[i]]^(2/alpha), alpha, uniroot(function(x) f2(A^(2/alpha), alpha,  x) - 0.1, c(1,10^10))$root),
                 f3(Alist[[i]]^(2/alpha),  alpha, uniroot(function(x) f3(A^(2/alpha), alpha, x) - 0.1, c(1,10^10))$root),
                 f4(Alist[[i]]^(2/alpha),  alpha, uniroot(function(x) f4(A^(2/alpha), alpha, x) - 0.1, c(1,10^10))$root))
  }
  indx1 <- (which(diffdiag <= 1e-10 & diff < 1e-10))
  indx2 <- (which(diffdiag <= 5 & diff < 1e-10))
  nmbr <- length(which(diff < 1e-10))
  return(list('nu1' = nu[indx1,], 'nu2' = nu[indx2,], 'nmbr' = nmbr))
}

alpha <- 4
nuA1 <- simResults(alpha, SigmaA1, A1)
nuA2 <- simResults(alpha, SigmaA2, A2)
nuA3 <- simResults(alpha, SigmaA3, A3)
c(nuA1$nmbr,nuA2$nmbr,nuA3$nmbr)
numbers <- c(nrow(nuA1$nu1),nrow(nuA2$nu1),nrow(nuA3$nu1),
             nrow(nuA1$nu2),nrow(nuA2$nu2),nrow(nuA3$nu2))
names <- c(expression(paste(A[1], " (38)")), expression(paste(A[1], " (68)")),
           expression(paste(A[2], " (12)")), expression(paste(A[2], " (58)")),
           expression(paste(A[3], " (16)")), expression(paste(A[3], " (72)")))

par(cex.lab=1.5,cex.axis=1.4,cex.main=1.75,mar=c(5,5,4,2))
boxplot(nuA1$nu1[,1],nuA1$nu2[,1],nuA2$nu1[,1],nuA2$nu2[,1],nuA3$nu1[,1],nuA3$nu2[,1],
        names = names, main = expression(paste("exponent measures of ", C[(sum)])), ylim = c(0.095,0.11))
abline(h = 0.1, lwd  = 2)

par(cex.lab=1.5,cex.axis=1.4,cex.main=1.75,mar=c(5,5,4,2))
boxplot(nuA1$nu1[,2],nuA1$nu2[,2],nuA2$nu1[,2],nuA2$nu2[,2],nuA3$nu1[,2],nuA3$nu2[,2],
        names = names, main = expression(paste("exponent measures of ", C[(prod)])), ylim = c(0.095,0.14))
abline(h = 0.1, lwd  = 2)

par(cex.lab=1.5,cex.axis=1.4,cex.main=1.75,mar=c(5,5,4,2))
boxplot(nuA1$nu1[,3],nuA1$nu2[,3],nuA2$nu1[,3],nuA2$nu2[,3],nuA3$nu1[,3],nuA3$nu2[,3],
        names = names, main = expression(paste("exponent measures of ", C[(min)])), ylim = c(0,0.35))
abline(h = 0.1, lwd  = 2)

par(cex.lab=1.5,cex.axis=1.4,cex.main=1.75,mar=c(5,5,4,2))
boxplot(nuA1$nu1[,4],nuA1$nu2[,4],nuA2$nu1[,4],nuA2$nu2[,4],nuA3$nu1[,4],nuA3$nu2[,4],
        names = names, main = expression(paste("exponent measures of ", C[(max)])), ylim = c(0.075,0.15))
abline(h = 0.1, lwd  = 2)



