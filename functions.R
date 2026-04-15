library(arrangements)
library(boot)
library(maps)
library(mev)

Fgpd <- function(data, y = NULL, qu = 0.95){ 
  n <- length(data)
  ud <- quantile(data,qu)
  quu <- mean(data <= ud)       
  indx <- which(data > ud)
  nindx <- c(1:n)[-indx]
  par <- fit.gpd(data, ud)$estimate
  if(is.null(y)){
    res <- rep(0,n)
    ecdfY <- ecdf(data)
    res[indx] <- 1 - (1-quu)*pmax(1 + par[2]*(pmax(data[indx] - ud,0))/par[1], 0)^(-1/par[2])
    
    res[nindx] <- ecdfY(data[nindx])
  } else{
    if(y <= ud){
      res <- length(which(data <= y))/n
    } else{
      res <- 1 - (1-quu)*pmax(1 + par[2]*(y - ud)/par[1],0)^(-1/par[2])
    }
  }
  return(res)
}

FgpdInv <- function(data, p, qu = 0.95) {
  u  <- as.numeric(quantile(data, qu))
  quu <- mean(data <= u)
  par <- fit.gpd(data, u)$estimate
  res <- rep(0, length(p))
  
  idx_lo <- (p <= quu)
  if (any(idx_lo)) {
    xs <- sort(data)
    k  <- ceiling(n * p[idx_lo] - 1e-12)  # subtract tiny eps to avoid k->k+1
    k  <- pmax(1, pmin(n, k))
    res[idx_lo] <- xs[k]
  }
  
  idx_hi <- (p > quu)
  if (any(idx_hi)) {
    ph <- p[idx_hi]
    res[idx_hi] <- u + (par[1]/par[2]) * (((1 - ph)/(1 - quu))^(-par[2]) - 1)
  }
  
  return(res)
}  

Hill <- function(datas, k) sum(log(datas[1:k]/datas[k+1]))/k

Hilleye <- function(data, knseq, smooth = FALSE, const = 0.01, eps = 0.3, h = 0.9, u = 2){ 
  n <- length(data)
  datas <- sort(data, decreasing = TRUE)
  hill <- sapply(knseq, function(k) Hill(datas, k))
  if(smooth){
    sequ <- min(knseq):floor(max(knseq/u))
    avhill <- sapply(sequ, function(k) mean(hill[k:(u * k)]))
    hill <- avhill
    knseq <- sequ
  }  
  alpha <- 1/hill
  w <- floor(const*n)
  smax <- length(knseq) - w
  skeye <- sapply(c(1:smax), function(k) mean(abs(alpha[k:(k+w)]-alpha[k]) <= eps))
  keye <- min(knseq) - 1 + which(skeye > h)[1]
  return(list('alpha' = alpha[keye], 'keye' = keye))
}

HillBoot <- function(data, kn){
  datas <-  sort(data, decreasing = TRUE)
  hill <- Hill(datas, kn)
  return(1/hill)
}

tailDepMatrix <- function(x, alpha, qu, mest = TRUE){ 
  n <- nrow(x)
  R <- apply(x, 1, function(i) sum(i^alpha)^(1/alpha))
  W <- t(sapply(c(1:n), function(i) x[i,]/R[i]))
  r <- quantile(R, qu)
  indx <- which(R > r)
  nexc <- length(indx)
  if(mest){
    m <- nexc*((r^alpha)/n)
  } else{
    m <- ncol(x)
  }
  Wvec <- (W[indx,])^(alpha/2)
  Sigma <- (m/nexc)*(t(Wvec) %*% Wvec)
  return(list('Sigma' = Sigma, 'A' = t(sqrt(m/nexc)*Wvec), 'm' = m))
}


chiEmp <-function(data, u){
  n <- nrow(data)
  datatr <- (apply(data,2,rank) - 0.5)/n   
  cu <- mean(apply(datatr, 1, min) >= u)
  return(cu/(1-u))
}

chiPairs <- function(data, lat, lon, u = 0.95){
  d <- ncol(data)
  latr <- lat/(180/pi)
  lonr <- lon/(180/pi)
  pairs <- combn(d, 2)
  chi <- tpdm <- dist <- matrix(1, ncol = d, nrow = d)
  for(i in 1:(d-1)){
    for(j in (i+1):d){
      chi[i,j] <- chi[j,i] <- chiEmp(data[,c(i,j)], u)
      dtemp <- 1.609344*3963*acos((sin(latr[i])*sin(latr[j])) 
                                  + cos(latr[i])*cos(latr[j])*cos(lonr[j]-lonr[i]))
      dist[i,j] <- dist[j,i] <- dtemp
    }
  }
  dist <- dist[lower.tri(dist)]
  chi <- chi[lower.tri(chi)]
  return(cbind(dist,chi))
}

# decompose based on a fixed path given by "order"
decomp <- function(Sigma, order = c(1:ncol(Sigma)), tolr = 1e-12){ 
  Sigma0 <- Sigma
  pstart <- ncol(Sigma)
  iall <- c(1:pstart)
  A <- matrix(0, nrow=pstart,ncol=pstart)
  for(p in pstart:2){
    iuser <- order[pstart - p + 1]
    iR <- which(iall == iuser)
    pairs <- combinations(c(1:p)[-iR], 2, replace = TRUE)
    tmp <- apply(pairs, 1, function(j) (Sigma[j[1],iR]*Sigma[j[2],iR])/Sigma[j[1],j[2]])/Sigma[iR,iR]
    D <- max(tmp, na.rm = TRUE)
    if(is.infinite(D)){
      A <- A[,1:(pstart-p), drop = FALSE] 
      break  
    } 
    tau <- Sigma[iR,]/sqrt(Sigma[iR,iR]*max(D,1))
    tau[iR] <- tau[iR]*max(D,1)
    if(p == pstart){
      A[,1] <- tau
    } else{
      A[-order[1:(pstart-p)],(pstart-p+1)] <- tau
    }
    Sigma <- pmax(Sigma[-iR,-iR] - tau[-iR] %*% t(tau[-iR]),0)
    iall <- iall[-iR]
    diff <- sum((Sigma0 - A[,1:(pstart-p+1)] %*% t(A[,1:(pstart-p+1)]))^2)
    if(diff < tolr){ #Check if < d columns suffice
      A <- A[,1:(pstart-p+1), drop = FALSE]        
      break
    }
  }
  if(ncol(A) == pstart){
    A[order[pstart],pstart] <- sqrt(Sigma)
  }
  if(sum(A[,ncol(A)]) < 1e-06){A <- A[,-ncol(A)]}
  return(A)
}

# decompose choosing the minimal value of D
decompmin <- function(Sigma, shuffle = FALSE, seed = 1, tolr = 1e-12){
  set.seed(seed)
  Sigma0 <- Sigma
  iR <- order <- c()
  pstart <- ncol(Sigma)
  iall <- c(1:pstart)
  A <- matrix(0, nrow=pstart,ncol=pstart)
  for(p in pstart:2){ 
    D <- rep(0, p)
    for(i in 1:p){
      pairs <- combinations(c(1:p)[-i], 2, replace = TRUE)
      D[i] <- max(apply(pairs, 1, function(j) (Sigma[j[1],i]*Sigma[j[2],i])/Sigma[j[1],j[2]])/Sigma[i,i])
    }  
    iset <- which(D < 1)
    if(shuffle && length(iset) > 1){
      imin <- iset[sample(c(1:length(iset)),1)]
    } else{
      imin <- which(D == min(D))[1]
    }
    if(is.infinite(D[imin])){
      A <- A[,1:(pstart-p), drop = FALSE] 
      break  
    } 
    tau <- Sigma[imin,]/sqrt(Sigma[imin,imin]*max(D[imin],1))
    tau[imin] <- tau[imin]*max(D[imin],1)
    
    iR <- append(iR,imin)
    if(p == pstart){
      order <- append(order, imin)
      A[,1] <- tau
    } else{
      order <- append(order,iall[iR[(pstart-p+1)]])
      A[-order[1:(pstart-p)],(pstart-p+1)] <- tau
    }
    Sigma <- Sigma[-imin,-imin] - tau[-imin] %*% t(tau[-imin])
    Sigma <- pmax(Sigma,0)
    iall <- iall[-imin]
    
    diff <- sum((Sigma0 - A[,1:(pstart-p+1)] %*% t(A[,1:(pstart-p+1)]))^2)
    if(diff < tolr){ #Check if < d columns suffice
      A <- A[,1:(pstart-p+1), drop = FALSE]        
      break
    }
  }
  order <- append(order, c(1:pstart)[-order[1:(pstart-1)]])
  A[order[pstart],pstart] <- sqrt(Sigma)
  return(list('A' = A, 'path' = order))
}

# Try arbitrary paths until you find an exact decomposition
decompOne <- function(Sigma, maxsim = 10^5, tolr = 1e-12, startseed = 0){
  count <- 0 
  flag <- TRUE
  while(count < maxsim & flag){
    count <- count + 1
    temp <- decompmin(Sigma, shuffle = TRUE, seed = count + startseed)
    if(sqrt(sum(abs(diag(Sigma) - diag(temp$A %*% t(temp$A)))^2)) < tolr){
      flag <- FALSE
    }
  }
  return(list('A' = temp$A, 'path' = temp$path))
}

MaxLinearProba <- function(x, A, alpha, region = c("min", "max", "sum"), v = NULL, prob = 0){
  if(region == "min"){
    return(sum(apply(A,2,function(i) min(i/x)^alpha)) - prob)
  } else if(region == "max"){
    return(sum(apply(A,2,function(i) max(i/x)^alpha)) - prob)
  } else if(region == "sum"){
    return((x^(-alpha))*sum(apply(A, 2, function(i) sum(v*i)^alpha)) - prob)
  }
}

MaxLinearSim <- function(A, N, seed = 1, alpha = 2){
  d <- nrow(A)
  q <- ncol(A)
  set.seed(seed)
  
  U <- matrix(runif(N*q), nrow = N, ncol = q)
  Z <- (-log(U))^(-1/alpha)
  Y <- matrix(0, nrow = N, ncol = d)
  
  for(l in 1:q){
    tmp <- tcrossprod(Z[, l], A[, l])  
    Y <- pmax(Y, tmp)
  }
  return(Y)
}


######### Chen's functions (adapted)

Hill2 <- function(k,x){
  y = log(sort(x,decreasing = TRUE)[1:(k+1)])
  mean(y[1:k]-y[k+1])
}

stat_single = function(kj,x,index, threshold,gamma,omega){ #gamma = average of component-wise hill estimators
  true_index = which(index==1)
  sums = 0 
  for(i in true_index){
    sums = sums + omega[i]*(log(x[i]) - log(threshold) - gamma) 
  }
  res = abs(sums)/gamma/sqrt(kj)
  return(res)
}

test = function(vk,vX,vindex,v_threshold,gamma,omega, p){
  stat_vec = numeric(p)
  for(j in 1:p){
    stat_vec[j] = stat_single(vk[j],vX[,j],vindex[,j],v_threshold[j],gamma,omega)
  } 
  return(max(stat_vec))
}

test_bootstrap = function(vk,vX,B=2000){
  n = nrow(vX)
  p = ncol(vX)
  
  vindex = matrix(0,nrow=n,ncol=p)
  v_threshold = numeric(p)
  for(j in 1:p){
    v_threshold[j] = sort(vX[,j],decreasing = TRUE)[vk[j]+1]
    vindex[,j] = vX[,j]>v_threshold[j]
  }
  
  vgamma_hat = numeric(p)
  for(j in 1:p){
    vgamma_hat[j] = Hill2(vk[j], vX[,j])
  }
  gamma = mean(vgamma_hat)
  
  res_boots = numeric(B)
  for(s in 1:B){
    res_boots[s] = test(vk,vX,vindex,v_threshold,gamma,rnorm(n), p)
  }
  test_orginal = test(vk,vX,vindex,v_threshold,gamma,rep(1,n), p)
  
  pval = mean(test_orginal<=res_boots)
  if(pval <=0.05) print("reject")
  as.numeric(pval)
}
