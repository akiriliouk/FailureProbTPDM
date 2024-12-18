library(arrangements)
library(boot)
library(maps)
library(mev)

Fgpd <- function(data, y = NULL, qu = 0.95){ 
  n <- length(data)
  ud <- quantile(data,qu)
  indx <- which(data > ud)
  nindx <- c(1:n)[-indx]
  par <- fit.gpd(data, ud)$estimate
  if(is.null(y)){
    res <- rep(0,n)
    ecdfY <- ecdf(data)
    res[indx] <- 1 - (1-qu)*pmax(1 + par[2]*(pmax(data[indx] - ud,0))/par[1], 0)^(-1/par[2])
    
    res[nindx] <- ecdfY(data[nindx])
  } else{
    if(y <= ud){
      res <- length(which(data <= y))/n
    } else{
      res <- 1 - (1-qu)*pmax(1 + par[2]*(y - ud)/par[1],0)^(-1/par[2])
    }
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


