## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
choose_variable <- function(k, m){
  q <- m - k
  M <- 2**q
  pp <- numeric(q+1)
  
  for (i in 1:(q+1)){
    pp[i] <- choose(q, i-1)
  }
  
  s <- numeric(M)
  s[1] <- k
  
  for (i in 2:(q+1)){
    temp <- sum(s > 0)
    s[(temp+1): (temp+pp[i])] <- rep((i+k-1), pp[i])
  }
  
  S <- matrix(nrow = M, ncol = m)
  
  for (i in 1:M) {
    S[i, 1:k] <- 1:k
  }
  
  tS <- 1
  
  for (i in 1:q) {
    Temp <- t(combn((k+1):m, i))
    S[(tS+1):(tS+nrow(Temp)),(k+1): (k+i)] <- Temp
    tS <- tS + nrow(Temp)
  }
  
  output <- list(s=s, S=S)
  return(output)
}
choose_variable(2,8)

## -----------------------------------------------------------------------------
library(MASS)
set.seed(100)

TSLS_SAIC <- function(G, T, M, Mfixed){
  C <- 2**(M - Mfixed)
  beta <- 0.1
  Rf2 <- 0.1
  pi <- numeric(M)
  ehat <- matrix(0, T, C-1)
  thetahatc <- numeric(C-1)
  AICc <- numeric(C-1)
  expc <- numeric(C-1)
  thetahat <- numeric(G)
  
  s <- choose_variable(Mfixed,M)$s
  S <- choose_variable(Mfixed,M)$S
  
  sumM <- 0
  
  for (i in 1:M) {
    sumM <- sumM + (1-i/(M+1))**8
  }
  
  cM <- (Rf2/(1-Rf2)/sumM)**0.5
  
  for (i in 1:M) {
    pi[i] <- cM*(1-i/(M+1))**4
  }
  
  for (g in 1:G) {
    Ome <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
    mu <- c(0,0)
    eu <- mvrnorm(T, mu, Ome)
    e <- eu[,1]
    u <- eu[,2]
    veta <- mvrnorm(T, mu, diag(2))
    v <- veta[,1]
    eta <- veta[,2]
    Z <- mvrnorm(T, numeric(M), diag(M))
    X <- v + eta
    Y <- Z %*% pi + u + eta
    y <- beta * Y + e
    batX <- Y
    
    for (c in 2:C) {
      pc <- s[c]
      a <- S[c, 1:pc]
      Zc <- Z[, a]
      PZc <- Zc %*% solve(t(Zc) %*% Zc) %*% t(Zc)
      PZcX <- PZc %*% batX
      thetahatc[c-1] <- solve(t(batX) %*% PZcX) %*% t(PZcX) %*% y
      ehat[, c-1]<- y - PZcX %*% thetahatc[c-1]
      RSSc <- t(ehat[, c-1]) %*% ehat[, c-1]
      qc <- 2 + length(a)
      AICc[c-1] <- 2 * qc + T * log(RSSc/T)
      expc[c-1] <- exp(-AICc[(c-1)]/2)
    }
    
    sumexp <- sum(expc)
    omega <- exp(-AICc/2)/sumexp
    thetahat[g] <- t(omega) %*% thetahatc
  }
  
  MAD <- median(abs(thetahat-0.1))
  MB <- median(thetahat-0.1)
  IDR = sort(thetahat)[round(length(thetahat)*0.9)]-sort(thetahat)[round(length(thetahat)*0.1)]
  
  return(c(MAD,MB,IDR))
}
TSLS_SAIC(500, 100, 10, 2)

