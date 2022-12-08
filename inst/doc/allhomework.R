## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----include=FALSE------------------------------------------------------------
#colorFunc <- "heat.colors"
colorFunc <- "rainbow"

## ----fig.cap = "The Maunga Whou volcano", echo=FALSE--------------------------
image(volcano, col = get(colorFunc)(200))

## ----echo = FALSE, results='asis'---------------------------------------------
library(knitr)
kable(mtcars[1:7, ], caption = "A knitr kable.")

## ----airmiles, echo=FALSE-----------------------------------------------------
plot(airmiles)

## -----------------------------------------------------------------------------
set.seed(12345)
n <- 1000
u <- runif(n)
a <- 2
b <- 2
x <- b/{(1-u)^{1/a}}
hist(x, breaks = "Scott",prob = TRUE)
y <- seq(0, 40, .01)
lines(y, a*b^a/{y^(a+1)}) 

## -----------------------------------------------------------------------------
grBeta <- function(n,a,b){
  j <- k <- 0
  y <- numeric(n)
  while (k < n) {
    u <- runif(1)
    j <- j+1
    x <- runif(1)
    if (x^(a-1)*(1-x)^(b-1) > u){
      k <- k+1
      y[k] <- x
    }
  }
  return(y)
}

## -----------------------------------------------------------------------------
y <- grBeta(1000,3,2)
hist(y, breaks = "Scott", prob = TRUE, ylim = c(0,2))
z <- seq(0,10,0.01)
lines(z,12*z^2*(1-z))

## -----------------------------------------------------------------------------
n <- 1000
r <- 4
beta <- 2
lambda <- rgamma(n, r, beta)
x <- rexp(n, rate = lambda)

## -----------------------------------------------------------------------------
hist(x, breaks = "Scott",prob = TRUE)
y <- seq(0, 10, .01)
lines(y, r*beta^r*(beta+y)^(-r-1))

## -----------------------------------------------------------------------------
#a function to implement the fast sorting algorithm
fast_sorting<-function(x){
  len<-length(x)
  if(len==0||len==1){return(x)
  }else{
    a<-x[1]
    y<-x[-1]
    lower<-y[y<a]
    upper<-y[y>=a]
    return(c(fast_sorting(lower),a,fast_sorting(upper)))}
}

## -----------------------------------------------------------------------------
#Calculate computation time averaged over 100 simulations
time_count <- function(n,times){
  a<-matrix(data=NA,ncol=n,nrow=times)
  for (i in 1:times) {
    a[i,]=sample(1:n)
  }
  sum=rep(0,times)
  for (i in 1:times) {
    start=Sys.time()
    t<-fast_sorting(a[i, ])
    end=Sys.time()
    sum[i]=end-start
  }
  m=mean(sum)
  return(m)
}

## -----------------------------------------------------------------------------
x<-rep(0,5)
x[1]=time_count(10000,100)
for (i in 2:5) {
  x[i]=time_count(2*(i-1)*10000,100)
}
t<-rep(0,5)
t[1]=10000*log(10000)
for (i in 2:5) {
  t[i]=2*(i-1)*10000*log(2*(i-1)*10000)
}

## -----------------------------------------------------------------------------
#regression and graphically show the results
library(ggplot2)
data<-data.frame(x=x,t=t)
ggplot(data = data,aes(x=t,y=x))+geom_point(color="black")+geom_smooth(method = lm)

## -----------------------------------------------------------------------------
set.seed(22090)
m <- 1e4
#use the simple Monte Carlo method to estimate θ
smc <- replicate(1000, expr ={
  x <- runif(m, min = 0, max = 1)
  mean(exp(x))
} )

#use the antithetic variate approach to estimate θ
anti <- replicate(1000, expr ={
  y <- runif(m/2)
  u <- 1 - y
  mean((exp(y)+exp(u))/2)
})

#compute the variance of each estimator
var1 <- var(smc)
var2 <- var(anti)
print(c(var1, var2))

# Compute the empirical estimate of the percent reduction in variance using the antithetic variate.
(var1-var2)/var1

#compare the result of simulation with the theoretical value
print(c(mean(smc),mean(anti),exp(1)-exp(0)))


## -----------------------------------------------------------------------------
x <- seq(1,10,0.01)
y <- x^2 * exp(-x^2/2)/sqrt(2*pi)
plot(x, y, type = "l", ylim=c(0,1))
f1 <- 2*dnorm(x, 1)
lines(x, f1, lty = 2)
f2 <- dgamma(x - 1, 7/4, 2)
lines(x, f2, lty = 3)
legend("topright", inset = 0.02, legend = c("g(x)", "f1", "f2"), lty = 1:3)

## -----------------------------------------------------------------------------
plot(x, y/f1, type = "l", lty = 2)
lines(x, y/f2, lty = 3)
legend("topright", inset = 0.02, legend = c("g/f1", "g/f2"), lty = 2:3)

## -----------------------------------------------------------------------------
set.seed(22090)
M <- 10000
k <- 5
m <- M/k
g <- function(x) exp(-x)/(1 + x^2)
f <- function(x) (k/(1 - exp(-1))) * exp(-x)
va <- numeric(k)
si <- numeric(k)
for (i in 1:k){
  u <- runif(m, (i-1)/k, i/k)
  x <- -log(1 - (1 - exp(-1)) * u)
  gf <- g(x)/f(x)
  si[i] <- mean(gf)
  va[i] <- var(gf)
}
sum(si)
mean(va)
sqrt(mean(va))

## -----------------------------------------------------------------------------
set.seed(22090)
M <- 10000
k <- 1
m <- M/k
g <- function(x) exp(-x)/(1 + x^2)
f <- function(x) (k/(1 - exp(-1))) * exp(-x)
va <- numeric(k)
si <- numeric(k)
for (i in 1:k){
  u <- runif(m, (i-1)/k, i/k)
  x <- -log(1 - (1 - exp(-1)) * u)
  gf <- g(x)/f(x)
  si[i] <- mean(gf)
  va[i] <- var(gf)
}
sum(si)
mean(va)
sqrt(mean(va))

## -----------------------------------------------------------------------------
set.seed(22090)
n <- 100
CI <- replicate(10000, expr={
  # generate a random sample X1,...,Xn
  x <- rlnorm(n)
  # transform X to a normal variable Y
  y <- log(x)
  ybar <- mean(y)
  se <- sd(y)/sqrt(n)
  # construct a 95% confidence interval for the parameter µ
  ybar + se * qnorm(c(0.025,0.975))
})
lcl <- CI[1, ]
ucl <- CI[2, ]
mean(lcl < 0 & ucl >0)

## -----------------------------------------------------------------------------
count5test <- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0)
return(as.integer(max(c(outx, outy)) > 5))
}

## -----------------------------------------------------------------------------
# repeat the simulation and compute the F test
set.seed(22090)
sigma1 <- 1
sigma2 <- 1.5
m <- 10000
power <- mean(replicate(m, expr={
  x <- rnorm(20, 0, sigma1)
  y <- rnorm(20, 0, sigma2)
  count5test(x, y)
}))
Ftest <- mean(replicate(m, expr={
  x <- rnorm(20, 0, sigma1)
  y <- rnorm(20, 0, sigma2)
  Fp <- var.test(x, y)$p.value
  as.integer(Fp <= 0.055)
}))
c(power,Ftest)

## -----------------------------------------------------------------------------
#Compare the power of the Count Five test and F test
set.seed(22090)
sigma1 <- 1
sigma2 <- 1.5
m <- 10000
for (n in c(20, 50, 100, 200, 500, 1000)){
  compare <- replicate(m, expr = {
    x <- rnorm(n, 0, sigma1)
    y <- rnorm(n, 0, sigma2)
    c5t <- count5test(x, y)
    Fp <- var.test(x, y)$p.value
    Ftest <- as.integer(Fp <= 0.055)
    c(c5t, Ftest)
  })
  cat(n, rowMeans(compare), "\n")
}

## -----------------------------------------------------------------------------
Xbar <- 6510
Ybar <- 6760
m <- 10000
p <- 1/(m+m)*(Xbar+Ybar)
Z <- (Xbar-Ybar)/sqrt(p*(1-p))*sqrt(m/2)/10000
Z

## -----------------------------------------------------------------------------
library(boot)
set.seed(22090)
times <- aircondit
x <- as.matrix(times)
B <- 1e4
#sample generation function
stat <- function(x, i) return(1/mean(x[i, ]))
#data analysis function
b <- function(stat, B){
  boot(x, statistic = stat, B)
}
b(stat, B)

## -----------------------------------------------------------------------------
detach(package:boot)
rm(list = ls())

## -----------------------------------------------------------------------------
library(boot)
set.seed(22090)
times <- aircondit
x <- as.matrix(times)
B <- 1e4
#sample generation function
stat2 <- function(x, i) return(mean(x[i, ]))
#data analysis function
b <- function(stat, B){
  boot(x, statistic = stat, B)
}
boot.out <- b(stat2, B)
boot.out

## -----------------------------------------------------------------------------
boot.ci(boot.out, conf = 0.95, type = c("norm","basic", "perc", "bca"))

## -----------------------------------------------------------------------------
hist(boot.out$t, main = "")

## -----------------------------------------------------------------------------
detach(package:boot)
rm(list = ls())

## -----------------------------------------------------------------------------
b<-1;n<-10;m<-1e3;mu<-0;sigma<-1
library(boot)
set.seed(22090)
boot.mean <- function(x,i) mean(x[i])
ci.norm<-ci.basic<-ci.perc<-matrix(NA,m,2)
for(i in 1:m){
  x = rnorm(n, mean = mu,sd = sigma)
  de = boot(data = x,statistic = boot.mean,R=1e4)
  ci <- boot.ci(de,type=c("norm","basic","perc"))
  ci.norm[i,]<-ci$norm[2:3]
  ci.basic[i,]<-ci$basic[4:5]
  ci.perc[i,]<-ci$percent[4:5]
}
cat('norm =',mean(ci.norm[,1]<=mu & ci.norm[,2]>=mu),
    'basic =',mean(ci.basic[,1]<=mu & ci.basic[,2]>=mu),
    'perc =',mean(ci.perc[,1]<=mu & ci.perc[,2]>=mu))

knitr::kable(data.frame("type"=c("norm","basic","perc"),
                      "p_miss_on_left"=c(mean(ci.norm[,1]>mu),
                                         mean(ci.basic[,1]>mu),
                                         mean(ci.perc[,1]>mu)),
                      "p_miss_on_right"=c(mean(ci.norm[,2]<mu),
                                        mean(ci.basic[,2]<mu),
                                        mean(ci.perc[,2]<mu))))

## -----------------------------------------------------------------------------
# Data generation
library(bootstrap)
attach(scor)

# Construct the function of Obtaining the jackknife estimates
jack_est<-function(data){
  x <- as.matrix(data)
  rowl <- nrow(x)
  theta.jack <- numeric(rowl)
  lambda <- eigen(cov(x))$values
  lambda1 <- max(eigen(cov(x))$values)
  thetahat <- lambda1/sum(lambda)
  for (i in 1:rowl){
    x.jack <- x[-i, ]
    v <- eigen(cov(x.jack))$values
    theta.jack[i] <- max(v/sum(v))
  }
  bias <- (rowl - 1) * (mean(theta.jack) - thetahat)
  se <- sqrt((rowl - 1)/rowl * sum((theta.jack - mean(theta.jack))^2))
  knitr::kable(data.frame("object"=c("est","bias","se"), "value"=c(thetahat,bias,se)))
}

# apply the function
jack_est(scor)

## -----------------------------------------------------------------------------
detach(package:bootstrap)
rm(list = ls())

## -----------------------------------------------------------------------------
# Data generation
library(DAAG)
attach(ironslag)

# Construct the function of comparing 4 different models using leave-two-out cross validation
compare_def <- function(data){
  data1<-data
  colnames(data1)<-c("x0","y0")
  n <- length(data1$x0)
  e1 <- e2 <- e3 <- e4 <- numeric(n*(n-1)/2)
  k<-0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      k<-k+1
      y <- data1$y0[-c(i,j)]
      x <- data1$x0[-c(i,j)]
      J1 <- lm(y ~ x)
      yhat1 <- J1$coef[1] + J1$coef[2] * data1$x0[c(i,j)]
      e1[k] <- sum((data1$y0[c(i,j)] - yhat1)^2)
      J2 <- lm(y ~ x + I(x^2))
      yhat2 <- J2$coef[1] + J2$coef[2] * data1$x0[c(i,j)]+ J2$coef[3] * data1$x0[c(i,j)]^2
      e2[k] <- sum((data1$y0[c(i,j)] - yhat2)^2)
      J3 <- lm(log(y) ~ x)
      logyhat3 <- J3$coef[1] + J3$coef[2] * data1$x0[c(i,j)]
      yhat3 <- exp(logyhat3)
      e3[k] <- sum((data1$y0[c(i,j)] - yhat3)^2)
      J4 <- lm(log(y) ~ log(x))
      logyhat4 <- J4$coef[1] + J4$coef[2] * log(data1$x0[c(i,j)])
      yhat4 <- exp(logyhat4)
      e4[k] <- sum((data1$y0[c(i,j)] - yhat4)^2)
    }
  }
l <- c(mean(e1),mean(e2),mean(e3),mean(e4))
cat("The estimates for prediction error for the four models(Linear,Quadratic,Exponential,Log-Log) are, respectively \n",l,"\n")
}
# Apply the function
compare_def(ironslag)


## -----------------------------------------------------------------------------
L2 <- lm(magnetic ~ chemical + I(chemical^2))
L2

## -----------------------------------------------------------------------------
detach(package:DAAG)
rm(list = ls())

## -----------------------------------------------------------------------------
# Data generation 1
set.seed(22090)
library(MASS)
sample1 <- function(n, mu, Sigma){
  n <- 50
  mu <- c(0,0)
  Sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  x <- mvrnorm(n, mu, Sigma)
  return (x)
}

# Construct the function of the bivariate Spearman rank correlation test
spear.test <- function (x,y){
  stest <- cor.test(x, y, method = "spearman")
  n <- length(x)
  res <- replicate(N, expr = {
    k <- sample(1:n)
    cor.test(x, y[k], method = "spearman")$estimate
  })
  res1 <- c(stest$estimate, res)
  pval <- mean(as.integer(stest$estimate <= res1))
  return(list(spear.rho = stest$estimate, pvalue = pval))
}

# Apply the function and compare the results(sample1)
N <- 500
samp <- sample1(n, mu, Sigma)
cor.test(samp[, 1], samp[, 2], method = "spearman")
spear.test(samp[, 1], samp[, 2])

# Data generation 2
set.seed(22090)
library(MASS)
sample2 <- function(n, mu, Sigma){
  n <- 50
  mu <- c(0,0)
  Sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  x <- exp(mvrnorm(n, mu, Sigma))
  return (x)
}

# Apply the function and compare the results(sample2)
samp <- sample2(n, mu, Sigma)
cor.test(samp[, 1], samp[, 2], method = "spearman")
spear.test(samp[, 1], samp[, 2])

## -----------------------------------------------------------------------------
detach(package:MASS)
rm(list = ls())

## -----------------------------------------------------------------------------
# Construct the function
rw.Metropolis <- function(sigma, x0, N){
  x <-  numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= exp(abs(x[i-1]) - abs(y)))
      x[i] <- y
    else{
      x[i] <- x[i-1]
      k <- k+1
    }
  }
  return(list(x = x, k = k))
}

## -----------------------------------------------------------------------------
set.seed(22090)
N <- 2000
sigma <- c(.05, .5, 2, 16)
x0 <- rnorm(1,0,1)
rw1 <- rw.Metropolis(sigma[1], x0, N)
rw2 <- rw.Metropolis(sigma[2], x0, N)
rw3 <- rw.Metropolis(sigma[3], x0, N)
rw4 <- rw.Metropolis(sigma[4], x0, N)
print(c(rw1$k, rw2$k, rw3$k, rw4$k))
knitr::kable(data.frame("σ"=c(0.05, 0.5, 2 ,4), 
                        "acceptance rates"=1-c(rw1$k, rw2$k, rw3$k, rw4$k)/N,
                        "rejection rates"=c(rw1$k, rw2$k, rw3$k, rw4$k)/N))

## -----------------------------------------------------------------------------
plot(rw1$x, type="l", main = "σ = 0.05")
plot(rw2$x, type="l", main = "σ = 0.5")
plot(rw3$x, type="l", main = "σ = 2")
plot(rw4$x, type="l", main = "σ = 4")

## -----------------------------------------------------------------------------
m <- 100
y1 <- rw1$x[(m + 1):N]
y2 <- rw2$x[(m + 1):N]
y3 <- rw3$x[(m + 1):N]
y4 <- rw4$x[(m + 1):N]

## -----------------------------------------------------------------------------
p <- ppoints(100)
y <- qexp(p, 1)
z <- c(-rev(y), y)
Q1 <- quantile(y1, p)
qqplot(z, Q1, cex = 0.5)
abline(0, 1)
Q2 <- quantile(y2, p)
qqplot(z, Q2, cex = 0.5)
abline(0, 1)
Q3 <- quantile(y3, p)
qqplot(z, Q3, cex = 0.5)
abline(0, 1)
Q4 <- quantile(y4, p)
qqplot(z, Q4, cex = 0.5)
abline(0, 1)

## -----------------------------------------------------------------------------
#initialize constants and parameters
N <- 5000 #length of chain
burn <- 1000 #burn-in length
X <- matrix(0, N, 2) #the chain, a bivariate sample
rho <- 0.9
mu1 <- 0
mu2 <- 0
sigma1 <- 1
sigma2 <- 1
s1 <- sqrt(1-rho^2)*sigma1
s2 <- sqrt(1-rho^2)*sigma2

###### generate the chain #####

X[1, ] <- c(mu1, mu2) #initialize

for (i in 2:N) {
  x2 <- X[i-1, 2]
  m1 <- mu1 + rho * (x2 - mu2) * sigma1/sigma2
  X[i, 1] <- rnorm(1, m1, s1)
  x1 <- X[i, 1]
  m2 <- mu2 + rho * (x1 - mu1) * sigma2/sigma1
  X[i, 2] <- rnorm(1, m2, s2)
}

b <- burn + 1
x <- X[b:N, ] #discard a suitable burn-in sample
X <- x[, 1]
Y <- x[, 2]
L <- lm(Y ~ X)
L
summary(L)

## -----------------------------------------------------------------------------
# compare sample statistics to parameters
colMeans(x)
cov(x)
cor(x)
plot(X,Y, main="", cex=.5)

## -----------------------------------------------------------------------------
plot(L$fit, L$res, cex = 0.5)
abline(h = 0)
qqnorm(L$res, cex = 0.5)
qqline(L$res)

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
  # psi[i,j] is the statistic psi(X[i,1:j]) for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi) #row means
  B <- n * var(psi.means) #between variance est.
  psi.w <- apply(psi, 1, "var") #within variances
  W <- mean(psi.w) #within est.
  v.hat <- W*(n-1)/n + (B/n) #upper variance est.
  r.hat <- v.hat / W #G-R statistic
  return(r.hat)
}

## -----------------------------------------------------------------------------
x0 <- 1.2 #parameter of proposal distribution
k <- 4 #number of chains to generate
n <- 15000 #length of chains
b <- 1000 #burn-in length

#different variances of the proposal distribution
sigma <- c(.05, .5, 2, 4)

#generate the chains
X <- matrix(0, nrow=k, ncol=n)
for (i in 1:k)
  X[i, ] <- rw.Metropolis(sigma[i], x0, N)$x

#compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
  psi[i,] <- psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi))

#plot psi for the four chains
for (i in 1:k)
  plot(psi[i, (b+1):n], type="l",xlab=i, ylab=bquote(psi))

#plot the sequence of R-hat statistics
rhat <- rep(0, n)
for (j in (b+1):n)
  rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
# Construct the function
perm1<-function(x,y,m,med=1,t){
  if(med==1){#$\alpha=\beta=0$，M is not related to X, Y is not related to M
    n <- length(x)
    beta.hat1 <- beta.se1 <- p.val11 <- beta.hat2 <- beta.se2 <- p.val12 <- p <- numeric(t)
    coe1 <- summary(lm(m~x))$coef
    coe2 <- summary(lm(y~x+m))$coef
    beta.hat1[1] <- coe1[2,1]
    beta.se1[1] <- coe1[2,2] 
    beta.hat2[1] <- coe2[2,1]
    beta.se2[1] <- coe2[2,2]
    p[1]=beta.hat1[1] * beta.hat2[1] / beta.se1[1] / beta.se2[1]
    for(i in 2:t){
      xt <- x[sample(1:n)]
      mt <- m[sample(1:n)]
      coe1 <- summary(lm(m~xt))$coef
      coe2 <- summary(lm(y~x+mt))$coef
      beta.hat1[i] <- coe1[2,1]
      beta.se1[i] <- coe1[2,2] 
      beta.hat2[i] <- coe2[2,1]
      beta.se2[i] <- coe2[2,2] 
      p[i]=beta.hat1[i] * beta.hat2[i] / beta.se1[i] / beta.se2[i]
    }
    pvalue <- 2*min(length(p[p<p[1]]),length(p[p>p[1]]))/(t-1)
    return(pvalue)
  }
  if(med==2){#$\alpha=0,\beta=1$
    n <- length(x)
    beta.hat1 <- beta.se1 <- p.val11 <- beta.hat2 <- beta.se2 <- p.val12 <- p <- numeric(t)
    coe1 <- summary(lm(m~x))$coef
    coe2 <- summary(lm(y~x+m))$coef
    beta.hat1[1] <- coe1[2,1]
    beta.se1[1] <- coe1[2,2] 
    beta.hat2[1] <- coe2[2,1]
    beta.se2[1] <- coe2[2,2]
    p[1]=beta.hat1[1] * beta.hat2[1] / beta.se1[1] / beta.se2[1]
    for(i in 2:t){
      xt <- x[sample(1:n)]
      coe1 <- summary(lm(m~xt))$coef
      coe2 <- summary(lm(y~x+m))$coef
      beta.hat1[i] <- coe1[2,1]
      beta.se1[i] <- coe1[2,2] 
      beta.hat2[i] <- coe2[2,1]
      beta.se2[i] <- coe2[2,2] 
      p[i]=beta.hat1[i] * beta.hat2[i] / beta.se1[i] / beta.se2[i]
    }
    pvalue <- 2*min(length(p[p<p[1]]),length(p[p>p[1]]))/(t-1)
    return(pvalue)
  }
  if(med==3){#$\alpha=1,\beta=0$
    n <- length(x)
    beta.hat1 <- beta.se1 <- p.val11 <- beta.hat2 <- beta.se2 <- p.val12 <- p <- numeric(t)
    coe1 <- summary(lm(m~x))$coef
    coe2 <- summary(lm(y~x+m))$coef
    beta.hat1[1] <- coe1[2,1]
    beta.se1[1] <- coe1[2,2] 
    beta.hat2[1] <- coe2[2,1]
    beta.se2[1] <- coe2[2,2]
    p[1]=beta.hat1[1] * beta.hat2[1] / beta.se1[1] / beta.se2[1]
    for(i in 2:t){
      mt <- m[sample(1:n)]
      coe1 <- summary(lm(m~x))$coef
      coe2 <- summary(lm(y~x+mt))$coef
      beta.hat1[i] <- coe1[2,1]
      beta.se1[i] <- coe1[2,2] 
      beta.hat2[i] <- coe2[2,1]
      beta.se2[i] <- coe2[2,2] 
      p[i]=beta.hat1[i] * beta.hat2[i] / beta.se1[i] / beta.se2[i]
    }
    pvalue <- 2*min(length(p[p<p[1]]),length(p[p>p[1]]))/(t-1)
    return(pvalue)
  }
}

## -----------------------------------------------------------------------------
# Apply the function
library(alr4)
t <- 1e3
data <- water
x <- data$APMAM
y <- data$APSAB
m <- data$BSAAM
perm1(x,y,m,1,t)

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
set.seed(22090)

# Construct the function
computing_alpha <- function(N, b1, b2, b3, f0){
  g <- function(alpha){
    tmp <- exp(-alpha-b1*x1-b2*x2-b3*x3)
    p <- 1/(1+tmp)
    mean(p) - f0
  }
  solution <- uniroot(g, c(-15,0))
  alpha <- solution$root
  return(alpha)
}

# Apply the function
N <- 1e6
i <- 1
valpha <- numeric(4)
b1 <- 0
b2 <- 1
b3 <- -1
x1 <- rpois(N, lambda = 1)
x2 <- rexp(N, 1)
x3 <- sample(0:1,N,replace=TRUE)
for (f0 in c(0.1, 0.01, 0.001, 0.0001)){
  valpha[i] <- computing_alpha(N, b1, b2, b3, f0)
  i <- i+1
}
valpha

## -----------------------------------------------------------------------------
# Draw the scatter diagram 
x <- numeric(4)
j <- 1
for (f0 in c(0.1, 0.01, 0.001, 0.0001)){
  x[j] <- -log10(f0)
  j <- j+1
}
plot(x, valpha, type = "l")

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
# construct the likelihood function
LL <- function(lambda, u, v ){
  n <- length(u)
  f <- numeric(n)
  for (i in 1:n){
    f[i] <- exp(-u[i]*lambda)- exp(-v[i]*lambda)
  }
  return(sum(log(f)))
}

## -----------------------------------------------------------------------------
u <- c(11, 8, 27, 13, 16, 0, 23, 10, 24, 2)
v <- c(12, 9, 28, 14, 17, 1, 24, 11, 25, 3)
opt <- optimize(LL, interval = c(0.01,1), u=u, v=v, maximum = T)
opt

## -----------------------------------------------------------------------------
# Construct the iteration function 
lambda1<-function(lambda,u,v){
  n<-length(u)
  return(n/(sum((exp(-lambda*u)*(u+1/lambda)-exp(-lambda*v)*(v+1/lambda))/(exp(-lambda*u)-exp(-lambda*v)))))
}

## -----------------------------------------------------------------------------
x <- mean((u+v)/2)
while (abs(x-lambda1(x,u,v))>=0.000001) {
  x <- lambda1(x,u,v)
}
x

## -----------------------------------------------------------------------------
1 == "1" 
-1 < FALSE
"one" < 2

## -----------------------------------------------------------------------------
b <- c(1:12)
dim(b)

## -----------------------------------------------------------------------------
# data frame with 0 rows
data.frame(a = integer(), b = logical())

# data frame with 0 columns and 3 rows
data.frame(row.names = 1:3)  

# data frame with 0 columns and 0 rows
data.frame()


## -----------------------------------------------------------------------------
scale01 <- function(x) {
rng <- range(x, na.rm = TRUE)
(x - rng[1]) / (rng[2] - rng[1])
}

## -----------------------------------------------------------------------------
df <- data.frame(lapply(iris, function(x) if (is.numeric(x)) scale01(x) else x))
head(df, 30)

## -----------------------------------------------------------------------------
vapply(cars, sd, numeric(1))

## -----------------------------------------------------------------------------
vapply(iris[vapply(iris, is.numeric, logical(1))],
       sd, 
       numeric(1))

## -----------------------------------------------------------------------------
library(Rcpp)
sourceCpp(code='
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix gibbsC(int N, int burn, double mu1, double mu2,double sigma1,double sigma2,double rho) {
  NumericMatrix mat(N, 2);
  double x = mu1, y = mu2;
  mat(1,0) = x;
  mat(1,1) = y;
  double s1 = sqrt(1-pow(rho,2))*sigma1;
  double s2 = sqrt(1-pow(rho,2))*sigma2;
  for(int i = 1; i < N; i++) {
    for(int j = 0; j < burn; j++) {
      y = mat(j-1, 2);
      double m1 = 0 + rho * (y - 0) * sigma1/sigma2;
      x = rnorm(1, m1, s1)[0];
      double m2 = 0 + rho * (x - 0) * sigma2/sigma1;
      y = rnorm(1, m2, s2)[0];
    }
    mat(i, 0) = x;
    mat(i, 1) = y;
  }
  return(mat);
}  
')

## -----------------------------------------------------------------------------
gibbs_R <- function(N, burn, mu1, mu2, sigma1, sigma2, rho){
  X <- matrix(0, N, 2)
  s1 <- sqrt(1-rho^2)*sigma1
  s2 <- sqrt(1-rho^2)*sigma2
  X[1, ] <- c(mu1, mu2) #initialize

  for (i in 2:N) {
    x2 <- X[i-1, 2]
    m1 <- mu1 + rho * (x2 - mu2) * sigma1/sigma2
    X[i, 1] <- rnorm(1, m1, s1)
    x1 <- X[i, 1]
    m2 <- mu2 + rho * (x1 - mu1) * sigma2/sigma1
    X[i, 2] <- rnorm(1, m2, s2)
  }

  b <- burn + 1
  x <- X[b:N, ]
  return(x)
}

## -----------------------------------------------------------------------------
par(mfrow = c(1, 2))
x2 <- gibbsC(5000, 1000, 0, 0, 1, 1, 0.9)
X2 <- x2[, 1]
Y2 <- x2[, 2]
qqplot(X2, Y2, main="QQ Plot_C")

## -----------------------------------------------------------------------------

x1 <- gibbs_R(5000, 1000, 0, 0, 1, 1, 0.9)
X1 <- x1[, 1]
Y1 <- x1[, 2]
qqplot(X1, Y1, main="QQ Plot_R")

