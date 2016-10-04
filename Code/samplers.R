# binomhm <- function(N0,Y0,N,thetai,thetaj,thetaij,sigma2,tau2r,tau2c,tau2c0,thetaj0,thetaj00,ridx0,cidx0,dMr,dMc,nmc,whichkeep,whichkeepi,dispint,hstr) {
#   THETAIJ <- matrix(0,length(whichkeep),nmc)
#   THETAI <- matrix(0,length(whichkeepi),nmc)
#   THETAJ <- matrix(0,length(thetaj),nmc)
#   SIGMA <- matrix(0,nmc,1)
#   #ridx0 <- factor(ridx0)
#   #cidx0 <- factor(cidx0)
#   kappa <- Y0-N0/2
#   t1 <- proc.time()
#   for (tm in 1:nmc) {
#     if (tm%%dispint==0) {
#       print(tm)
#       t2 <- proc.time()
#       print(t2-t1)
#       t1 <- proc.time()
#     }
#     
#     # sample polya-gamma
#     omega <- rpg(N,N0,thetaij)
#     
#     if (hstr==2) {
#       theta0 <- thetai[ridx0] + thetaj[cidx0]
#     }
#     if (hstr==3) {
#       theta0 <- thetai[ridx0]
#     }
#     
#     # sample thetaij
#     if (hstr==0) {
#       #sigma2 <- 100
#       sij <- (omega + sigma2^(-1))^(-1)
#       mij <- sij * (kappa)
#       thetaij <- mij+sqrt(sij)*zrnorm(N)
#       THETAIJ[,tm] <- thetaij[whichkeep]
#     }
#     if (hstr!=0) {
#       sij <- (omega + (sigma2)^(-1))^(-1)
#       mij <- sij * (kappa + ((sigma2)^(-1))*(theta0))
#       thetaij <- mij+sqrt(sij)*zrnorm(N)
#       THETAIJ[,tm] <- thetaij[whichkeep]
#     }
#     
#     if (hstr==2) {
#       # sample sigma2
#       res <- thetaij - (thetai[ridx0]+thetaj[cidx0])
#       ssr <- sum(res^2)
#       sigma2 <- rgamma(1,(1+N)/2,(1+ssr)/2)
#       SIGMA[tm] <- sigma2
#       
#       # sample thetai and thetaj
#       #res <- tapply(thetaij - thetaj[cidx0],ridx0,sum)
#       res <- t(t(thetaij-thetaj[cidx0])%*%dMr)
#       si <- (Nc/sigma2+1/tau2r)^(-1)
#       mi <- si*(res/sigma2)
#       thetai <- mi+sqrt(si)*zrnorm(Nr)
#       thetai <- as.matrix(thetai)
#       THETAI[,tm] <- thetai[whichkeepi]
#       
#       ssr <- sum(thetai^2)
#       tau2r <- rgamma(1,(Nr-1)/2,(ssr)/2)
#       
#       #res <- tapply(thetaij - thetai[ridx0],cidx0,sum)
#       res <- t(t(thetaij-thetai[ridx0])%*%dMc)
#       sj <- (Nr/sigma2+1/tau2c)^(-1)
#       mj <- sj*(res/sigma2+thetaj0/tau2c)
#       thetaj <- mj+sqrt(sj)*zrnorm(Nc)    
#       thetaj <- as.matrix(thetaj)
#       THETAJ[,tm] <- thetaj
#       
#       ssr <- sum((thetaj-thetaj0)^2)
#       tau2c <- rgamma(1,(Nc-1)/2,(ssr)/2)
#       
#       sj <- (Nc/tau2c+1/tau2c0)^(-1)
#       mj <- sj*(sum(thetaj)/tau2c+thetaj00/tau2c0)
#       thetaj0 <- mj + sqrt(sj)*zrnorm(1) 
#     }
#     if (hstr==3) {
#       # sample sigma2
#       res <- thetaij - (thetai[ridx0]+thetaj[cidx0])
#       ssr <- sum(res^2)
#       sigma2 <- rgamma(1,(N-1)/2,(ssr)/2)
#       SIGMA[tm] <- sigma2      
#       
#       # sample thetai and thetaj
#       #res <- tapply(thetaij - thetaj[cidx0],ridx0,sum)
#       res <- t(t(thetaij)%*%dMr)
#       si <- (Nc/sigma2+1/tau2r)^(-1)
#       mi <- si*(res/sigma2+thetaj0/tau2r)
#       thetai <- mi+sqrt(si)*zrnorm(Nr)
#       thetai <- as.matrix(thetai)
#       THETAI[,tm] <- thetai[whichkeepi]
#       
#       ssr <- sum((thetai-thetaj0)^2)
#       tau2r <- rgamma(1,(Nr-1)/2,(ssr)/2)
#       
#       sj <- (Nr/tau2c+1/tau2c0)^(-1)
#       mj <- sj*(sum(thetaj)/tau2c+thetaj00/tau2c0)
#       thetaj0 <- mj + sqrt(sj)*zrnorm(1) 
#     }     
#   }
#   
#  return(list(THETAIJ=THETAIJ,THETAI=THETAI,THETAJ=THETAJ,SIGMA=SIGMA))
#   
# }

binomhm_simple <- function(N,Y,kappa,mu,mu0,mu00,tau20,Tau2,n,nmc,dispint) {
  MU <- matrix(0,nmc,n)
  MU0 <- matrix(0,nmc,1)
  TAU <- matrix(0,nmc,1)
  t1 <- proc.time()
  for (tm in 1:nmc) {
    if (tm%%dispint==0) {
      print(tm)
      t2 <- proc.time()
      print(t2-t1)
      t1 <- proc.time()
    }
    
    omega <- rpg(n,N,mu)
    
    sij <- (omega + Tau2^(-1))^(-1)
    mij <- sij * (kappa + ((Tau2)^(-1))*mu0)
    mu <- mij+sqrt(sij)*rnorm(n)
    MU[tm,] <- mu
    
    
    res <- mu - mu0
    ssr <- sum(res^2)
    Tau2 <- rgamma(1,(n-1)/2,ssr/2)
    
    s <- (n/Tau2+1/tau20)^(-1)
    m <- s*(sum(mu)/Tau2+mu00/tau20)
    mu0 <- rnorm(1,m,sqrt(s))
    
    TAU[tm] <- Tau2
    MU[tm,] <- mu
    MU0[tm] <- mu0
  }
  return(list(TAU=TAU,MU=MU,MU0=MU0))
}


makedummy <- function(df) {
  n <- nrow(df)
  nlevels <- sapply(df, nlevels)
  i <- rep(seq_len(n), ncol(df))
  j <- unlist(lapply(df, as.integer)) +
    rep(cumsum(c(0, head(nlevels, -1))), each = n)
  x <- 1
  M <- sparseMatrix(i = i, j = j, x = x)
  return(M)
} 

# albert and chib sampler for probit
albert_chib <- function(YS,NS,nmc,mu0,Tau2,disp_int) {
  ny <- length(YS)
  MU <- matrix(0,nmc,ny)
  
  MU[1,] <- mu0
  t0 <- proc.time()
  for (t in 2:nmc) {
    if (t%%disp_int==0) {
      print(t)
      t1 <- proc.time()
      print(t1-t0)
      t0 <- proc.time()
    }
    # the Robert accept-reject sampler is very slow in this case, so I've used Geweke
    #w1 <- rtmvnorm(N-Y,mu,1,lower=-Inf,upper=0)
    #w2 <- rtmvnorm(Y,mu,1,lower=0,upper=Inf)
    
    # first for Y=1
    for (j in 1:ny) {    
      Y <- YS[j]
      N <- NS[j]
      mu <- MU[t-1,j]
      w1 <- rtruncnorm(N-Y,a=-Inf,b=0,mu,1)
      w2 <- rtruncnorm(Y,a=0,b=Inf,mu,1)
      wbar <- (sum(w1)+sum(w2))/N
      s <- (N+1/Tau2)^(-1)
      m <- s*(N*wbar)
      mu <- rnorm(1,m,sqrt(s))
      MU[t,j] <- mu
    }
    
  }
  
  PHI.AC <- matrix(0,ny,100)
  ES.AC <- matrix(0,ny,1)
  PHI.AC.L <- matrix(0,ny,nmc-1)  
  for (j in 1:ny) {
    #f <- ar(MU[,j],order.max=1)
    #PHI.AC[j] <- f$ar
    PHI.AC[j,] <- autocorr(as.mcmc(MU[,j]),lags = seq(100))
    ES.AC[j] <- effectiveSize(as.mcmc(MU[,j]))
    PHI.AC.L[j,] <- autocorr(as.mcmc(MU[,j]),lags = seq(nmc-1)) 
  }
  ret <- list()
  ret$MU <- MU
  ret$PHI <- PHI.AC
  ret$ES <- ES.AC
  ret$PHI.L <- PHI.AC.L  
  return(ret)
}


# binompg_hmc <- function(n,kappa,mu,B,dohmc,dispint,nmc) {
#   MU <- matrix(0,nmc,1)
#   acc <- matrix(0,nmc,1)
#   t1 <- proc.time()
#   for (tm in 1:nmc) {
#     if (tm%%dispint==0) {
#       print(tm)
#       t2 <- proc.time()
#       print(t2-t1)
#       t1 <- proc.time()
#     }
#     
#     omega <- rpg(1,n,mu)
#     
#     s <- (omega + B^(-1))^(-1)
#     mn <- s * kappa
#     mu <- rnorm(1,mn,sqrt(s))
#     #MU[tm] <- mu
#     
#     if (dohmc) {
#       # hmc step
#       # sample p(0)
#       #m <- (omega + B^(-1))^(-1)
#       m <- 1
#       #m <- 10*n
#       th <- log(n)
#       p0 <- rnorm(1,0,sqrt(m))
#       # compute constants
#       C1 <- p0
#       a <- kappa
#       b <- (omega + B^(-1))
#       C2 <- mu - a/b 
#       # compute state at time th
#       pt <- C1*cos(th*sqrt(b)/sqrt(m))-C2*sqrt(b*m)*sin(th*sqrt(b)/sqrt(m))
#       mut <- a/b + C2*cos(th*sqrt(b)/sqrt(m))+C1*sqrt(1/(b*m))*sin(th*sqrt(b)/sqrt(m))
#       # calculate U and K
#       U0 <- (mu-a/b)^2/(2/b)
#       Ut <- (mut-a/b)^2/(2/b)
#       K0 <- (p0^2)/(2*m)
#       Kt <- (pt^2)/(2*m)
#       ar <- exp(-Ut+U0-Kt+K0)
#       ar <- min(ar,1)
#       ac <- (runif(1)<ar)
#       acc[tm] <- ac
#       if (ac) {
#         mu <- mut
#       } 
#     } else {
#       acc[tm] <- NA  
#     }
#     
#     
#     MU[tm] <- mu
#   }
#   return(list(MU=MU,acc=acc))
# }
# 
# 
# binompg_inflate <- function(n,y,kappa,mu,b,B,dispint,nmc,d.x,d.y) {
#   MU <- matrix(0,nmc,1)
#   acc <- matrix(0,nmc,1)
#   t1 <- proc.time()
#   for (tm in 1:nmc) {
#     if (tm%%dispint==0) {
#       print(tm)
#       t2 <- proc.time()
#       print(t2-t1)
#       t1 <- proc.time()
#     }
#     
#     mu0 <- mu
#     # construct proposal
#     omega <- rpg(1,ceiling(log(n+1)),mu)
#     
#     s <- (omega + B^(-1))^(-1)
#     #mn <- s * kappa
#     
#     mu <- rnorm(1,mn,sqrt(s))
#     #MU[tm] <- mu
# 
#     # compute the likelihood ratio
#     p0 <- exp(mu0)/(1+exp(mu0))
#     pp <- exp(mu)/(1+exp(mu))
#     lr <- y*(log(pp)-log(p0)) + (n-y)*(log(1-pp)-log(1-p0))
#     pr <- -(mu-b)^2/(2*B)+(mu0-b)^2/(2*B)
#     mf1 <- -(mu0-mn)^2/(2*s) + lpgdens(omega,d.x,d.y,ceiling(log(n+1)),mu)
#     mf2 <- -(mu-mn)^2/(2*s) + lpgdens(omega,d.x,d.y,ceiling(log(n+1)),mu0)
#     ar <- lr + pr + mf1 - mf2
#     ac <- runif(1) < exp(ar)
#     
#         
#     acc[tm] <- ac
#     if (ac) {
#       mu <- mu
#     } else {
#       mu <- mu0
#     }
#     
#     
#     MU[tm] <- mu
#   }
#   return(list(MU=MU,acc=acc))
# }

# lpgdens <- function(om,d.x,d.y,b,c) {
#   ld0 <- log(d.y[min(which(d.x>om))])
#   ld <- ld0 - c^2/2*om + b * log(cosh(c/2))
# }




logmcmc <- function(niter,y,n,tau)
{
  out   = rep(NA,niter)
  out[1] = 1
  for(i in 2:niter)
  {
    prop = out[i-1] + rnorm(1)
    acc  = targ(prop,y,n,tau) - targ(out[i-1],y,n,tau)
    if(-rexp(1) < acc)
    {
      out[i] = prop
    } else
    {
      out[i] = out[i-1]
    }
  }
  return(out)
}

targ = function(x,y,n,tau) 
{
  y*x - n*log(1+ exp(x)) - x^2/(2*tau)
}




