rm(list=ls(all=T))
#setwd('~/Dropbox (Personal)/Gibbsmixing/Code/')
if (length(grep("sailfish",Sys.info()["nodename"]))>0) {
  setwd('/ytmp/jej/Projects/gibbsmixing/Code/')
} else {
  #setwd('~/Dropbox (Personal)/Gibbsmixing/Code/')
  setwd('~/Dropbox/Gibbsmixing/Code/Final-AoS')
}
#setwd(script.dir)
require('tmvtnorm')
require('truncnorm')
require('coda')
require('BayesLogit')
require('rstan')
require('xtable')
require('reshape2')
require('ggplot2')
require('R.matlab')
require('ggmcmc')
source('samplers.R')

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

run_samplers <- T

set.seed(5297121)

calcac <- function(x) {
  y <- as.mcmc(x)
  return(autocorr(y,lags=c(1,seq(from=5,to=100,by=5))))
}

calces <- function(x) {
  y <- as.mcmc(x)
  return(effectiveSize(y))
}

# real data example : maxpoint conversions
mpdat <- readMat('Data/convdat.mat')
Ntr <- mpdat$Ntr
Ytr <- as.matrix(mpdat$Ytr)
Ytr <- Ytr[Ntr>0,]
Ntr <- Ntr[Ntr>0]
# Ytr <- Ytr[1:1000,]
# Ntr <- Ntr[1:1000]
Y <- Ytr[,46]
N <- Ntr
n <- length(Ntr)
tmp <- log(((Y+1)/N)/(1-(Y+1)/N))
tmp[N==1] <- mean(tmp[N>1])
mu <- tmp
mu0 <- mean(tmp)
mu00 <- -12
tau20 <- 49
Tau2 <- 4
nmc <- 10000
dispint <- 100
kappa <- Y-N/2  

if (run_samplers) {
fit <- binomhm_simple(N,Y,kappa,mu,mu0,mu00,tau20,Tau2,n,nmc,dispint)
AR <- matrix(0,n,21)
ES <- matrix(0,n,1)


for (j in 1:n) {
  if (j%%1000==0) {
    print(j)
  }
  #tmp <- ar(fit$MU[5000:10000,j],aic=FALSE,order.max=1)
  #AR[j] <- tmp$ar
  tmp <- autocorr(as.mcmc(fit$MU[5000:10000,j]),lags = c(1,seq(from=5,to=100,by=5)))
  AR[j,] <- c(tmp)
  ES[j] <- effectiveSize(fit$MU[5000:10000,j])
} 
  MU0 <- fit$MU0
  TAU <- fit$TAU
  MUsml <- fit$MU[,seq(from=1,to=n,by=100)]

  AR <- data.frame(AR)
  names(AR) <- paste('lag-',as.character(c(1,seq(from=5,to=100,by=5))),sep='')
  AR$site <- seq(n)

  save(AR,MU0,TAU,MUsml,ES,file='Outputs/PG_maxpt.RData')
  #save.image(file='Ouputs/PG_maxpt.RData')
} else {
  load('Outputs/PG_maxpt.RData')
}

ARl <- melt(AR,id='site')
names(ARl) <- c('site','lag','autocorrelation')
ARl$lag <- gsub("lag-","",ARl$lag)
ARl$lag <- as.factor(as.numeric(ARl$lag))

ESscl <- ES/5001
ESscl <- data.frame(ESscl)
names(ESscl) <- c('EffectiveSize')

png('Figures/PG_maxpt_acf_box.png',width=450,height=300)
ggplot(ARl,aes(y=autocorrelation,x=lag)) + geom_boxplot(outlier.shape=NA) + ggtitle('Polya Gamma') + scale_y_continuous(limits=c(-.25,1)) +
  theme(text=element_text(size=24)) + scale_x_discrete(breaks=c(1,seq(from=10,to=100,by=10)))
dev.off()

png('Figures/PG_maxpt_es_ecdf.png',width=450,height=300)
#ggplot(ESscl,aes(x=EffectiveSize)) + geom_histogram() + scale_x_continuous(limits=c(0,1)) + ggtitle('PG -- Effective size (proportion of samples)')
ggplot(ESscl,aes(x=EffectiveSize)) + stat_ecdf(geom="step",size=1.25) + scale_x_continuous(limits=c(0,1.5)) + ggtitle('Polya Gamma') + 
  theme(text=element_text(size=24),axis.title.y=element_blank())
dev.off()

png('Figures/PG_maxpt_es_hist.png',width=450,height=300)
#ggplot(ESscl,aes(x=EffectiveSize)) + geom_histogram() + scale_x_continuous(limits=c(0,1)) + ggtitle('PG -- Effective size (proportion of samples)')
ggplot(ESscl,aes(x=EffectiveSize)) + geom_histogram() + scale_x_continuous(limits=c(-0.1,1.5)) + ggtitle('Polya Gamma') + theme(text=element_text(size=24),axis.title.y=element_blank())
dev.off()


# png('Figures/PG_maxpt_hist.png',width=450,height=300)
# hist(AR,50,col='gray',main='PG',xlim=c(-.1,1))
# dev.off()
# 
# png('Figures/PG_maxpt_traces.png',width=600,height=600)
# par(mfrow=c(2,2),oma=rep(1.5,4),mar=rep(1.5,4))
# fake.ggplot(fit$MU[5000:10000,5],cex=.5,ylab='',xlab='t')
# fake.ggplot(fit$MU[5000:10000,5005],cex=.5,ylab='',xlab='t')
# fake.ggplot(fit$MU[5000:10000,15005],cex=.5,ylab='',xlab='t')
# fake.ggplot(fit$MU[5000:10000,50000],cex=.5,ylab='',xlab='t')
# dev.off()


##########################

if (run_samplers) {
smod <- "data {
  int<lower=0> n; //num sites
  int<lower=0> Y[n]; // Y vector
  int<lower=1> N[n]; // N vector
}

parameters {
real mu[n];  // vector of mus
real mu0; // grand mean
real<lower=0> Tau;   // std dev
}

model {
//Priors
Tau~uniform(0,100); 
mu0~normal(-12,7);
for(i in 1:n) {
  mu[i] ~ normal(0,Tau);
}
// Likelihood
for(i in 1:n) {
  Y[i] ~ binomial(N[i],inv_logit(mu[i]));
} 
}
"
#subset <- seq(from=1,to=n,by=100)
subset <- seq(n)
dat = list(N=N[subset],Y=Y[subset],n=length(subset))
nch <- 1
hlm = stan(model_name="Hierarchical Binomial Model", model_code = smod, data=dat , iter = nmc, chains = nch, verbose = TRUE, control=list(adapt_delta=0.9),algorithm="NUTS")
save.image(file='Outputs/hbm_hmc.RData')
} else {
  load('Outputs/hbm_hmc.RData')
}

tmp <- ggs(hlm)
gc()
df <- data.frame(tmp)
rm(tmp)
acs<-with(df,tapply(value,list(as.factor(Parameter),as.factor(Chain)),calcac))
gc()
acu <- unlist(acs)
acu <- t(matrix(acu,21,nch*length(unique(df$Parameter))))

ess <- with(df,tapply(value,list(as.factor(Parameter),as.factor(Chain)),calces))
ess <- c(ess)

AR <- data.frame(acu)
names(AR) <- paste('lag-',as.character(c(1,seq(from=5,to=100,by=5))),sep='')
AR$site <- seq(dim(AR)[1])

ES.HMC <- ess
ESscl <- ES.HMC/5000
ESscl <- data.frame(ESscl)
names(ESscl) <- c('EffectiveSize')

save(PHI.HMC,AR,ESscl,df,file='Outputs/HMC_maxpt.RData')

# acs <- acs/4
# ess<-
# 
# 
# #coef <- extract(hlm)
# #MU.HMC <- coef$mu[seq(from=1,to=20000,by=4),]
# PHI.HMC <- matrix(0,length(subset),21)
# ES.HMC <- matrix(0,length(subset),1)
# for (j in 1:length(subset)) {
#   print(j)
#   #f <- ar(MU.HMC[,j],order.max=1,aic=FALSE)
#   fs <- matrix(0,21,4)
#   ess <- matrix(0,4,1)
#   for (l in 1:4) {
#     mcmcobj <- as.mcmc(df$value[(df$Chain==l) & (df$Parameter==paste('mu.',j,'.',sep=''))])
#     f <- autocorr(mcmcobj,lags = c(1,seq(from=5,to=100,by=5)))
#     fs[,l] <- f
#     ess[l] <- effectiveSize(mcmcobj)
#   }
#   PHI.HMC[j,] <- apply(fs,1,mean)
#   ES.HMC[j] <- mean(ess)
# }
# 
# AR <- data.frame(PHI.HMC)

ARl <- melt(AR,id='site')
names(ARl) <- c('site','lag','autocorrelation')
ARl$lag <- gsub("lag-","",ARl$lag)
ARl$lag <- as.factor(as.numeric(ARl$lag))

png('Figures/HMC_maxpt_acf_box_full.png',width=450,height=300)
ggplot(ARl,aes(y=autocorrelation,x=lag)) + geom_boxplot(outlier.shape=NA) + theme(text=element_text(size=24)) + scale_y_continuous(limits=c(-.25,1)) + ggtitle('Hamiltonian Monte Carlo') +
scale_x_discrete(breaks=c(1,seq(from=10,to=100,by=10)))
#ggplot(ARl,aes(y=autocorrelation,x=lag)) + geom_violin() + theme(axis.text.x=element_text(angle=90)) + ggtitle('HMC -- Autocorrelations') + scale_y_continuous(limits=c(-.25,1))
dev.off()

png('Figures/HMC_maxpt_es_ecdf_full.png',width=450,height=300)
#ggplot(ESscl,aes(x=EffectiveSize)) + geom_histogram() + scale_x_continuous(limits=c(0,1.5)) + ggtitle('HMC -- Effective size (proportion of samples)')
ggplot(ESscl,aes(x=EffectiveSize)) + stat_ecdf(geom="step",size=1.25) + theme(text=element_text(size=24),axis.title.y=element_blank()) +  ggtitle('Hamiltonian Monte Carlo')
dev.off()

png('Figures/HMC_maxpt_es_hist_full.png',width=450,height=300)
ggplot(ESscl,aes(x=EffectiveSize)) + geom_histogram() + scale_x_continuous(limits=c(-0.1,1.5)) + theme(text=element_text(size=24),axis.title.y=element_blank()) +  ggtitle('Hamiltonian Monte Carlo')
dev.off()

# png('Figures/HMC_maxpt_hist.png',width=450,height=300)
# hist(PHI.HMC,30,col='gray',main='HMC',xlim=c(-.1,1))
#dev.off()

#datmh <- readMat('Outputs/binomhmmh_simple_small.mat')
for (j in 1:12) {
  dattmp <- readMat(paste('Outputs/binomhmmh_mu',j,'.mat',sep=''))
  if (j==1) {
    MU.MH <- dattmp$MU[5001:10000,]
  } else {
    MU.MH <- cbind(MU.MH,dattmp$MU[5001:10000,])
  }
}

#datmh <- readMat('Outputs/binomhmmh_mu.mat')
#PHI.MH <- datmh$AR
#MU.MH <- datmh$MUsml[5001:10000,]
#nsml <- 594
gc()
nsml <- dim(MU.MH)[2]

PHI.MH <- matrix(0,nsml,21)
ES.MH <- matrix(0,nsml,1)
for (j in 1:nsml) {
  if (j%%100==0) {
    print(j)
  }
  #f <- ar(MU.HMC[,j],order.max=1,aic=FALSE)
  f <- autocorr(as.mcmc(MU.MH[,j]),lags = c(1,seq(from=5,to=100,by=5)))
  PHI.MH[j,] <- c(f)
  es <- effectiveSize(MU.MH[,j])
  ES.MH[j] <- es
}
AR <- data.frame(PHI.MH)
names(AR) <- paste('lag-',as.character(c(1,seq(from=5,to=100,by=5))),sep='')
AR$site <- seq(length(subset))
ARl <- melt(AR,id='site')
names(ARl) <- c('site','lag','autocorrelation')
ARl$lag <- gsub("lag-","",ARl$lag)
ARl$lag <- as.factor(as.numeric(ARl$lag))

png('Figures/MH_maxpt_acf_box.png',width=450,height=300)
ggplot(ARl,aes(y=autocorrelation,x=lag)) + geom_boxplot(outlier.shape=NA) + theme(text=element_text(size=24)) + scale_y_continuous(limits=c(-.25,1)) + ggtitle('Metropolis') +
scale_x_discrete(breaks=c(1,seq(from=10,to=100,by=10)))
dev.off()
ESscl <- ES.MH/5000
ESscl <- data.frame(ESscl)
names(ESscl) <- c('EffectiveSize')
png('Figures/MH_maxpt_es_ecdf.png',width=450,height=300)
ggplot(ESscl,aes(x=EffectiveSize)) + stat_ecdf(geom="step",size=1.25) + scale_x_continuous(limits=c(0,1.5)) + ggtitle('Metropolis') + theme(text=element_text(size=24),axis.title.y=element_blank())
#ggplot(ESscl,aes(x=EffectiveSize)) + geom_histogram() + scale_x_continuous(limits=c(0,1)) + ggtitle('MH -- Effective size (proportion of samples)')
dev.off()


png('Figures/MH_maxpt_es_hist.png',width=450,height=300)
ggplot(ESscl,aes(x=EffectiveSize)) + geom_histogram() + scale_x_continuous(limits=c(-0.1,1.5)) + ggtitle('Metropolis') + theme(text=element_text(size=24),axis.title.y=element_blank())
#ggplot(ESscl,aes(x=EffectiveSize)) + geom_histogram() + scale_x_continuous(limits=c(0,1)) + ggtitle('MH -- Effective size (proportion of samples)')
dev.off()


# png('Figures/MH_maxpt_hist.png',width=450,height=300)
# hist(PHI.MH,50,col='gray',main='MH',xlim=c(-.1,1))
# dev.off()

# png('Figures/MH_maxpt_traces.png',width=600,height=600)
# par(mfrow=c(2,2),oma=rep(1.5,4),mar=rep(1.5,4))
# fake.ggplot(MU.MH[5000:10000,5],cex=.5,ylab='',xlab='t')
# fake.ggplot(MU.MH[5000:10000,47],cex=.5,ylab='',xlab='t')
# fake.ggplot(MU.MH[5000:10000,151],cex=.5,ylab='',xlab='t')
# fake.ggplot(MU.MH[5000:10000,387],cex=.5,ylab='',xlab='t')
# dev.off()
# 
# TAU <- datmh$TAU
# MU0 <- datmh$MU0
# save(PHI.MH,MU.MH,TAU,MU0,file='Outputs/MH_maxpt.RData')





