library(tidyverse)
library(splines)
library(nimble)
library(nimbleEcology)
library(coda)

folderpath = "~FitNmix/Data/ECN/"
filepath1 = "ProcessedBlueTitEcn.txt"

# The plan here is to fit the GAM n mixture models. 
# For model selection, I will first try GlobalGamNmix, SiteSpecificShrinkGamNmix, SiteSpecificFreeGamNmix 
# then for the best model I will try one knot lower and higher than the Fewster recommendation, which is N years * 0.3
# model selection will be based on WAIC 

# !Note! we only 'profile' the models here to save computation time, and we, therefore, only run one chain. For a proper analysis we would need more samples from the posterior
# and to spend more time evaluating chain mixing across 4 or so chains


# load data 
brdt = read.table(paste0(folderpath,filepath1))

# make spline matrix
years =length(unique(brdt$year))
numknots = round(years*0.3) # years * 0.3 Fewster recommendation 
knotsl=quantile(1:years,probs = seq(0,1,length.out=numknots))
knotsl[-c(1,numknots)] # gets rid of first and last vals
sp = bs(1:years,knots=knotsl[-c(1,numknots)],degree=3,intercept = F) # cubic
plot(NULL,xlim=c(1,years),ylim=c(0,1),ylab="basis",xlab="year")
for(i in 1:ncol(sp)){lines(1:years,sp[,i])}


# build list
SimListData1 = list(
  count = as.matrix(brdt[,11:16]),
  splineMat = sp[1:length(unique(brdt$year)),1:ncol(sp)]
)

# set up data 
SimCon1 = list(
  n = nrow(brdt),
  nsite = length(unique(brdt$SITECODE)),
  nyear = years,
  visits = brdt$visits,
  site = as.integer(as.factor(brdt$SITECODE)),
  year =  brdt$year - min(brdt$year) + 1,
  nspline = ncol(sp)
)

# initial values for parameters
Inits1 = list(
  N = apply(brdt[,11:16], 1, max), 
  beta1 = rep(0.01,ncol(sp)), 
  beta0 = rep(0.01,length(unique(brdt$SITECODE))),
  delta0 = 0.01
  
)


globalgammod <- nimbleCode({
  
  # Priors
  for(j in 1:nsite) {             
    beta0[j] ~ dnorm(2,sd=1)  #  Intercepts on ecology
  }
  
  for(s in 1:nspline) {
    beta1[s] ~ dnorm(0, sd=0.05)  # slope on ecology
  }
  
  delta0 ~ dnorm(0,sd=1.6) # intercept on observation
  
  yearweights[1:nyear] <- splineMat[1:nyear,1:nspline] %*%  beta1[1:nspline]
  logit(p) <- delta0  # this is for the observation model (on the logit scale)
  
  # Ecological model for true abundance
  for(i in 1:n) {
    log(lambda[i]) <- beta0[site[i]] + sum(yearweights[1:year[i]])
    N[i] ~ dpois(lambda[i])
    
    # Observation model for replicated counts - still only 1 p
    count[i,1:visits[i]] ~ dNmixture_s(lambda=lambda[i],p=p,Nmin=-1,Nmax=-1,len=visits[i])
    
  }
  
  # Derived quantities
  for(j in 1:nsite){        
    for(k in 1:nyear){    
      N.pred[j, k] <- exp(beta0[j] + sum(yearweights[1:k]))
    }
  }
})



# convert to nimble code
globalgam = nimbleModel(globalgammod, constants = SimCon1,data=SimListData1,inits = Inits1 )


watchlist =c("beta0","beta1","delta0","p","N.pred")

mcmcconf=configureMCMC(globalgam, enableWAIC = TRUE) 
mcmcconf$addMonitors(watchlist)
mcmcbuild = buildMCMC(mcmcconf) 
cglobalgam = compileNimble(globalgam) 
cglobalgamMCMC1 = compileNimble(mcmcbuild, project = globalgam) 
globalsamples = runMCMC(cglobalgamMCMC1,  niter = 22000,
                        nburnin = 2000,
                        thin = 10,
                        nchains=1,
                        samplesAsCodaMCMC = T, WAIC = TRUE)


G1 = summary(globalsamples$samples)
G1 = G1$statistics
G1 = rownames_to_column(data.frame(G1))
Nestimates = filter(G1, grepl("N",rowname))

Nestimates$siteind = rep(1:11,21)
Nestimates$yearind = rep(1:21,each=11)

# match the data
resultdtw = brdt
resultdtw$siteind = as.integer(as.factor(resultdtw$SITECODE))
resultdtw$yearind =  resultdtw$year - min(resultdtw$year) + 1

Nestimates = dplyr::select(Nestimates,Mean,siteind,yearind)
resultdtw=left_join(resultdtw,Nestimates,by=c("siteind","yearind"))

resultdtw=rename(resultdtw, globalmean =Mean)

ggplot(resultdtw,aes(year,meanind))+
  geom_point()+
  geom_point(aes(year,globalmean),color="red")+
  geom_line()+
  facet_wrap(~SITECODE)+
  theme_classic()


##
# 2) Site level GAM with shrinkage to global fit and a site level intercept 
# here we only have to change inits

# initial values for parameters
Inits2 = list(
  N = apply(brdt[,11:16], 1, max), 
  beta1 = matrix(0.01,nrow = ncol(sp),ncol=length(unique(brdt$SITECODE))), 
  beta0 = rep(0.01,length(unique(brdt$SITECODE))),
  delta0 = 0.01,
  hyperbeta = rep(0.01,ncol(sp)),
  sigmabeta =rep(0.01,ncol(sp))
  
)


shrinkgammod <- nimbleCode({
  
  # Priors
  for(j in 1:nsite) {             
    beta0[j] ~ dnorm(2,sd=1)  #  Intercepts on ecology
  }
  
  # hyper parameters for beta1
  for(s in 1:nspline) {
    hyperbeta[s] ~ dnorm(0,sd=0.05)
    sigmabeta[s] ~ dexp(5)
    for(l in 1:nsite) {
      beta1[s,l] ~ dnorm(hyperbeta[s], sigmabeta[s])  # slope on ecology
    }
  }
  
  delta0 ~ dnorm(0,sd=1.6) # intercept on observation
  logit(p) <- delta0  # this is for the observation model (on the logit scale)
  
  
  # define the coef matrix 
  for(m in 1:nsite){
    yearweights[m,1:nyear] <- splineMat[1:nyear,1:nspline] %*%  beta1[1:nspline,m]
  }
  
  # Ecological model for true abundance
  for(i in 1:n) {
    log(lambda[i]) <- beta0[site[i]] + sum(yearweights[site[i],1:year[i]])
    
    # Observation model for replicated count
    count[i,1:visits[i]] ~ dNmixture_s(lambda=lambda[i],p=p,Nmin=-1,Nmax=-1,len=visits[i])
    
  }
  
  # Derived quantities
  for(j in 1:nsite){        
    for(k in 1:nyear){    
      N.pred[j, k] <- exp(beta0[j] + sum(yearweights[j,1:k]))
    }
  }
})


# convert to nimble code
shrinkgam = nimbleModel(shrinkgammod, constants = SimCon1,data=SimListData1,inits = Inits2 )

mcmcconf2=configureMCMC(shrinkgam, enableWAIC = TRUE) 
mcmcconf2$addMonitors(watchlist)
mcmcbuild2 = buildMCMC(mcmcconf2) 
cshrinkgam = compileNimble(shrinkgam) 
cshrinkgamMCMC = compileNimble(mcmcbuild2, project = shrinkgam) 
shrinkgamsamples = runMCMC(cshrinkgamMCMC,  niter = 22000,
                           nburnin = 2000,
                           thin = 10,
                           nchains=1,
                           samplesAsCodaMCMC = T, WAIC = TRUE)


G2 = summary(shrinkgamsamples$samples)
G2 = G2$statistics
G2 = rownames_to_column(data.frame(G2))
# extract the N estimates 
Nestimates = filter(G2, grepl("N",rowname))

Nestimates$siteind = rep(1:11,21)
Nestimates$yearind = rep(1:21,each=11)


Nestimates = dplyr::select(Nestimates,Mean,siteind,yearind)
resultdtw=left_join(resultdtw,Nestimates,by=c("siteind","yearind"))

resultdtw=rename(resultdtw, shrinkmean =Mean)


ggplot(resultdtw,aes(year,meanind))+
  geom_point()+
  geom_point(aes(year,globalmean),color="red")+
  geom_point(aes(year,shrinkmean),color="blue")+
  geom_line()+
  facet_wrap(~SITECODE)+
  theme_classic()



# 3) Site level GAM with no shrinkage and site level intercept
Inits3 = list(
  N = apply(brdt[,11:16], 1, max), 
  beta1 = matrix(0.01,nrow = ncol(sp),ncol=length(unique(brdt$SITECODE))), # 6 knots so 7 splines
  beta0 = rep(0.01,length(unique(brdt$SITECODE))),
  delta0 = 0.01
)


freegammod <- nimbleCode({
  
  # Priors
  for(j in 1:nsite) {             
    beta0[j] ~ dnorm(2,sd=1)  #  Intercepts on ecology
  }
  
  # hyper parameters for beta1
  for(s in 1:nspline) {
    for(l in 1:nsite) {
      beta1[s,l] ~ dnorm(0,sd=0.05)  # slope on ecology
    }
  }
  
  delta0 ~ dnorm(0,sd=1.6) # intercept on observation
  logit(p) <- delta0  # this is for the observation model (on the logit scale)
  
  
  # define the coef matrix 
  for(m in 1:nsite){
    yearweights[m,1:nyear] <- splineMat[1:nyear,1:nspline] %*%  beta1[1:nspline,m]
  }
  
  # Ecological model for true abundance
  for(i in 1:n) {
    log(lambda[i]) <- beta0[site[i]] + sum(yearweights[site[i],1:year[i]])
    
    # Observation model for replicated count
    count[i,1:visits[i]] ~ dNmixture_s(lambda=lambda[i],p=p,Nmin=-1,Nmax=-1,len=visits[i])
    
  }
  
  # Derived quantities
  for(j in 1:nsite){      
    for(k in 1:nyear){
      N.pred[j, k] <- exp(beta0[j] + sum(yearweights[j,1:k]))
    }
  }
})

# convert to nimble code
freegam = nimbleModel(freegammod, constants = SimCon1,data=SimListData1,inits = Inits3 )

mcmcconf3=configureMCMC(freegam, enableWAIC = TRUE) 
mcmcconf3$addMonitors(watchlist)
mcmcbuild3 = buildMCMC(mcmcconf3)
cfreegam = compileNimble(freegam) 
cfreegamMCMC3 = compileNimble(mcmcbuild3, project = freegam) 
freegamsamples = runMCMC(cfreegamMCMC3 ,  niter = 22000,
                         nburnin = 2000,
                         thin = 10,
                         nchains=1,
                         samplesAsCodaMCMC = T,WAIC = TRUE)


G3 = summary(freegamsamples$samples)
G3 = G3$statistics
G3 = rownames_to_column(data.frame(G3))
# extract the N estimates 
Nestimates = filter(G3, grepl("N",rowname))

Nestimates$siteind = rep(1:11,21)
Nestimates$yearind = rep(1:21,each=11)


Nestimates = dplyr::select(Nestimates,Mean,siteind,yearind)
resultdtw=left_join(resultdtw,Nestimates,by=c("siteind","yearind"))

resultdtw=rename(resultdtw, freemean =Mean)

ggplot(resultdtw,aes(year,meanind,color="Mean"))+
  geom_point(size=2,alpha=0.8)+
  geom_point(aes(year,globalmean,color="Global"),size=2,alpha=0.8)+
  geom_point(aes(year,shrinkmean,color="Shrink"),size=2,alpha=0.8)+
  geom_point(aes(year,freemean,color="Free"),size=2,alpha=0.8)+
  facet_wrap(~SITECODE)+
  theme_classic()  +
  theme(axis.text = element_text(size=14),axis.title = element_text(size=18))+
  ylab("Population estimate")+
  xlab("Year")+
  scale_color_manual(name="Model",values=c(Mean = 'black',Global='red', Shrink='blue', Free='darkblue'))


# compare 
globalsamples$WAIC
shrinkgamsamples$WAIC 
freegamsamples$WAIC

effectiveSize(shrinkgamsamples$samples) # low effective sample size though - need to do more chains and more samples!



# 4) Try other knot numbers - using the shrinkage gam 
## try other knots

# one lower
dev.off()

numknots = 6
knotsl=quantile(1:years,probs = seq(0,1,length.out=numknots))
knotsl[-c(1,numknots)] # gets rid of first and last vals
sp6 = bs(1:years,knots=knotsl[-c(1,numknots)],degree=3,intercept = F) # cubic
plot(NULL,xlim=c(1,years),ylim=c(0,1),ylab="basis",xlab="year")
for(i in 1:ncol(sp6)){lines(1:years,sp6[,i])}



SimListData1low = list(
  count = as.matrix(brdt[,11:16]),
  splineMat = sp6[1:length(unique(brdt$year)),1:ncol(sp6)]
)
# set up data 
SimCon1low = list(
  n = nrow(brdt),
  nsite = length(unique(brdt$SITECODE)),
  nyear = years,
  visits = brdt$visits,
  site = as.integer(as.factor(brdt$SITECODE)),
  year =  brdt$year - min(brdt$year) + 1,
  nspline = ncol(sp6)
)

Inits2low = list(
  N = apply(brdt[,11:16], 1, max), 
  beta1 = matrix(0.01,nrow = ncol(sp6),ncol=length(unique(brdt$SITECODE))), 
  beta0 = rep(0,length(unique(brdt$SITECODE))),
  delta0 = 0,
  hyperbeta = rep(0,ncol(sp6)),
  sigmabeta =rep(0,ncol(sp6))
  
)

shrinkgamlow = nimbleModel(shrinkgammod, constants = SimCon1low,data=SimListData1low,inits = Inits2low )

mcmcconflow=configureMCMC(shrinkgamlow, enableWAIC = TRUE) 
mcmcconflow$addMonitors(watchlist)
mcmcbuildlow = buildMCMC(mcmcconflow) 
cshrinkgamlow = compileNimble(shrinkgamlow) 
cshrinkgamlowMCMC = compileNimble(mcmcbuildlow, project = shrinkgamlow) 
shrinkgamlowsamples = runMCMC(cshrinkgamlowMCMC,  niter = 22000,
                              nburnin = 2000,
                              thin = 10,
                              nchains=1,
                              samplesAsCodaMCMC = T, WAIC = TRUE)



G2low = summary(shrinkgamlowsamples$samples)
G2low = G2low$statistics
G2low = rownames_to_column(data.frame(G2low))
# extract the N estimates 
Nestimates = filter(G2low, grepl("N",rowname))

Nestimates$siteind = rep(1:11,21)
Nestimates$yearind = rep(1:21,each=11)

Nestimates = dplyr::select(Nestimates,Mean,siteind,yearind)
resultdtw=left_join(resultdtw,Nestimates,by=c("siteind","yearind"))

resultdtw=rename(resultdtw, shrinkmeanlow =Mean)
# add on new 


ggplot(resultdtw,aes(year,meanind))+
  geom_point()+
  geom_point(aes(year,shrinkmean),color="blue")+
  geom_point(aes(year,shrinkmeanlow),color="lightblue")+
  geom_line()+
  facet_wrap(~SITECODE)+
  theme_classic()


# one higher
numknots = 8
knotsl=quantile(1:years,probs = seq(0,1,length.out=numknots))
knotsl[-c(1,numknots)] # gets rid of first and last vals
sp8 = bs(1:years,knots=knotsl[-c(1,numknots)],degree=3,intercept = F) # cubic
plot(NULL,xlim=c(1,years),ylim=c(0,1),ylab="basis",xlab="year")
for(i in 1:ncol(sp8)){lines(1:years,sp8[,i])}


SimListData1high = list(
  count = as.matrix(brdt[,11:16]),
  splineMat = sp8[1:length(unique(brdt$year)),1:ncol(sp8)]
)
# set up data 
SimCon1high = list(
  n = nrow(brdt),
  nsite = length(unique(brdt$SITECODE)),
  nyear = years,
  visits = brdt$visits,
  site = as.integer(as.factor(brdt$SITECODE)),
  year =  brdt$year - min(brdt$year) + 1,
  nspline = ncol(sp8)
)

Inits2high = list(
  N = apply(brdt[,11:16], 1, max), 
  beta1 = matrix(0.01,nrow = ncol(sp8),ncol=length(unique(brdt$SITECODE))), 
  beta0 = rep(0,length(unique(brdt$SITECODE))),
  delta0 = 0,
  hyperbeta = rep(0,ncol(sp8)),
  sigmabeta =rep(0,ncol(sp8))
  
)

shrinkgamhigh = nimbleModel(shrinkgammod, constants = SimCon1high,data=SimListData1high,inits = Inits2high )

mcmcconfhigh=configureMCMC(shrinkgamhigh, enableWAIC = TRUE) 
mcmcconfhigh$addMonitors(watchlist)
mcmcbuildhigh = buildMCMC(mcmcconfhigh) 
cshrinkgamhigh = compileNimble(shrinkgamhigh) 
cshrinkgamhighMCMC = compileNimble(mcmcbuildhigh, project = shrinkgamhigh) 
shrinkgamhighsamples = runMCMC(cshrinkgamhighMCMC,  niter = 22000,
                               nburnin = 2000,
                               thin = 10,
                               nchains=1,
                               samplesAsCodaMCMC = T, WAIC = TRUE)



G2high = summary(shrinkgamhighsamples$samples)
G2high = G2high$statistics
G2high = rownames_to_column(data.frame(G2high))
# extract the N estimates 
Nestimates = filter(G2high, grepl("N",rowname))

Nestimates$siteind = rep(1:11,21)
Nestimates$yearind = rep(1:21,each=11)

Nestimates = dplyr::select(Nestimates,Mean,siteind,yearind)
resultdtw=left_join(resultdtw,Nestimates,by=c("siteind","yearind"))

resultdtw=rename(resultdtw, shrinkmeanhigh =Mean)

# add on new 

ggplot(resultdtw,aes(year,meanind))+
  geom_point()+
  geom_point(aes(year,shrinkmean),color="blue")+
  geom_point(aes(year,shrinkmeanlow),color="lightblue")+
  geom_point(aes(year,shrinkmeanhigh),color="darkblue")+
  geom_line()+
  facet_wrap(~SITECODE)+
  theme_classic()


# check waic - the high does better but not by much 
shrinkgamlowsamples$WAIC
shrinkgamsamples$WAIC
shrinkgamhighsamples$WAIC # for birds in this sample seems like 8 is better 


# this will produce a ton of plots - some of these seem fine, others look sketchy
# I think there is a need to up sampling and thinning 
plot(shrinkgamsamples$samples)

