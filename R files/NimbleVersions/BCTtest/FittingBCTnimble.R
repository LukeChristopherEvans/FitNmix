library(tidyverse)
library(splines)
library(nimble)
library(nimbleEcology)
library(coda)

folderpath = "/home/lukee/Insync/vs917256@reading.ac.uk/OneDrive Biz/FitNmix/Data/BCT/"
filepath1 = "FieldCommonPipWide.txt"

# The plan here is to fit the model across 100 sites
# A change from the previous model is that we will have detection changed per site 
# and we have parameters on the detector model 

# !Note! we only 'profile' the models here to save computation time, and we, therefore, only run one chain. For a proper analysis we would need more samples from the posterior
# and to spend more time evaluating chain mixing across 4 or so chains


# load data 
batdt = read.table(paste0(folderpath,filepath1))

# sample 100 sites 
set.seed(1234)
sitesample = sample(unique(batdt$SiteCode),100)
btdt=filter(batdt,SiteCode %in% sitesample)

# make the spline matrix 
years =length(unique(btdt$CountYear))
numknots = round(years*0.3) 
knotsl=quantile(1:years,probs = seq(0,1,length.out=numknots))
knotsl[-c(1,numknots)] # gets rid of first and last vals
sp = bs(1:years,knots=knotsl[-c(1,numknots)],degree=3,intercept = F) # cubic

# 1) Global GAM with site level intercept 
# build list
SimListData1 = list(
  count = as.matrix(btdt[,8:11]),
  spotsmatrix = as.matrix(btdt[,17:20]),
  splineMat = sp[1:length(unique(btdt$CountYear)),1:ncol(sp)],
  detectormatrix = as.matrix(btdt[,13:16])
)
# set up data 
SimCon1 = list(
  n = nrow(btdt),
  nsite = length(unique(btdt$SiteCode)),
  nyear = length(unique(btdt$CountYear)),
  visits = btdt$visitcount,
  site = as.integer(as.factor(btdt$SiteCode)),
  year =  btdt$CountYear - min(btdt$CountYear) + 1,
  nspline = ncol(sp),
  maxvisits =4,
  ndetectors =  length(1:max(btdt[,13:16]))
)

# initial values for parameters
Inits1 = list(
  N = apply(btdt[,8:11], 1, max)+1, 
  beta1 = rep(0.01,ncol(sp)), 
  beta0 = rep(0,length(unique(btdt$SiteCode))),
  delta0 = 0,
  delta1 = rep(0.01,length(1:max(btdt[,13:16]))),
  delta2 = 0.01
)



globalgammod <- nimbleCode({
  
  # Priors
  for(j in 1:nsite) {             
    beta0[j] ~ dnorm(2,1)  #  Intercepts on ecology
  }
  
  for(s in 1:nspline) {
    beta1[s] ~ dnorm(0, 0.05)  # slope on ecology
  }

  for(m in 1:ndetectors) {
    delta1[m] ~ dnorm(0, 1.6)  # slope on observation
  }
  
  delta0 ~ dnorm(0,1.6) # intercept on observation
  delta2 ~ dnorm(0,1.6)  
  
  yearweights[1:nyear] <- splineMat[1:nyear,1:nspline] %*%  beta1[1:nspline]

  # Ecological model for true abundance
  for(i in 1:n) {
    logit(p[i,1:visits[i]]) <- delta0 + delta2 * spotsmatrix[i,1:visits[i]] 
    log(lambda[i]) <- beta0[site[i]] + sum(yearweights[1:year[i]])
    N[i] ~ dpois(lambda[i])
    
    # Observation model for replicated counts - note we use dNmixture_v not s as we have vector of p's
    count[i,1:visits[i]] ~ dNmixture_v(lambda=N[i],p=p[i,1:visits[i]],Nmin=-1,Nmax=-1,len=visits[i])
    
  }
  

  # Derived quantities
  for(j in 1:nsite){        
    for(k in 1:nyear){    
      N.pred[j, k] <- exp( beta0[j] + sum(yearweights[1:k]) )
    }
  }
})



# convert to nimble code
globalgam = nimbleModel(globalgammod, constants = SimCon1,data=SimListData1,inits = Inits1 )

# I drop p here because it is ragged 
watchlist =c("beta0","beta1","delta0","delta1","delta2","N.pred") # delta 1 not fit here

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

Nestimates$siteind = rep(1:100,24)
Nestimates$yearind = rep(1:24,each=100)

# match the data
resultdtw = btdt
resultdtw$siteind = as.integer(as.factor(resultdtw$SiteCode))
resultdtw$yearind =  resultdtw$CountYear - min(resultdtw$CountYear) + 1

Nestimates = dplyr::select(Nestimates,Mean,SD,siteind,yearind)
resultdtw=left_join(resultdtw,Nestimates,by=c("siteind","yearind"))

resultdtw=rename(resultdtw, globalmean =Mean,globalsd=SD)

# get average and sum for a comparable non-modeled index 
meanind = function(x) mean(x[x!=-1])
sumind = function(x) sum(x[x!=-1])
avg=c()
sums=c()
for(i in 1:nrow(resultdtw)) {
  avg=c(avg,meanind(btdt[i,8:11]))
  sums=c(sums,sumind(btdt[i,8:11]))
}

resultdtw$avg= avg
resultdtw$sum= sums

filter(resultdtw, SiteCode %in% sample(unique(resultdtw$SiteCode),20)) %>%
ggplot(aes(CountYear,avg))+
  geom_point()+
  geom_point(aes(CountYear,globalmean),color="red")+
  geom_point(aes(CountYear,globalmean+globalsd),color="darkred")+
  geom_point(aes(CountYear,globalmean-globalsd),color="darkred")+
  facet_wrap(~SiteCode)+
  theme_classic()


# Now I can look at the observational model 
Destimates = filter(G1, grepl("delta",rowname))

plot(Destimates$Mean,1:nrow(Destimates),xlab=c("parameter estimate"),ylab="",col="red",pch=19,xlim = c(min(Destimates$Mean-Destimates$SD),max(Destimates$Mean+Destimates$SD)+1)  )
arrows(x0=Destimates$Mean-Destimates$SD , y0=1:nrow(Destimates), x1=Destimates$Mean+Destimates$SD, y1=1:nrow(Destimates), code=3, angle=90, length=0.1)
abline(v=0,lwd=1.2)
text(label=Destimates$rowname,x=rep(max(Destimates$Mean+Destimates$SD)+1,nrow(Destimates)),y=1:nrow(Destimates),cex=0.5 )



# We can also do detectors  - much slower when fitting this because of the second loop
# it's possible there is a faster way to write the model 
globalgammod1 <- nimbleCode({
  
  # Priors
  for(j in 1:nsite) {             
    beta0[j] ~ dnorm(2,1)  #  Intercepts on ecology
  }
  
  for(s in 1:nspline) {
    beta1[s] ~ dnorm(0, 0.05)  # slope on ecology
  }
  
  for(m in 1:ndetectors) {
    delta1[m] ~ dnorm(0, 1.6)  # slope on observation
  }
  
  delta0 ~ dnorm(0,1.6) # intercept on observation
  delta2 ~ dnorm(0,1.6)  
  
  yearweights[1:nyear] <- splineMat[1:nyear,1:nspline] %*%  beta1[1:nspline]
  
  # Ecological model for true abundance
  for(i in 1:n) {
    for(j in 1:visits[i]){ # there may be a faster way to do this
    logit(p[i,j]) <- delta0 + delta2 * spotsmatrix[i,j] + delta1[detectormatrix[i, j]]
    }
    log(lambda[i]) <- beta0[site[i]] + sum(yearweights[1:year[i]])
    N[i] ~ dpois(lambda[i])
    
    
    # Observation model for replicated counts 
    count[i,1:visits[i]] ~ dNmixture_v(lambda=N[i],p=p[i,1:visits[i]],Nmin=-1,Nmax=-1,len=visits[i])
    
  }
  
  
  # Derived quantities
  for(j in 1:nsite){        
    for(k in 1:nyear){    
      N.pred[j, k] <- exp( beta0[j] + sum(yearweights[1:k]) )
    }
  }
})



# convert to nimble code
globalgam1 = nimbleModel(globalgammod1, constants = SimCon1,data=SimListData1,inits = Inits1 )


mcmcconf1=configureMCMC(globalgam1, enableWAIC = TRUE) 
mcmcconf1$addMonitors(watchlist)
mcmcbuild1 = buildMCMC(mcmcconf1) 
cglobalgam1 = compileNimble(globalgam1) 
cglobalgamMCMC11 = compileNimble(mcmcbuild1, project = globalgam1) 
globalsamples1 = runMCMC(cglobalgamMCMC11,  niter = 22000,
                        nburnin = 2000,
                        thin = 10,
                        nchains=1,
                        samplesAsCodaMCMC = T, WAIC = TRUE)



G1.1 = summary(globalsamples1$samples)
G1.1 = G1.1$statistics
G1.1 = rownames_to_column(data.frame(G1.1))
Nestimates = filter(G1.1, grepl("N",rowname))

Nestimates$siteind = rep(1:100,24)
Nestimates$yearind = rep(1:24,each=100)

# match the data
Nestimates = dplyr::select(Nestimates,Mean,siteind,yearind)
resultdtw=left_join(resultdtw,Nestimates,by=c("siteind","yearind"))

resultdtw=rename(resultdtw, globalmean1 =Mean,globalsd1=SD)


filter(resultdtw, SiteCode %in% sample(unique(resultdtw$SiteCode),20)) %>%
  ggplot(aes(CountYear,avg))+
  geom_point()+
  geom_point(aes(CountYear,globalmean1),color="red")+
  geom_point(aes(CountYear,globalmean1+globalsd1),color="darkred")+
  geom_point(aes(CountYear,globalmean1-globalsd1),color="darkred")+
  facet_wrap(~SiteCode)+
  theme_classic()


# Now I can look at the observational model 
# shows that the number of sections completed increases bat detections (makes sense)
# and variation in the bat detectors
Destimates1 = filter(G1.1, grepl("delta",rowname))

plot(Destimates1$Mean,1:nrow(Destimates1),xlab=c("parameter estimate"),ylab="",col="red",pch=19,xlim = c(min(Destimates1$Mean-Destimates1$SD),max(Destimates1$Mean+Destimates1$SD)+1)  )
arrows(x0=Destimates1$Mean-Destimates1$SD , y0=1:nrow(Destimates1), x1=Destimates1$Mean+Destimates1$SD, y1=1:nrow(Destimates1), code=3, angle=90, length=0.1)
abline(v=0,lwd=1.2)
text(label=Destimates1$rowname,x=rep(max(Destimates1$Mean+Destimates1$SD)+1,nrow(Destimates1)),y=1:nrow(Destimates1),cex=0.5 )

# what are the detectors - note not all are fit in the sample
dts = 1:max(btdt[,13:16])
notin = dts[!(dts %in% as.matrix(btdt[,13:16]))]

# blues aren't fit in this sample
points(Destimates1$Mean[notin+1],(1:nrow(Destimates1))[notin+1],col="blue",pch=19 )

# we can also check the names of these detectors
# load data 
filepath2 = "Detectorcode.txt"
detectcode = read.table(paste0(folderpath,filepath2))

# you line up the id by the delta 1 number 
# but you should probably take this with a pinch of salt as it is only on the global model and with 100 sites



# 2) Site level spline with shrinkage to global fit and a site level intercept 


# initial values for parameters

Inits2 = list(
  N = apply(btdt[,8:11], 1, max)+1, 
  beta1 = matrix(0.01,nrow=ncol(sp),ncol=length(unique(btdt$SiteCode))), 
  beta0 = rep(0,length(unique(btdt$SiteCode))),
  
  hyperbeta = rep(0.01,ncol(sp)),
  sigmabeta =rep(0.01,ncol(sp)),
  
  delta0 = 0,
  delta1 = rep(0.01,length(1:max(btdt[,13:16]))),
  delta2 = 0.01
)


nmixsplineshrink <- nimbleCode({
  
  # Priors
  for(j in 1:nsite) {             
    beta0[j] ~ dnorm(2,1)  #  Intercepts on ecology
  }
  
  # hyper parameters for beta1
  for(s in 1:nspline) {
    hyperbeta[s] ~ dnorm(0,0.05)
    sigmabeta[s] ~ dexp(5)
    for(l in 1:nsite) {
      beta1[s,l] ~ dnorm(hyperbeta[s], sigmabeta[s])  # slope on ecology
    }
  }
  
  delta0 ~ dnorm(0,1.6) # intercept on observation
  delta2 ~ dnorm(0,1.6)  
  for(m in 1:ndetectors) {
    delta1[m] ~ dnorm(0, 1.6)  # slope on observation
  }
  
  # define the coef matrix 
  for(m in 1:nsite){
    yearweights[m,1:nyear] <- splineMat[1:nyear,1:nspline] %*%  beta1[1:nspline,m]
  }
  
  # Ecological model for true abundance
  for(i in 1:n) {
    for(j in 1:visits[i]){ 
      logit(p[i,j]) <- delta0 + delta2 * spotsmatrix[i,j] + delta1[detectormatrix[i, j]]
    }
    log(lambda[i]) <- beta0[site[i]] + sum(yearweights[site[i],1:year[i]])
    N[i] ~ dpois(lambda[i])
    
    
    # Observation model for replicated counts 
    count[i,1:visits[i]] ~ dNmixture_v(lambda=N[i],p=p[i,1:visits[i]],Nmin=-1,Nmax=-1,len=visits[i])
    
  }
  
  # Derived quantities
  for(j in 1:nsite){        
    for(k in 1:nyear){    
      N.pred[j, k] <- exp(beta0[j] + sum(yearweights[j,1:k]))
    }
  }
})


shrinkgam = nimbleModel(nmixsplineshrink, constants = SimCon1,data=SimListData1,inits = Inits2 )

mcmcconf2=configureMCMC(shrinkgam, enableWAIC = TRUE) 
mcmcconf2$addMonitors(c("beta0","beta1","delta0" ,"delta1","delta2", "N.pred"))
mcmcbuild2 = buildMCMC(mcmcconf2) 
cshrinkgam = compileNimble(shrinkgam) 
cshrinkgamMCMC1 = compileNimble(mcmcbuild2, project = shrinkgam) 
shrinksamples = runMCMC(cshrinkgamMCMC1,  niter = 22000,
                        nburnin = 2000,
                        thin = 10,
                        nchains=1,
                        samplesAsCodaMCMC = T, WAIC = TRUE)



G2 = summary(shrinksamples$samples)
G2 = G2$statistics
G2 = rownames_to_column(data.frame(G2))
Nestimates = filter(G2, grepl("N",rowname))

Nestimates$siteind = rep(1:100,24)
Nestimates$yearind = rep(1:24,each=100)

Nestimates = dplyr::select(Nestimates,Mean,SD,siteind,yearind)
resultdtw=left_join(resultdtw,Nestimates,by=c("siteind","yearind"))

resultdtw=rename(resultdtw, shrinkmean =Mean,shrinksd=SD)

# plot
filter(resultdtw, SiteCode %in% sample(unique(resultdtw$SiteCode),20)) %>%
  ggplot(aes(CountYear,avg))+
  geom_point()+
  geom_point(aes(CountYear,globalmean),color="red")+
  geom_point(aes(CountYear,shrinkmean),color="blue")+
  facet_wrap(~SiteCode)+
  theme_classic()


## 3) Site level GAM with no shrinkage and site level intercept
Inits3 = list(
  N = apply(btdt[,8:11], 1, max)+1, 
  beta1 = matrix(0.01,nrow=ncol(sp),ncol=length(unique(btdt$SiteCode))), 
  beta0 = rep(0,length(unique(btdt$SiteCode))),
  delta0 = 0,
  delta1 = rep(0.01,length(1:max(btdt[,13:16]))),
  delta2 = 0.01
)



nmixsplinefree <- nimbleCode({
  
  # Priors
  for(j in 1:nsite) {             
    beta0[j] ~ dnorm(2,1)  #  Intercepts on ecology
  }
  
  # hyper parameters for beta1
  for(s in 1:nspline) {
    for(l in 1:nsite) {
      beta1[s,l] ~ dnorm(0, 0.05)  # slope on ecology
    }
  }
  
  delta0 ~ dnorm(0,1.6) # intercept on observation
  delta2 ~ dnorm(0,1.6)  
  for(m in 1:ndetectors) {
    delta1[m] ~ dnorm(0, 1.6)  # slope on observation
  }
  
  # define the coef matrix 
  for(m in 1:nsite){
    yearweights[m,1:nyear] <- splineMat[1:nyear,1:nspline] %*%  beta1[1:nspline,m]
  }
  
  # Ecological model for true abundance
  for(i in 1:n) {
    for(j in 1:visits[i]){ 
      logit(p[i,j]) <- delta0 + delta2 * spotsmatrix[i,j] + delta1[detectormatrix[i, j]]
    }
    log(lambda[i]) <- beta0[site[i]] + sum(yearweights[site[i],1:year[i]])
    N[i] ~ dpois(lambda[i])
    
    
    # Observation model for replicated counts 
    count[i,1:visits[i]] ~ dNmixture_v(lambda=N[i],p=p[i,1:visits[i]],Nmin=-1,Nmax=-1,len=visits[i])
    
  }
  
  # Derived quantities
  for(j in 1:nsite){        
    for(k in 1:nyear){    
      N.pred[j, k] <- exp(beta0[j] + sum(yearweights[j,1:k]))
    }
  }
})

freegam = nimbleModel(nmixsplinefree, constants = SimCon1,data=SimListData1,inits = Inits3 )

mcmcconf3=configureMCMC(freegam, enableWAIC = TRUE) 
mcmcconf3$addMonitors(watchlist)
mcmcbuild3 = buildMCMC(mcmcconf3) 
cfreegam = compileNimble(freegam) 
cfreegamMCMC1 = compileNimble(mcmcbuild3, project = freegam) 
freesamples = runMCMC(cfreegamMCMC1,  niter = 22000,
                        nburnin = 2000,
                        thin = 10,
                        nchains=1,
                        samplesAsCodaMCMC = T, WAIC = TRUE)


G3 = summary(freesamples$samples)
G3 = G3$statistics
G3 = rownames_to_column(data.frame(G3))
Nestimates = filter(G3, grepl("N",rowname))

Nestimates$siteind = rep(1:100,24)
Nestimates$yearind = rep(1:24,each=100)

Nestimates = dplyr::select(Nestimates,Mean,SD,siteind,yearind)
resultdtw=left_join(resultdtw,Nestimates,by=c("siteind","yearind"))

resultdtw=rename(resultdtw, freemean =Mean,freesd=SD)

# plot
filter(resultdtw, SiteCode %in% sample(unique(resultdtw$SiteCode),20)) %>%
  ggplot(aes(CountYear,avg,color="Mean"))+
  geom_point(size=2,alpha=0.8)+
  geom_point(aes(CountYear,globalmean,color="Global"),size=2,alpha=0.8)+
  geom_point(aes(CountYear,shrinkmean,color="Shrink"),size=2,alpha=0.8)+
  geom_point(aes(CountYear,freemean,color="Free"),size=2,alpha=0.8)+
  facet_wrap(~SiteCode)+
  theme_classic() +
  theme(axis.text = element_text(size=14),axis.title = element_text(size=18))+
  ylab("Population estimate")+
  xlab("Year")+
  scale_color_manual(name="Model",values=c(Mean = 'black',Global='red', Shrink='blue', Free='darkblue'))




# free is best so we can look at the detectors 
Destimates3 = filter(G3, grepl("delta",rowname))

plot(Destimates3$Mean,1:nrow(Destimates3),xlab=c("Parameter estimate"),ylab="",col="red",pch=19,xlim = c(min(Destimates3$Mean-Destimates3$SD),max(Destimates3$Mean+Destimates3$SD)+1),cex=1.2,cex.lab=1.2  )
arrows(x0=Destimates3$Mean-Destimates3$SD , y0=1:nrow(Destimates3), x1=Destimates3$Mean+Destimates3$SD, y1=1:nrow(Destimates3), code=3, angle=90, length=0.1)
abline(v=0,lwd=1.2)
text(label=Destimates3$rowname,x=rep(max(Destimates3$Mean+Destimates3$SD)+1,nrow(Destimates3)),y=1:nrow(Destimates3),cex=0.5 )

# blues aren't fit in this sample
points(Destimates3$Mean[notin+1],(1:nrow(Destimates3))[notin+1],col="blue",pch=19,cex=1.2 )


# compare models 
globalsamples$WAIC;shrinksamples$WAIC ;freesamples$WAIC # free samples seems best 


# 4) Try other knot numbers - using the best model
## try other knots
numknots = 6
knotsl=quantile(1:years,probs = seq(0,1,length.out=numknots))
knotsl[-c(1,numknots)] # gets rid of first and last vals
sp6 = bs(1:years,knots=knotsl[-c(1,numknots)],degree=3,intercept = F) # cubic
plot(NULL,xlim=c(1,years),ylim=c(0,1),ylab="basis",xlab="year")
for(i in 1:ncol(sp6)){lines(1:years,sp6[,i])}

SimListData6 = list(
  count = as.matrix(btdt[,8:11]),
  spotsmatrix = as.matrix(btdt[,17:20]),
  splineMat = sp6[1:length(unique(btdt$CountYear)),1:ncol(sp6)],
  detectormatrix = as.matrix(btdt[,13:16])
)
# set up data 
SimCon6 = list(
  n = nrow(btdt),
  nsite = length(unique(btdt$SiteCode)),
  nyear = length(unique(btdt$CountYear)),
  visits = btdt$visitcount,
  site = as.integer(as.factor(btdt$SiteCode)),
  year =  btdt$CountYear - min(btdt$CountYear) + 1,
  nspline = ncol(sp6),
  maxvisits =4,
  ndetectors =  length(1:max(btdt[,13:16]))
)

Inits6 = list(
  N = apply(btdt[,8:11], 1, max)+1, 
  beta1 = matrix(0.01,nrow=ncol(sp6),ncol=length(unique(btdt$SiteCode))), 
  beta0 = rep(0,length(unique(btdt$SiteCode))),
  delta0 = 0.01,
  delta1 = rep(0.01,length(1:max(btdt[,13:16]))),
  delta2 = 0.01
)

freegam6 = nimbleModel(nmixsplinefree, constants = SimCon6,data=SimListData6,inits = Inits6 )

mcmcconf3.6=configureMCMC(freegam6, enableWAIC = TRUE) 
mcmcconf3.6$addMonitors(watchlist)
mcmcbuild3.6 = buildMCMC(mcmcconf3.6) 
cfreegam6 = compileNimble(freegam6) 
cfreegam6MCMC1 = compileNimble(mcmcbuild3.6, project = freegam6) 
freesamples6 = runMCMC(cfreegam6MCMC1,  niter = 22000,
                      nburnin = 2000,
                      thin = 10,
                      nchains=1,
                      samplesAsCodaMCMC = T, WAIC = TRUE)


# one up
numknots = 8
knotsl=quantile(1:years,probs = seq(0,1,length.out=numknots))
knotsl[-c(1,numknots)] # gets rid of first and last vals
sp8 = bs(1:years,knots=knotsl[-c(1,numknots)],degree=3,intercept = F) # cubic
plot(NULL,xlim=c(1,years),ylim=c(0,1),ylab="basis",xlab="year")
for(i in 1:ncol(sp8)){lines(1:years,sp8[,i])}

SimListData8 = list(
  count = as.matrix(btdt[,8:11]),
  spotsmatrix = as.matrix(btdt[,17:20]),
  splineMat = sp8[1:length(unique(btdt$CountYear)),1:ncol(sp8)],
  detectormatrix = as.matrix(btdt[,13:16])
)
# set up data 
SimCon8 = list(
  n = nrow(btdt),
  nsite = length(unique(btdt$SiteCode)),
  nyear = length(unique(btdt$CountYear)),
  visits = btdt$visitcount,
  site = as.integer(as.factor(btdt$SiteCode)),
  year =  btdt$CountYear - min(btdt$CountYear) + 1,
  nspline = ncol(sp8),
  maxvisits =4,
  ndetectors =  length(1:max(btdt[,13:16]))
)

Inits8 = list(
  N = apply(btdt[,8:11], 1, max)+1, 
  beta1 = matrix(0.01,nrow=ncol(sp8),ncol=length(unique(btdt$SiteCode))), 
  beta0 = rep(0,length(unique(btdt$SiteCode))),
  delta0 = 0.01,
  delta1 = rep(0.01,length(1:max(btdt[,13:16]))),
  delta2 = 0.01
)

freegam8 = nimbleModel(nmixsplinefree, constants = SimCon8,data=SimListData8,inits = Inits8 )

mcmcconf3.8=configureMCMC(freegam8, enableWAIC = TRUE) 
mcmcconf3.8$addMonitors(watchlist)
mcmcbuild3.8 = buildMCMC(mcmcconf3.8) 
cfreegam8 = compileNimble(freegam8) 
cfreegam8MCMC1 = compileNimble(mcmcbuild3.8, project = freegam8) 
freesamples8 = runMCMC(cfreegam8MCMC1,  niter = 22000,
                       nburnin = 2000,
                       thin = 10,
                       nchains=1,
                       samplesAsCodaMCMC = T, WAIC = TRUE)


# Gam with 8 knots preferred 
freesamples6$WAIC ; freesamples$WAIC ; freesamples8$WAIC
