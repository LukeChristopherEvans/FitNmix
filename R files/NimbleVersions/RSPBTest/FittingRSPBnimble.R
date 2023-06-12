library(tidyverse)
library(splines)
library(nimble)
library(nimbleEcology)
library(coda)

#vignette("Introduction_to_nimbleEcology")

folderpath = "~FitNmix/Data/RSPB/"
filepath1 = "BlueTitWide.txt"

# !Note! we only 'profile' the models here to save computation time, and we, therefore, only run one chain. For a proper analysis we would need more samples from the posterior
# and to spend more time evaluating chain mixing across 4 or so chains

#! LUKE check whether we want N[i] ~ dpois(lambda[i]) dNmixture_v(lambda=N[i] ; or  dNmixture_v(lambda=lambda[i]

# load data 
birddt = read.table(paste0(folderpath,filepath1))

# sample 100 sites 
set.seed(1234)
sitesample = sample(unique(birddt$Gridref),100)
brdt=filter(birddt,Gridref %in% sitesample)

# make the spline matrix 
years =length(unique(brdt$Year))
numknots = round(years*0.3) 
knotsl=quantile(1:years,probs = seq(0,1,length.out=numknots))
knotsl[-c(1,numknots)] # gets rid of first and last vals
sp = bs(1:years,knots=knotsl[-c(1,numknots)],degree=3,intercept = F) # cubic


# build list
SimListData1 = list(
  count = as.matrix(brdt[,9:16]),
  visitmatrix = as.matrix(brdt[,18:25]),
  bandmatrix = as.matrix(brdt[,18:25]),
  splineMat = sp[1:length(unique(brdt$Year)),1:ncol(sp)]
)
# set up data 
SimCon1 = list(
  n = nrow(brdt),
  nsite = length(unique(brdt$Gridref)),
  nyear = length(unique(brdt$Year)),
  visits = brdt$visitcount,
  site = as.integer(as.factor(brdt$Gridref)),
  year =  brdt$Year - min(brdt$Year) + 1,
  nspline = ncol(sp),
  maxvisits =8
)

# initial values for parameters
Inits1 = list(
  N = apply(brdt[,9:16], 1, max)+1, 
  beta1 = rep(0.01,ncol(sp)), 
  beta0 = rep(0,length(unique(brdt$Gridref))),
  delta0 = 0,
  delta1 = 0.01,
  delta2 = 0.01
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
  delta1 ~ dnorm(0,sd=1.6) # slope on obs
  delta2 ~ dnorm(0,sd=1.6)  
  
  yearweights[1:nyear] <- splineMat[1:nyear,1:nspline] %*%  beta1[1:nspline]
  
  # Ecological model for true abundance
  for(i in 1:n) {
    logit(p[i,1:visits[i]]) <- delta0 + delta1 * visitmatrix[i,1:visits[i]]  + delta2 *  bandmatrix[i,1:visits[i]] 
    log(lambda[i]) <- beta0[site[i]] + sum(yearweights[1:year[i]])
    
    # Observation model for replicated counts 
    count[i,1:visits[i]] ~ dNmixture_v(lambda=lambda[i],p=p[i,1:visits[i]],Nmin=-1,Nmax=-1,len=visits[i])
    
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
watchlist =c("beta0","beta1","delta0","delta1","delta2","N.pred") # delta 1 not fit!

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

Nestimates$siteind = rep(1:100,27)
Nestimates$yearind = rep(1:27,each=100)

# match the data
resultdtw = brdt
resultdtw$siteind = as.integer(as.factor(resultdtw$Gridref))
resultdtw$yearind =  resultdtw$Year - min(resultdtw$Year) + 1

Nestimates = dplyr::select(Nestimates,Mean,SD,siteind,yearind)
resultdtw=left_join(resultdtw,Nestimates,by=c("siteind","yearind"))

resultdtw=rename(resultdtw, globalmean =Mean,globalsd=SD)

# get average and sum for a comparable non-modeled index 
meanind = function(x) mean(x[x!=-1])
sumind = function(x) sum(x[x!=-1])
avg=c()
sums=c()
for(i in 1:nrow(resultdtw)) {
  avg=c(avg,meanind(brdt[i,9:16]))
  sums=c(sums,sumind(brdt[i,9:16]))
}

resultdtw$avg= avg
resultdtw$sum= sums

filter(resultdtw, Gridref %in% sample(unique(resultdtw$Gridref),20)) %>%
  ggplot(aes(Year,avg))+
  geom_point()+
  geom_point(aes(Year,globalmean),color="red")+
  geom_point(aes(Year,globalmean+globalsd),color="darkred")+
  geom_point(aes(Year,globalmean-globalsd),color="darkred")+
  facet_wrap(~Gridref)+
  theme_classic()


# Now I can look at the observational model - increase distance band decreases detectability - visits no estimated effect
Destimates = filter(G1, grepl("delta",rowname))

plot(Destimates$Mean,1:nrow(Destimates),xlab=c("parameter estimate"),ylab="",col="red",pch=19,xlim = c(min(Destimates$Mean-Destimates$SD),max(Destimates$Mean+Destimates$SD)+1)  )
arrows(x0=Destimates$Mean-Destimates$SD , y0=1:nrow(Destimates), x1=Destimates$Mean+Destimates$SD, y1=1:nrow(Destimates), code=3, angle=90, length=0.1)
abline(v=0,lwd=1.2)
text(label=Destimates$rowname,x=rep(max(Destimates$Mean+Destimates$SD)+1,nrow(Destimates)),y=1:nrow(Destimates),cex=0.5 )


# 2) Site level spline with shrinkage to global fit and a site level intercept 


# initial values for parameters
Inits2 = list(
  N = apply(brdt[,9:16], 1, max)+1, 
  beta1 = matrix(0.01,nrow=ncol(sp),ncol=length(unique(brdt$Gridref))), 
  beta0 = rep(0.01,length(unique(brdt$Gridref))),
  
  hyperbeta = rep(0.01,ncol(sp)),
  sigmabeta =rep(0.01,ncol(sp)),
  
  delta0 = 0.01,
  delta1 = 0.01,
  delta2 = 0.01
)


nmixsplineshrink <- nimbleCode({
  
  # Priors
  for(j in 1:nsite) {             
    beta0[j] ~ dnorm(2,1)  #  Intercepts on ecology
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
  delta2 ~ dnorm(0,sd=1.6)  
  delta1 ~ dnorm(0,sd=1.6)

  
  # define the coef matrix 
  for(m in 1:nsite){
    yearweights[m,1:nyear] <- splineMat[1:nyear,1:nspline] %*%  beta1[1:nspline,m]
  }
  
  # Ecological model for true abundance
  for(i in 1:n) {
    logit(p[i,1:visits[i]]) <- delta0 + delta1 * visitmatrix[i,1:visits[i]]  + delta2 *  bandmatrix[i,1:visits[i]] 
    log(lambda[i]) <- beta0[site[i]] + sum(yearweights[site[i],1:year[i]])

    # Observation model for replicated counts 
    count[i,1:visits[i]] ~ dNmixture_v(lambda=lambda[i],p=p[i,1:visits[i]],Nmin=-1,Nmax=-1,len=visits[i])
    
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

Nestimates$siteind = rep(1:100,27)
Nestimates$yearind = rep(1:27,each=100)

Nestimates = dplyr::select(Nestimates,Mean,SD,siteind,yearind)
resultdtw=left_join(resultdtw,Nestimates,by=c("siteind","yearind"))

resultdtw=rename(resultdtw, shrinkmean =Mean,shrinksd=SD)

# plot
filter(resultdtw, Gridref %in% sample(unique(resultdtw$Gridref),20)) %>%
  ggplot(aes(Year,avg))+
  geom_point()+
  geom_point(aes(Year,globalmean),color="red")+
  geom_point(aes(Year,shrinkmean),color="blue")+
  facet_wrap(~Gridref)+
  theme_classic()


## 3) Site level GAM with no shrinkage and site level intercept

Inits3 = list(
  N = apply(brdt[,9:16], 1, max)+1, 
  beta1 = matrix(0.01,nrow=ncol(sp),ncol=length(unique(brdt$Gridref))), 
  beta0 = rep(0.01,length(unique(brdt$Gridref))),
  
  delta0 = 0.01,
  delta1 = 0.01,
  delta2 = 0.01
)


nmixsplinefree <- nimbleCode({
  
  # Priors
  for(j in 1:nsite) {             
    beta0[j] ~ dnorm(2,sd=1)  #  Intercepts on ecology
  }
  
  # hyper parameters for beta1
  for(s in 1:nspline) {
    for(l in 1:nsite) {
      beta1[s,l] ~ dnorm(0,sd= 0.05)  # slope on ecology
    }
  }
  
  delta0 ~ dnorm(0,sd=1.6) # intercept on observation
  delta2 ~ dnorm(0,sd=1.6)  
  delta1 ~ dnorm(0,sd=1.6)  # slope on observation
  
  
  # define the coef matrix 
  for(m in 1:nsite){
    yearweights[m,1:nyear] <- splineMat[1:nyear,1:nspline] %*%  beta1[1:nspline,m]
  }
  
  # Ecological model for true abundance
  for(i in 1:n) {
    logit(p[i,1:visits[i]]) <- delta0 + delta1 * visitmatrix[i,1:visits[i]]  + delta2 *  bandmatrix[i,1:visits[i]] 
    log(lambda[i]) <- beta0[site[i]] + sum(yearweights[site[i],1:year[i]])

    # Observation model for replicated counts 
    count[i,1:visits[i]] ~ dNmixture_v(lambda=lambda[i],p=p[i,1:visits[i]],Nmin=-1,Nmax=-1,len=visits[i])
    
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


Nestimates$siteind = rep(1:100,27)
Nestimates$yearind = rep(1:27,each=100)

Nestimates = dplyr::select(Nestimates,Mean,SD,siteind,yearind)
resultdtw=left_join(resultdtw,Nestimates,by=c("siteind","yearind"))

resultdtw=rename(resultdtw, freemean =Mean,freesd=SD)

# plot

filter(resultdtw, Gridref %in% sample(unique(resultdtw$Gridref),20)) %>%
  ggplot(aes(Year,avg,color="Mean"))+
  geom_point(size=2,alpha=0.8)+
  geom_point(aes(Year,globalmean,color="Global"),size=2,alpha=0.8)+
  geom_point(aes(Year,shrinkmean,color="Shrink"),size=2,alpha=0.8)+
  geom_point(aes(Year,freemean,color="Free"),size=2,alpha=0.8)+
  facet_wrap(~Gridref)+
  theme_classic() +
  theme(axis.text = element_text(size=14),axis.title = element_text(size=18))+
  ylab("Population estimate")+
  xlab("Year")+
  scale_color_manual(name="Model",values=c(Mean = 'black',Global='red', Shrink='blue', Free='darkblue'))




# compare models 
globalsamples$WAIC;shrinksamples$WAIC ;freesamples$WAIC # shrink samples seems to work best (but remeber chains and effective sample size caveats)



# detectability on best model - distance reduce detectability - possible higher detectability on late
Destimates2 = filter(G2, grepl("delta",rowname))

plot(Destimates2$Mean,1:nrow(Destimates2),xlab=c("Parameter estimate"),ylab="",col="red",pch=19,xlim = c(min(Destimates2$Mean-Destimates2$SD),max(Destimates2$Mean+Destimates2$SD)+1),cex=1.2,cex.lab=1.2   )
arrows(x0=Destimates2$Mean-Destimates2$SD , y0=1:nrow(Destimates2), x1=Destimates2$Mean+Destimates2$SD, y1=1:nrow(Destimates2), code=3, angle=90, length=0.1)
abline(v=0,lwd=1.2)
text(label=Destimates2$rowname,x=rep(max(Destimates2$Mean+Destimates2$SD)+1,nrow(Destimates2)),y=1:nrow(Destimates2),cex=0.5 )



# 4) Try other knot numbers - using the best model
## try other knots

# one lower
dev.off()

numknots = 6
knotsl=quantile(1:years,probs = seq(0,1,length.out=numknots))
knotsl[-c(1,numknots)] # gets rid of first and last vals
sp6 = bs(1:years,knots=knotsl[-c(1,numknots)],degree=3,intercept = F) # cubic
plot(NULL,xlim=c(1,years),ylim=c(0,1),ylab="basis",xlab="year")
for(i in 1:ncol(sp6)){lines(1:years,sp6[,i])}


SimListData6 = list(
  count = as.matrix(brdt[,9:16]),
  visitmatrix = as.matrix(brdt[,18:25]),
  bandmatrix = as.matrix(brdt[,18:25]),
  splineMat = sp6[1:length(unique(brdt$Year)),1:ncol(sp6)]
)
# set up data 
SimCon6 = list(
  n = nrow(brdt),
  nsite = length(unique(brdt$Gridref)),
  nyear = length(unique(brdt$Year)),
  visits = brdt$visitcount,
  site = as.integer(as.factor(brdt$Gridref)),
  year =  brdt$Year - min(brdt$Year) + 1,
  nspline = ncol(sp6),
  maxvisits =8
)

Inits6 = list(
  N = apply(brdt[,9:16], 1, max)+1, 
  beta1 = matrix(0.01,nrow=ncol(sp6),ncol=length(unique(brdt$Gridref))), 
  beta0 = rep(0.01,length(unique(brdt$Gridref))),
  
  hyperbeta = rep(0.01,ncol(sp6)),
  sigmabeta =rep(0.01,ncol(sp6)),
  
  delta0 = 0.01,
  delta1 = 0.01,
  delta2 = 0.01
)


shrinkgam6 = nimbleModel(nmixsplineshrink, constants = SimCon6,data=SimListData6,inits = Inits6 )

mcmcconf6=configureMCMC(shrinkgam6, enableWAIC = TRUE) 
mcmcconf6$addMonitors(c("beta0","beta1","delta0" ,"delta1","delta2", "N.pred"))
mcmcbuild6 = buildMCMC(mcmcconf6) 
cshrinkgam6 = compileNimble(shrinkgam6) 
cshrinkgamMCMC6 = compileNimble(mcmcbuild6, project = shrinkgam6) 
shrinksamples6 = runMCMC(cshrinkgamMCMC6,  niter = 22000,
                        nburnin = 2000,
                        thin = 10,
                        nchains=1,
                        samplesAsCodaMCMC = T, WAIC = TRUE)



# one up
SimListData8 = list(
  count = as.matrix(brdt[,9:16]),
  visitmatrix = as.matrix(brdt[,18:25]),
  bandmatrix = as.matrix(brdt[,18:25]),
  splineMat = sp8[1:length(unique(brdt$Year)),1:ncol(sp8)]
)
# set up data 
SimCon8 = list(
  n = nrow(brdt),
  nsite = length(unique(brdt$Gridref)),
  nyear = length(unique(brdt$Year)),
  visits = brdt$visitcount,
  site = as.integer(as.factor(brdt$Gridref)),
  year =  brdt$Year - min(brdt$Year) + 1,
  nspline = ncol(sp8),
  maxvisits =8
)

Inits8 = list(
  N = apply(brdt[,9:16], 1, max)+1, 
  beta1 = matrix(0.01,nrow=ncol(sp8),ncol=length(unique(brdt$Gridref))), 
  beta0 = rep(0.01,length(unique(brdt$Gridref))),
  
  hyperbeta = rep(0.01,ncol(sp8)),
  sigmabeta =rep(0.01,ncol(sp8)),
  
  delta0 = 0.01,
  delta1 = 0.01,
  delta2 = 0.01
)


shrinkgam8 = nimbleModel(nmixsplineshrink, constants = SimCon8,data=SimListData8,inits = Inits8 )

mcmcconf8=configureMCMC(shrinkgam8, enableWAIC = TRUE) 
mcmcconf8$addMonitors(c("beta0","beta1","delta0" ,"delta1","delta2", "N.pred"))
mcmcbuild8 = buildMCMC(mcmcconf8) 
cshrinkgam8 = compileNimble(shrinkgam8) 
cshrinkgamMCMC8 = compileNimble(mcmcbuild8, project = shrinkgam8) 
shrinksamples8 = runMCMC(cshrinkgamMCMC8,  niter = 22000,
                         nburnin = 2000,
                         thin = 10,
                         nchains=1,
                         samplesAsCodaMCMC = T, WAIC = TRUE)




shrinksamples6$WAIC;shrinksamples$WAIC ;shrinksamples8$WAIC 