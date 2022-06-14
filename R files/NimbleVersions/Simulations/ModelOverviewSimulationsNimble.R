### Model fit simulation ###
# The goal here is to simulate very simple count data from a process similar to BTC and RSPB 
# counts and then fit and explore basic N mixture models using nimble

# !Note! we only 'profile' the models here to save computation time, and we, therefore, only run one chain. For a proper analysis we would need more samples from the posterior
# and to spend more time evaluating chain mixing across 4 or so chains

library(gstat)
library(tidyverse)
library(nimble)
library(coda)

## Here we test nimble versions of the n-mixture model

# Model choices:
# global linear trend with site level intercept


## Stage 1: set up naive simulation and run the base model

folderpath = "/home/lukee/Insync/vs917256@reading.ac.uk/OneDrive Biz/FitNmix"

# create landscape with spatial autocorrelation 
area =  expand.grid(1:100,1:100)
names(area) <- c("x","y")
# spatial correlation
a.corr = gstat(formula=z~1, locations = ~x+y,beta = 1,dummy =T, model=vgm(psill = 20,range=10,model='Exp'),nmax=20)


# generate basemap; positive integers only
set.seed(1234)
amap = predict(a.corr,newdata=area,nsim=1)
amap$sim1=amap$sim1 + abs(min(amap$sim1))

# plot
ggplot(amap,aes(x=x,y=y,fill=sim1))+
  geom_raster()

#generate data from samples 
xsamp<-round(runif(20,1,100))
ysamp<-round(runif(20,1,100))

ggplot(amap,aes(x=x,y=y,fill=sim1))+
  geom_raster()+
  geom_point(data=data.frame(xsamp,ysamp),aes(xsamp,ysamp),inherit.aes = F,color="red")

# use a linear increasing trend 
years =20
sitedata = matrix(NA,nrow = length(xsamp),ncol = years)
for(k in 1:ncol(sitedata)) { 
  for(i in 1:nrow(sitedata)) {
    sitedata[i,k] = amap[amap$x==xsamp[i] & amap$y==ysamp[i],]$sim1 + k * 0.2
  }
}

#sitedata row site, column sample 
plot(1:20,sitedata[1,]) # linear increasing trend 

# now sample from this multiple times using poisson error 
# and given detection probability 

p = 0.8 # detection prob


observdata = data.frame()
for(j in 1:ncol(sitedata)) { # each year
  for(i in 1:nrow(sitedata)) {
    visits = round(runif(1,min=1,max=5)) # between 1 and 5 visits to the site 
    counts = round(rpois(visits,sitedata[i,j]) * p) # counts 
    td=data.frame(site=i,year=j+1990,count=counts,x=xsamp[i],y=ysamp[i])
    observdata=rbind(observdata,td)
    
  }
}



# counts across sites - trend somewhat obscured by sampling variation  (also apparent population growth and declines coming solely from sampling error!!)
ggplot(observdata,aes(year,count))+
  geom_point(size=2)+
  facet_wrap(~site)+
  theme_classic()+
  theme(axis.text = element_text(size=14),axis.title = element_text(size=18))+
  ylab("Count")+
  xlab("Year")


# wide format for sending to nimble
observdata=observdata %>% group_by(site,year) %>% mutate(visit=row_number())
observdataw=observdata %>% pivot_wider(names_from=visit,values_from= count)
observdataw=replace_na(observdataw,list('1'=-1,'2'=-1,'3'=-1,'4'=-1,'5'=-1)) # pad with -1
# now its row
observdataw$visits=apply(observdataw[,5:9] >=0,1,sum)


# In nimble we keep the model within R 

# build list
SimListData1 = list(
  count = as.matrix(observdataw[,5:9]),
  year = observdataw$year - min(observdataw$year) + 1
)
# set up data 
SimCon1 = list(
  n = nrow(observdataw),
  nsite = length(unique(observdataw$site)),
  nyear = length(unique(observdataw$year)),
  visits = observdataw$visits,
  site = observdataw$site
)

# initial values for parameters
Inits1 = list(
  N = apply(observdataw[,5:9], 1, max), 
  beta1 = 0,
  beta0 = rep(0,length(unique(observdataw$site))),
  delta0 = 0
)

# Best to use '<-' assignment in nimble
nmix1 <- nimbleCode({
  
  # Priors
  for(j in 1:nsite) {             
    beta0[j] ~ dnorm(2,1)  #  Intercepts on ecology
   }
   beta1 ~ dnorm(0, 0.05)  # slope on ecology
   delta0 ~ dnorm(0,1.6) # intercept on observation
  
    # Ecological model for true abundance
   for(i in 1:n) {
    log(lambda[i]) <- beta0[site[i]] + beta1 * year[i]
    N[i] ~ dpois(lambda[i])
  
    # Observation model for replicated counts
    for(k in 1:visits[i]) {
      logit(p) <- delta0  # this is for the observation model (on the logit scale)
      count[i,k] ~ dbin(p , N[i])
       }
  
    }
  
  # Derived quantities
  for(j in 1:nsite){        
    for(k in 1:nyear){    
      N.pred[j, k] <- exp( beta0[j] + beta1 * k)
    }
  }
})

# convert to nimble code
runnmix1 = nimbleModel(nmix1, constants = SimCon1,data=SimListData1,inits = Inits1 )


watchlist =c("beta0","beta1","delta0","p","N.pred")

mcmcconf=configureMCMC(runnmix1) 
mcmcconf$addMonitors(watchlist)
mcmcbuild = buildMCMC(mcmcconf) # build it 
crun = compileNimble(runnmix1) # you must compile the original model
crunMCMC = compileNimble(mcmcbuild, project = runnmix1) #note runnmix1  not crun! This seems strange to me!
samples = runMCMC(crunMCMC,  niter = 22000,
                              nburnin = 2000,
                              thin = 10,
                              nchains=1,
                              samplesAsCodaMCMC = T)


# stuff to evaluate
ab = summary(samples)
ab= ab$statistics
ab=rownames_to_column(data.frame(ab))
Nestimates = filter(ab, grepl("N",rowname))


# add on the prediction
resultdtw = observdataw # make a copy of original data 
resultdtw$Nest= Nestimates$Mean

# quick function for grabbing data from the ragged array and applying a function
finds=function(f,x,y){
  z=c()
  for(i in 1:y) {
    z =c(z, pull(x[1,i]))
  }
  f(z)
}

# get average and sum for a comparable non-modeled index 
avg=c()
sums=c()
for(i in 1:400) {
  avg=c(avg,finds(mean,resultdtw[i,5:9],resultdtw$visits[i]))
  sums=c(sums,finds(sum,resultdtw[i,5:9],resultdtw$visits[i]))
}

resultdtw$avg= avg
resultdtw$sum= sums

# have a look at the fit 
par(mfrow=c(2,2))
for (i in c(2,4,8,16)){
  plot(1:20,sitedata[i,],ylim=c(0,65),ylab=paste("pop size site",i),xlab="Year") # real data - site 8
  points(1:20,filter(resultdtw,site==i)$Nest,col="red") # N mixture
  points(1:20,filter(resultdtw,site==i)$avg,col="blue") # mean
  points(1:20,filter(resultdtw,site==i)$sum,col="green") # sum
}

# put this back
par(mfrow=c(1,1))

# chain diagnostics -
#plot(samples)
#autocorr.plot(samples)
#gelman.diag(samples)


## Stage 2: Features of the nimble model  

# Why does it overestimate abundance?

# We see above that the model is very effective at estimating 
# the trend but overestimates the abundance. This is because a classic nmixture model is
# set up in a way that we imagine we have a fixed population (not a density - an actual number of individuals) per site and per year 
# that we sample from and which cannot be smaller than the maximum observed value (in the range minpop:k).
# Effectively the lambda in the simulation is really more of a latent density and not a latent population size - causing the mismatch. Or another way of saying
# it, is that the lambda parameter shows what is the likely mean number I would observe in a sample given a fixed population in the local region.

# We can see this feature of the model here in an example with an actual fixed population size: 
# We imagine we have set sized populations, growing by one individual per year, observed while moving around randomly in a 1x1 grid (i.e. a poisson process).
# we observe half the square and so have the possibility of seeing half of these 
# animals on average (lambda is 0.5 * N). I also set the observation prob centered on 80%
# we'll complete the simulation at just two sites 

# Here we run the simulations to demonstrate the effect in nimble 

# heres the sampling scheme for 1 
Site1N = 10 # starting pop site 1
Site2N = 20 # starting pop site 2

# example sample at t=0
rx = runif(Site1N,0,1)
ry = runif(Site1N,0,1)

plot(rx,ry,xlim=c(0,1),ylim=c(0,1))
abline(v=0.5)
text(0.1,0.9,labels=paste("count =", sum(rx<0.5)))


# Now we run two multi visit simulations: the first with sampling half the square and the second 
# with sampling all the square - we will see that the model estimates the same population size (and not the density)
# because it is trying to estimate the overall fixed population
FixedSizeDataset1 = data.frame()

for(i in 1:years) {
  
  N1possibles  = c(sum(runif(Site1N,0,1) < 0.5),sum(runif(Site1N,0,1) < 0.5),sum(runif(Site1N,0,1) < 0.5)) #sample from half square
  N2possibles  = c(sum(runif(Site2N,0,1) < 0.5),sum(runif(Site2N,0,1) < 0.5),sum(runif(Site2N,0,1) < 0.5))
  
  N1observed = round(N1possibles *  rnorm(3,p,0.05)) # obs error
  N2observed = round(N2possibles *  rnorm(3,p,0.05))
  
  # build dataframe of obs
  thisyear= rbind(data.frame(site=1,year=i+1990,wholesquare =Site1N,halfsquare=Site1N*0.5, x1=N1observed[1],x2=N1observed[2],x3=N1observed[3]), c(site=2,year=i+1990,wholesquare =Site2N,halfsquare=Site2N*0.5,x1=N2observed[1],x2=N2observed[2],x3=N2observed[3]))
  FixedSizeDataset1 = rbind.data.frame(FixedSizeDataset1,thisyear)
  
  # grow both pops
  Site1N = Site1N + 1
  Site2N = Site2N + 1
  
}

#rest for simulation 2
Site1N = 10
Site2N = 20

FixedSizeDataset2 = data.frame()
for(i in 1:years) {
  
  N1possibles  = c(Site1N,Site1N,Site1N) # sample from the whole square
  N2possibles  = c(Site2N,Site2N,Site2N)
  
  N1observed = round(N1possibles * rnorm(3,p,0.05))
  N2observed = round(N2possibles * rnorm(3,p,0.05))
  
  thisyear= rbind(data.frame(site=1,year=i+1990,wholesquare =Site1N,halfsquare=Site1N*0.5, x1=N1observed[1],x2=N1observed[2],x3=N1observed[3]), c(site=2,year=i+1990,wholesquare =Site2N,halfsquare=Site2N*0.5,x1=N2observed[1],x2=N2observed[2],x3=N2observed[3]))
  FixedSizeDataset2 = rbind.data.frame(FixedSizeDataset2,thisyear)
  
  
  Site1N = Site1N + 1
  Site2N = Site2N + 1
  
}

##  we need to reset the data 
# build list
SimListData2.1 = list(
  count = as.matrix(FixedSizeDataset1[,5:7]),
  year = FixedSizeDataset1$year - min(FixedSizeDataset1$year) + 1
)
# set up data 
SimCon2.1 = list(
  n = nrow(FixedSizeDataset1),
  nsite = length(unique(FixedSizeDataset1$site)),
  nyear = length(unique(FixedSizeDataset1$year)),
  visits = rep(3,nrow(FixedSizeDataset1)),
  site = FixedSizeDataset1$site
)

# initial values for parameters
Inits2.1 = list(
  N = apply(FixedSizeDataset1[,5:7], 1, max), 
  beta1 = 0,
  beta0 = rep(0,length(unique(FixedSizeDataset1$site))),
  delta0 = 0
)

# second sim
SimListData2.2 = list(
  count = as.matrix(FixedSizeDataset2[,5:7]),
  year = FixedSizeDataset2$year - min(FixedSizeDataset2$year) + 1
)
# set up data 
SimCon2.2 = list(
  n = nrow(FixedSizeDataset2),
  nsite = length(unique(FixedSizeDataset2$site)),
  nyear = length(unique(FixedSizeDataset2$year)),
  visits = rep(3,nrow(FixedSizeDataset2)),
  site = FixedSizeDataset2$site
)

# initial values for parameters
Inits2.2 = list(
  N = apply(FixedSizeDataset2[,5:7], 1, max), 
  beta1 = 0,
  beta0 = rep(0,length(unique(FixedSizeDataset2$site))),
  delta0 = 0
)


watchlist =c("beta0","beta1","delta0","p","N.pred")

runnmix2.1 = nimbleModel(nmix1, constants = SimCon2.1,data=SimListData2.1,inits = Inits2.1 )

# mcmc part
mcmcconf2.1=configureMCMC(runnmix2.1) 
mcmcconf2.1$addMonitors(watchlist)
mcmcconf2.1$enableWAIC =T
mcmcbuild2.1 = buildMCMC(mcmcconf2.1) 
crun2.1 = compileNimble(runnmix2.1)
crunMCMC2.1 = compileNimble(mcmcbuild2.1, project = runnmix2.1) 
set.seed(1234)
sim2.1 = runMCMC(crunMCMC2.1,  niter = 22000,
                  nburnin = 2000,
                  thin = 10,
                  nchains=1,
                  samplesAsCodaMCMC = T)


runnmix2.2 = nimbleModel(nmix1, constants = SimCon2.2,data=SimListData2.2,inits = Inits2.2 )

# mcmc part
mcmcconf2.2=configureMCMC(runnmix2.2) 
mcmcconf2.2$addMonitors(watchlist)
mcmcbuild2.2 = buildMCMC(mcmcconf2.2) 
crun2.2 = compileNimble(runnmix2.2) 
crunMCMC2.2 = compileNimble(mcmcbuild2.2, project = runnmix2.2) 
set.seed(1234)
sim2.2 = runMCMC(crunMCMC2.2,  niter = 22000,
                 nburnin = 2000,
                 thin = 10,
                 nchains=1,
                 samplesAsCodaMCMC = T)


#get data 
ab2.1 = summary(sim2.1)
ab2.1 = ab2.1$statistics
ab2.1=rownames_to_column(data.frame(ab2.1))
Nestimates2.1 = filter(ab2.1, grepl("N",rowname))
FixedSizeDataset1$Nest= Nestimates2.1$Mean

ab2.2 = summary(sim2.2)
ab2.2 = ab2.2$statistics
ab2.2=rownames_to_column(data.frame(ab2.2))
Nestimates2.2 = filter(ab2.2, grepl("N",rowname))
FixedSizeDataset2$Nest= Nestimates2.2$Mean

# get average for an alternative index 
FixedSizeDataset1$avg = apply(FixedSizeDataset1[,5:7],1,mean)
FixedSizeDataset2$avg = apply(FixedSizeDataset2[,5:7],1,mean)


# Now look at the estimates - the key thing is that the model estimates the whole square 
# not the local density in the sample which should be 50% of the whole square 
par(mfrow=c(1,1))
plot(1:20,filter(FixedSizeDataset1,site==1)$wholesquare,ylim=c(0,50),ylab=paste("pop size site",1),xlab="Year") # data from the whole square
points(1:20,filter(FixedSizeDataset1,site==1)$halfsquare,col="grey") # half square (Expected lambda)
points(1:20,filter(FixedSizeDataset1,site==1)$Nest,col="red") # N mixture
points(1:20,filter(FixedSizeDataset1,site==1)$avg,col="blue") # mean

plot(1:20,filter(FixedSizeDataset1,site==2)$wholesquare,ylim=c(0,50),ylab=paste("pop size site",1),xlab="Year") # data from the whole square
points(1:20,filter(FixedSizeDataset1,site==2)$halfsquare,col="grey") # half square (Expected lambda)
points(1:20,filter(FixedSizeDataset1,site==2)$Nest,col="red") # N mixture
points(1:20,filter(FixedSizeDataset1,site==2)$avg,col="blue") # mean

# and because of this it puts the 50% sampling into the observation error (roughly 0.5*0.8)
plot(sim2.1[,45],main="p")

# here is the estimates when sampling the whole square (exactly the same N estimates)
plot(1:20,filter(FixedSizeDataset2,site==1)$wholesquare,ylim=c(0,50),ylab=paste("pop size site",1),xlab="Year") # data from the whole square
points(1:20,filter(FixedSizeDataset2,site==1)$Nest,col="red") # N mixture
points(1:20,filter(FixedSizeDataset2,site==1)$avg,col="blue") # mean

plot(1:20,filter(FixedSizeDataset2,site==2)$wholesquare,ylim=c(0,50),ylab=paste("pop size site",1),xlab="Year") # real data - site 8
points(1:20,filter(FixedSizeDataset2,site==2)$Nest,col="red") # N mixture
points(1:20,filter(FixedSizeDataset2,site==2)$avg,col="blue") # mean

# but now better estimate of the 0.8 observation error 
plot(sim2.2[,45],main="p")


# What to do?

# A second approach is to model the population as a density (not the finite pop size)
# and to have the marginalizing bounded between the observed count per visit respectively  
# I implement this in the nimble model below 

nmix2 <- nimbleCode({
  # Priors
  for(j in 1:nsite) {             
    beta0[j] ~ dnorm(2,1)  
  }
  beta1 ~ dnorm(0, 0.05)  
  delta0 ~ dnorm(0,1.6) 
  
  # Ecological model for true abundance
  for(i in 1:n) {
    for(k in 1:visits[i]) {
    log(lambda[i,k]) <- beta0[site[i]] + beta1 * year[i]
     N[i,k] ~ dpois(lambda[i,k])
    # Observation model for replicated counts
      logit(p) <- delta0  
      count[i,k] ~ dbin(p , N[i,k])
    }
    
  }
  
  # Derived quantities
  for(j in 1:nsite){        
    for(k in 1:nyear){    
      N.pred[j, k] <- exp( beta0[j] + beta1 * k)
    }
  }
})

# build list
SimListData3 = list(
  count = as.matrix(FixedSizeDataset1[,5:7]),
  year = FixedSizeDataset1$year - min(FixedSizeDataset1$year) + 1
)
# set up data 
SimCon3 = list(
  n = nrow(FixedSizeDataset1),
  nsite = length(unique(FixedSizeDataset1$site)),
  nyear = length(unique(FixedSizeDataset1$year)),
  visits = rep(3,nrow(FixedSizeDataset1)),
  site = FixedSizeDataset1$site
)

# initial values for parameters
Inits3 = list(
  N = as.matrix(FixedSizeDataset1[,5:7]), 
  beta1 = 0,
  beta0 = rep(0,length(unique(FixedSizeDataset1$site))),
  delta0 = 0
)

runnmixdens = nimbleModel(nmix2, constants = SimCon3,data=SimListData3,inits = Inits3 )

# mcmc part
mcmcconf3=configureMCMC(runnmixdens) 
mcmcconf3$addMonitors(watchlist)
mcmcbuild3 = buildMCMC(mcmcconf3) 
crun3 = compileNimble(runnmixdens) 
crunMCMC3 = compileNimble(mcmcbuild3, project = runnmixdens) 
set.seed(1234)
sim3 = runMCMC(crunMCMC3,   niter = 22000,
               nburnin = 2000,
               thin = 10,
               nchains=1,
               samplesAsCodaMCMC = T)



ab3 = summary(sim3)
ab3 = ab3$statistics
ab3=rownames_to_column(data.frame(ab3))
Nestimates3 = filter(ab3, grepl("N",rowname))
FixedSizeDataset1$Nestdens= Nestimates3$Mean

# now this estimates the lambda as a density more correctly 
plot(1:20,filter(FixedSizeDataset1,site==1)$wholesquare,ylim=c(0,50),ylab=paste("pop size site",1),xlab="Year") # data from the whole square
points(1:20,filter(FixedSizeDataset1,site==1)$halfsquare,col="grey") # half square (Expected lambda)
points(1:20,filter(FixedSizeDataset1,site==1)$Nest,col="red") # original N mixture
points(1:20,filter(FixedSizeDataset1,site==1)$avg,col="blue") # mean
points(1:20,filter(FixedSizeDataset1,site==1)$Nestdens,col="green") # new n mixture

plot(1:20,filter(FixedSizeDataset1,site==2)$wholesquare,ylim=c(0,50),ylab=paste("pop size site",1),xlab="Year") # data from the whole square
points(1:20,filter(FixedSizeDataset1,site==2)$halfsquare,col="grey") # half square (Expected lambda)
points(1:20,filter(FixedSizeDataset1,site==2)$Nest,col="red") # original N mixture
points(1:20,filter(FixedSizeDataset1,site==2)$avg,col="blue") # mean
points(1:20,filter(FixedSizeDataset1,site==2)$Nestdens,col="green") # new n mixture

# lots of uncertainty here but peaks over correct value
plot(sim3[,45],main="p")


# speed comparison 
dentime = system.time(runMCMC(crunMCMC3,  niter = 50000,
                                nburnin = 2000,
                                thin = 10,
                                nchains=1,
                                samplesAsCodaMCMC = T))
nondentime = system.time(runMCMC(crunMCMC2.1,  niter = 50000,
                                 nburnin = 2000,
                                 thin = 10,
                                 nchains=1,
                                 samplesAsCodaMCMC = T))

dentime ; nondentime

coda::effectiveSize(simeco)



# As we use nimble we can also use nimble ecology 
# that provides an implementation of 'The summation over N uses the efficient method given by Meehan et al. (2020)'
# Meehan, T. D., Michel, N. L., & Rue, H. (2020). Estimating Animal Abundance with N-Mixture
# Models Using the R—INLA Package for R. Journal of Statistical Software, 95(2). https://doi.org/10.18637/jss.v095.i02

# Basically this is a fast version of the marginalising I do in stan where I specify an upper K - it gives
# very similar results to stan, but should be quick!

library(nimbleEcology)

nmixeco <- nimbleCode({
  
  # Priors
  for(j in 1:nsite) {             
    beta0[j] ~ dnorm(2,1)  
  }
  beta1 ~ dnorm(0, 0.05)  
  delta0 ~ dnorm(0,1.6) 
  
  # Ecological model for true abundance
  for(i in 1:n) {
    log(lambda[i]) <- beta0[site[i]] + beta1 * year[i]
    N[i] ~ dpois(lambda[i])
    
    # Observation model for replicated counts
      logit(p) <- delta0  
      count[i,1:visits[i]] ~ dNmixture_s(lambda=N[i],p=p,Nmin=-1,Nmax=-1,len=visits[i])
    
    
  }
  
  # Derived quantities
  for(j in 1:nsite){        
    for(k in 1:nyear){    
      N.pred[j, k] <- exp( beta0[j] + beta1 * k)
    }
  }
})


runnmixeco = nimbleModel(nmixeco, constants = SimCon2.1,data=SimListData2.1,inits = Inits2.1 )

# mcmc part
mcmcconfeco=configureMCMC(runnmixeco, WAIC = TRUE) 
mcmcconfeco$addMonitors(watchlist)
mcmcbuildeco = buildMCMC(mcmcconfeco) 
cruneco = compileNimble(runnmixeco) 
crunMCMCeco = compileNimble(mcmcbuildeco, project = runnmixeco) 
set.seed(1234)
simeco = runMCMC(crunMCMCeco, niter = 22000,
                 nburnin = 2000,
                 thin = 10,
                 nchains=1,
                 samplesAsCodaMCMC = T)

abeco = summary(simeco)
abeco  = abeco$statistics
abeco=rownames_to_column(data.frame(abeco))
Nestimateseco = filter(abeco, grepl("N",rowname))
FixedSizeDataset1$Nestdenseco= Nestimateseco$Mean

plot(1:20,filter(FixedSizeDataset1,site==1)$wholesquare,ylim=c(0,50),ylab=paste("pop size site",1),xlab="Year") # data from the whole square
points(1:20,filter(FixedSizeDataset1,site==1)$halfsquare,col="grey") # half square (Expected lambda)
points(1:20,filter(FixedSizeDataset1,site==1)$Nest,col="red") # original N mixture
points(1:20,filter(FixedSizeDataset1,site==1)$avg,col="blue") # mean
points(1:20,filter(FixedSizeDataset1,site==1)$Nestdens,col="green") # dens n mixture
points(1:20,filter(FixedSizeDataset1,site==1)$Nestdenseco,col="pink") # new n mixture


plot(1:20,filter(FixedSizeDataset1,site==2)$wholesquare,ylim=c(0,50),ylab=paste("pop size site",1),xlab="Year") # data from the whole square
points(1:20,filter(FixedSizeDataset1,site==2)$halfsquare,col="grey") # half square (Expected lambda)
points(1:20,filter(FixedSizeDataset1,site==2)$Nest,col="red") # original N mixture
points(1:20,filter(FixedSizeDataset1,site==2)$avg,col="blue") # mean
points(1:20,filter(FixedSizeDataset1,site==2)$Nestdens,col="green") # dens n mixture
points(1:20,filter(FixedSizeDataset1,site==2)$Nestdenseco,col="pink") # new n mixture


# nice output for report 
par(mfrow=c(1,2))
for(i in c(1,2)){
  plot(1:20,filter(FixedSizeDataset1,site==i)$wholesquare,ylim=c(0,70),ylab=paste("Population size site",i),xlab="Year",pch=19,col="grey",cex=1.2,cex.lab=1.2) 
  points(1:20,filter(FixedSizeDataset1,site==i)$Nest,col="red",pch=19,cex=1.2) # N mixture
  points(1:20,filter(FixedSizeDataset1,site==i)$avg,col="black",pch=19,cex=1.2) # mean
  points(1:20,filter(FixedSizeDataset1,site==i)$Nestdenseco,col="pink",pch=19,cex=1.2) # sum
}
legend(1, 67, legend=c("True value", "Mean of counts","Nmixture estimate","Nmixture estimate approx"),
       fill =c("grey", "black","red","pink"), lty=1:2, cex=0.8)


dentime2 = system.time(runMCMC(crunMCMCeco,  niter = 22000,
                               nburnin = 2000,
                               thin = 10,
                               nchains=1,
                               samplesAsCodaMCMC = T))

# interestingly slower but probably better estimates
dentime ; dentime2

# check effective sizes
effectiveSize(sim2.1) ; effectiveSize(sim2.2) ; effectiveSize(sim3); effectiveSize(simeco)
