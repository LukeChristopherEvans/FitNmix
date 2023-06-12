library(gstat)
library(tidyverse)
library(nimble)
library(nimbleEcology)

# Here we develop and explore GAM based fits as these are the 
# most realistic approach to modeling changes in population abundance 
# over time.
# 1) Global GAM with site level intercept 
# 2) Site level GAM with shrinkage to global fit and a site level intercept 
# 3) Site level GAM with no shrinkage and site level intercept
# 4) Site level GAM with geographic distance based shrinkage 


# create landscape with spatial autocorrelation (We use the spatial autocorrelation in model 4)
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


# we need introduce a wiggle in the population
xwig <- seq(0, pi * 2, length.out=20)
sinwig <- sin(xwig)

plot(1:20,sinwig*5) # change in numbers 
years = 20

# simulate samples 
sitedatawiggle = matrix(NA,nrow = length(xsamp),ncol = years)
for(k in 1:ncol(sitedatawiggle)) { 
  for(i in 1:nrow(sitedatawiggle)) {
    sitedatawiggle[i,k] = amap[amap$x==xsamp[i] & amap$y==ysamp[i],]$sim1 + (sinwig*5)[k]
  }
}

#sitedata row= site, column = sample 
plot(1:20,sitedatawiggle[13,]) #
min(sitedatawiggle) # check the min is above 0 (no negative pop sizes)

# run the sampling procedure
p =0.8
observdatawiggle = data.frame()
for(j in 1:ncol(sitedatawiggle)) { # each year
  for(i in 1:nrow(sitedatawiggle)) {
    visits = round(runif(1,min=1,max=5))
    counts = round(rpois(visits,sitedatawiggle[i,j]) * p)
    td=data.frame(site=i,year=j+1990,count=counts,x=xsamp[i],y=ysamp[i])
    observdatawiggle=rbind(observdatawiggle,td)
    
  }
}


# To model the 'wiggle' we need to use splines 
library(splines)
# code adapted from McElreath rethinking
numknots = 6 # start with 6 knots
knotsl=quantile(1:20,probs = seq(0,1,length.out=numknots))
knotsl[-c(1,numknots)] # gets rid of first and last vals

sp = bs(1:years,knots=knotsl[-c(1,numknots)],degree=3,intercept = F) # cubic
plot(NULL,xlim=c(1,20),ylim=c(0,1),ylab="basis",xlab="year")
for(i in 1:ncol(sp)){lines(1:20,sp[,i])}


# going to get max count per season as the minimum of the theoretical
# population size that year - this goes into the stan model 
minpopswig = observdatawiggle %>% group_by(site,year) %>%
  summarise(minpop = max(count))
observdatawiggle = left_join(observdatawiggle,minpopswig,by=c("site","year"))

# wide format 
observdatawiggle=observdatawiggle %>% group_by(site,year) %>% mutate(visit=row_number())
# these names are getting out of control
observdataWideWiggle=observdatawiggle %>% pivot_wider(names_from=visit,values_from= count)
observdataWideWiggle=replace_na(observdataWideWiggle,list('1'=-1,'2'=-1,'3'=-1,'4'=-1,'5'=-1)) # pad with -1
# now its row
observdataWideWiggle$visits=apply(observdataWideWiggle[,6:10] >=0,1,sum)

ggplot(observdatawigglexy,aes(year,count))+
  geom_point(size=2)+
  facet_wrap(~site)+
  theme_classic() +
  theme(axis.text = element_text(size=14),axis.title = element_text(size=18))+
  ylab("Count")+
  xlab("Year")



# 1) Global GAM with site level intercept 
# build list
SimListData1 = list(
  count = as.matrix(observdataWideWiggle[,6:10]),
  splineMat = sp[1:length(unique(observdataWideWiggle$year)),1:ncol(sp)]
)

# set up data 
SimCon1 = list(
  n = nrow(observdataWideWiggle),
  nsite = length(unique(observdataWideWiggle$site)),
  nyear = length(unique(observdataWideWiggle$year)),
  visits = observdataWideWiggle$visits,
  site = observdataWideWiggle$site,
  year = observdataWideWiggle$year - min(observdataWideWiggle$year) + 1,
  nspline = 7
)

# initial values for parameters
Inits1 = list(
  N = apply(observdataWideWiggle[,6:10], 1, max), 
  beta1 = rep(0.01,7), # 6 knots so 7 splines
  beta0 = rep(0,length(unique(observdataWideWiggle$site))),
  delta0 = 0

)


nmixspline <- nimbleCode({
  
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
nmixspline1 = nimbleModel(nmixspline, constants = SimCon1,data=SimListData1,inits = Inits1 )


watchlist =c("beta0","beta1","delta0","p","N.pred")

mcmcconf=configureMCMC(nmixspline1) 
mcmcconf$addMonitors(watchlist)
mcmcbuild = buildMCMC(mcmcconf) 
crun1 = compileNimble(nmixspline1) 
crunMCMC1 = compileNimble(mcmcbuild, project = nmixspline1) 
samples = runMCMC(crunMCMC1,  niter = 22000,
                  nburnin = 2000,
                  thin = 10,
                  nchains=1,
                  samplesAsCodaMCMC = T)


# stuff to evaluate
G1 = summary(samples)
G1 = G1$statistics
G1 = rownames_to_column(data.frame(G1))
Nestimates = filter(G1, grepl("N",rowname))


# add on the prediction
resultdtw = observdataWideWiggle # make a copy of original data 
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
  avg=c(avg,finds(mean,resultdtw[i,6:10],resultdtw$visits[i]))
  sums=c(sums,finds(sum,resultdtw[i,6:10],resultdtw$visits[i]))
}

resultdtw$avg= avg
resultdtw$sum= sums

# have a look at the fit 
par(mfrow=c(2,2))
for (i in c(2,4,8,16)){
  plot(1:20,sitedatawiggle[i,],ylim=c(0,100),ylab=paste("pop size site",i),xlab="Year") # real data - site 8
  points(1:20,filter(resultdtw,site==i)$Nest,col="red") # N mixture
  points(1:20,filter(resultdtw,site==i)$avg,col="blue") # mean
  points(1:20,filter(resultdtw,site==i)$sum,col="green") # sum
}

# put this back
par(mfrow=c(1,1))

# With nimble I keep getting a very low p parameter, which drives up the intercept of the fit beyond that observed with the stan
# model. I'm not sure why this is. I think there is a lot of scope for working through this, as nimble has a sampler that requires 
# a good amount more configuration than the stan HMC 



# 2) Site level spline with shrinkage to global fit and a site level intercept 

# To make this worthwhile I'm going to add variation 
# I'm going to add variation in the wiggle based on x y position  
# heres the idea
plot(sinwig * ((xsamp[1]-median(xsamp))  + (ysamp[1]-median(ysamp)) ))
points(sinwig * ((max(xsamp)-median(xsamp))*0.08  + (max(ysamp)-median(ysamp))*0.1 ))
points(sinwig * 5,col="blue")
points(sinwig * ((min(xsamp)-median(xsamp))*0.1  + (min(ysamp)-median(ysamp))*0.1 ),col="red")

sitedatawigglexy = matrix(NA,nrow = length(xsamp),ncol = years)
for(k in 1:ncol(sitedatawigglexy)) { 
  for(i in 1:nrow(sitedatawigglexy)) {
    sitedatawigglexy[i,k] = amap[amap$x==xsamp[i] & amap$y==ysamp[i],]$sim1 + (sinwig*((xsamp[i]-median(xsamp))*0.08 +(ysamp[i]-median(ysamp))*0.08 ))[k]
  }
}

# new vals
plot(NULL,xlim=c(0,20),ylim=c(0,40))
for(i in 1:20){lines(1:20,sitedatawigglexy[i,],col=i)}
# as compared to
plot(NULL,xlim=c(0,20),ylim=c(0,40))
for(i in 1:20){lines(1:20,sitedatawiggle[i,],col=i)}

## make the data
observdatawigglexy = data.frame()
for(j in 1:ncol(sitedatawigglexy)) { # each year
  for(i in 1:nrow(sitedatawigglexy)) {
    visits = round(runif(1,min=1,max=5))
    counts = round(rpois(visits,sitedatawigglexy[i,j]) * p)
    td=data.frame(site=i,year=j+1990,count=counts,x=xsamp[i],y=ysamp[i])
    observdatawigglexy=rbind(observdatawigglexy,td)
    
  }
}

minpopswig = observdatawigglexy %>% group_by(site,year) %>%
  summarise(minpop = max(count))
observdatawigglexy = left_join(observdatawigglexy,minpopswig,by=c("site","year"))

observdatawigglexy=observdatawigglexy %>% group_by(site,year) %>% mutate(visit=row_number())

ggplot(observdatawigglexy,aes(year,count))+
  geom_point(size=2)+
  facet_wrap(~site)+
  theme_classic() +
  theme(axis.text = element_text(size=14),axis.title = element_text(size=18))+
  ylab("Count")+
  xlab("Year")

# sorry about these names 
observdataWideWigglexy=observdatawigglexy %>% pivot_wider(names_from=visit,values_from= count)
observdataWideWigglexy=replace_na(observdataWideWigglexy,list('1'=-1,'2'=-1,'3'=-1,'4'=-1,'5'=-1)) # pad with -1
observdataWideWigglexy$visits=apply(observdataWideWigglexy[,6:10] >=0,1,sum)

# fit the mode1
# build list
SimListData2 = list(
  count = as.matrix(observdataWideWigglexy[,6:10]),
  splineMat = sp[1:length(unique(observdataWideWigglexy$year)),1:ncol(sp)]
)

# set up data 
SimCon2 = list(
  n = nrow(observdataWideWigglexy),
  nsite = length(unique(observdataWideWigglexy$site)),
  nyear = length(unique(observdataWideWigglexy$year)),
  visits = observdataWideWigglexy$visits,
  site = observdataWideWigglexy$site,
  year = observdataWideWigglexy$year - min(observdataWideWigglexy$year) + 1,
  nspline = 7
)

# initial values for parameters
Inits2 = list(
  N = apply(observdataWideWigglexy[,6:10], 1, max), 
  beta1 = matrix(0.01,nrow = 7,ncol=length(unique(observdataWideWigglexy$site))), # 6 knots so 7 splines
  beta0 = rep(0,length(unique(observdataWideWigglexy$site))),
  delta0 = 0,
  hyperbeta = rep(0,7),
  sigmabeta =rep(0,7)
  
)


nmixsplineshrink <- nimbleCode({
  
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
nmixsplineshrink2 = nimbleModel(nmixsplineshrink, constants = SimCon2,data=SimListData2,inits = Inits2 )

watchlist =c("beta0","beta1","delta0","p","N.pred")

mcmcconf2=configureMCMC(nmixsplineshrink2) 
mcmcconf2$addMonitors(watchlist)
mcmcbuild2 = buildMCMC(mcmcconf2) 
crun2 = compileNimble(nmixsplineshrink2) 
crunMCMC2 = compileNimble(mcmcbuild2, project = nmixsplineshrink2) 
samples2 = runMCMC(crunMCMC2,  niter = 22000,
                  nburnin = 2000,
                  thin = 10,
                  nchains=1,
                  samplesAsCodaMCMC = T)



G2 = summary(samples2)
G2 = G2$statistics
G2 = rownames_to_column(data.frame(G2))
# extract the N estimates 
Nestimates = filter(G2, grepl("N",rowname))

# add on the prediction
resultdtw2 = observdataWideWigglexy # make a copy of original data 
resultdtw2$Nest= Nestimates$Mean

avg=c()
sums=c()
for(i in 1:400) {
  avg=c(avg,finds(mean,observdataWideWigglexy[i,6:10],observdataWideWigglexy$visits[i]))
  sums=c(sums,finds(sum,observdataWideWigglexy[i,6:10],observdataWideWigglexy$visits[i]))
}

resultdtw2$avg= avg
resultdtw2$sum= sums

# have a look at the fit - here the model fits the change across sites with a shrink back to
# the overall mean 
par(mfrow=c(2,2))
for (i in c(2,4,8,16)){
  plot(1:20,sitedatawigglexy[i,],ylim=c(0,200),ylab=paste("pop size site",i),xlab="Year") # real data - site 8
  points(1:20,filter(resultdtw2,site==i)$Nest,col="red") # N mixture
  points(1:20,filter(resultdtw2,site==i)$avg,col="blue") # mean
  points(1:20,filter(resultdtw2,site==i)$sum,col="green") # sum
}


# 3) Site level GAM with no shrinkage and site level intercept

# we can use the same dataset for the non-shrinkage version

# fit the mode1
# build list
SimListData3 = list(
  count = as.matrix(observdataWideWigglexy[,6:10]),
  splineMat = sp[1:length(unique(observdataWideWigglexy$year)),1:ncol(sp)]
)

# set up data 
SimCon3 = list(
  n = nrow(observdataWideWigglexy),
  nsite = length(unique(observdataWideWigglexy$site)),
  nyear = length(unique(observdataWideWigglexy$year)),
  visits = observdataWideWigglexy$visits,
  site = observdataWideWigglexy$site,
  year = observdataWideWigglexy$year - min(observdataWideWigglexy$year) + 1,
  nspline = ncol(sp)
)

# initial values for parameters
Inits3 = list(
  N = apply(observdataWideWigglexy[,6:10], 1, max), 
  beta1 = matrix(0.01,nrow = 7,ncol=length(unique(observdataWideWigglexy$site))), # 6 knots so 7 splines
  beta0 = rep(0,length(unique(observdataWideWigglexy$site))),
  delta0 = 0
  
)


nmixsplinefree <- nimbleCode({
  
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
  

  
  # Ecological model for true abundance
  for(i in 1:n) {
    log(lambda[i]) <- beta0[site[i]] + sum( (splineMat[1:nyear,1:nspline] %*% beta1[1:nspline,site[i]])[1:year[i],1] )

    # Observation model for replicated count
    count[i,1:visits[i]] ~ dNmixture_s(lambda=lambda[i],p=p,Nmin=-1,Nmax=-1,len=visits[i])
    
  }
  
  # Derived quantities
  for(j in 1:nsite){      
    for(k in 1:nyear){
      N.pred[j, k] <- exp(beta0[j] + sum( (splineMat[1:nyear,1:nspline] %*% beta1[1:nspline,j])[1:k,1] ) )
    }
  }
})

# convert to nimble code
nmixsplinefree3 = nimbleModel(nmixsplinefree, constants = SimCon3,data=SimListData3,inits = Inits3 )

watchlist =c("beta0","beta1","delta0","p","N.pred")

mcmcconf3=configureMCMC(nmixsplinefree3) 
mcmcconf3$addMonitors(watchlist)
mcmcbuild3 = buildMCMC(mcmcconf3)
crun3 = compileNimble(nmixsplinefree3) 
crunMCMC3 = compileNimble(mcmcbuild3, project = nmixsplinefree3) 
samples3 = runMCMC(crunMCMC3,  niter = 22000,
                   nburnin = 2000,
                   thin = 10,
                   nchains=1,
                   samplesAsCodaMCMC = T)



G3 = summary(samples3)
G3 = G3$statistics
G3 = rownames_to_column(data.frame(G3))

# extract the N estimates 
Nestimates = filter(G3, grepl("N",rowname))

# add on the prediction
resultdtw3 = observdataWideWigglexy # make a copy of original data 
resultdtw3$Nest= Nestimates$Mean

avg=c()
sums=c()
for(i in 1:400) {
  avg=c(avg,finds(mean,observdataWideWigglexy[i,6:10],observdataWideWigglexy$visits[i]))
  sums=c(sums,finds(sum,observdataWideWigglexy[i,6:10],observdataWideWigglexy$visits[i]))
}

resultdtw3$avg= avg
resultdtw3$sum= sums

# have a look at the fit - this seems much more wobbly than the stan fit - again interesting 
par(mfrow=c(2,2))
for (i in c(2,4,8,16)){
  plot(1:20,sitedatawigglexy[i,],ylim=c(0,100),ylab=paste("pop size site",i),xlab="Year") # real data - site 8
  points(1:20,filter(resultdtw3,site==i)$Nest,col="red") # N mixture
  points(1:20,filter(resultdtw3,site==i)$avg,col="blue") # mean
  points(1:20,filter(resultdtw3,site==i)$sum,col="green") # sum
}



# 4) Global GAM with site level intercept - but now using a vector of detection per visit 


# initial values for parameters
Inits1 = list(
  N = apply(observdataWideWiggle[,6:10], 1, max), 
  beta1 = rep(0.01,7), # 6 knots so 7 splines
  beta0 = rep(0,length(unique(observdataWideWiggle$site))),
  delta0 = 0
)



nmixspline <- nimbleCode({
  
  # Priors
  for(j in 1:nsite) {             
    beta0[j] ~ dnorm(2,sd=1)  #  Intercepts on ecology
  }
  
  for(s in 1:nspline) {
    beta1[s] ~ dnorm(0, sd=0.05)  # slope on ecology
  }
  
  delta0 ~ dnorm(0,sd=1.6) # intercept on observation
  
  yearweights[1:nyear] <- splineMat[1:nyear,1:nspline] %*%  beta1[1:nspline]
  # this is for the observation model (on the logit scale)
  
  # Ecological model for true abundance
  for(i in 1:n) {
    logit(p[i,1:visits[i]]) <- delta0
    log(lambda[i]) <- beta0[site[i]] + sum(yearweights[1:year[i]])

    count[i,1:visits[i]] ~ dNmixture_v(lambda=lambda[i],p=p[i,1:visits[i]],Nmin=-1,Nmax=-1,len=visits[i])
  }
  
  # Derived quantities
  for(j in 1:nsite){        
    for(k in 1:nyear){    
      N.pred[j, k] <- exp(beta0[j] + sum(yearweights[1:k]))
    }
  }
})


# convert to nimble code
nmixspline1 = nimbleModel(nmixspline, constants = SimCon1,data=SimListData1,inits = Inits1 )



# 5) Here I'm going to run the shrink model and 
# plot a version similar to that in the model overview 
# I am going to double abundance so abundance doesn't hover around zero

Site1N = round(sitedatawigglexy[[4,1]] * 2) # starting pop site 1
Site2N = round(sitedatawigglexy[[14,1]] * 2) # starting pop site 2

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
  Site1N = round(sitedatawigglexy[[4,i]] * 2)
  Site2N = round(sitedatawigglexy[[14,i]] * 2) 
  
}



SimListDatap = list(
  count = as.matrix(FixedSizeDataset1[,5:7]),
  splineMat = sp[1:length(unique(FixedSizeDataset1$year)),1:ncol(sp)]
)

# set up data 
SimConp = list(
  n = nrow(FixedSizeDataset1),
  nsite = length(unique(FixedSizeDataset1$site)),
  nyear = length(unique(FixedSizeDataset1$year)),
  visits = rep(3,nrow(FixedSizeDataset1)),
  site = FixedSizeDataset1$site,
  year = FixedSizeDataset1$year - min(FixedSizeDataset1$year) + 1,
  nspline = ncol(sp)
)

# initial values for parameters
Initsp = list(
  N = apply(FixedSizeDataset1[,5:7], 1, max), 
  beta1 = matrix(0.01,nrow = 7,ncol=length(unique(FixedSizeDataset1$site))), # 6 knots so 7 splines
  beta0 = rep(0.01,length(unique(FixedSizeDataset1$site))),
  delta0 = 0.01
  
)

nmixsplinefreep = nimbleModel(nmixsplinefree, constants = SimConp,data=SimListDatap,inits = Initsp)

watchlist =c("beta0","beta1","delta0","p","N.pred")

mcmcconfp=configureMCMC(nmixsplinefreep) 
mcmcconfp$addMonitors(watchlist)
mcmcbuildp = buildMCMC(mcmcconfp) # build it 
crunp = compileNimble(nmixsplinefreep) # you must compile the original model
crunMCMCp = compileNimble(mcmcbuildp, project = nmixsplinefreep) 
samplesp = runMCMC(crunMCMCp,  niter = 22000,
                   nburnin = 2000,
                   thin = 10,
                   nchains=1,
                   samplesAsCodaMCMC = T)



Gp = summary(samplesp)
Gp = Gp$statistics
Gp = rownames_to_column(data.frame(Gp))

# extract the N estimates 
Nestimates = filter(Gp, grepl("N",rowname))
FixedSizeDataset1$Nest= Nestimates$Mean
FixedSizeDataset1$avg = apply(FixedSizeDataset1[,5:7],1,mean)

par(mfrow=c(1,2))
for(i in c(1,2)){
  plot(1:20,filter(FixedSizeDataset1,site==i)$wholesquare,ylim=c(0,70),ylab=paste("Population size site",i),xlab="Year",pch=19,col="grey",cex=1.2,cex.lab=1.2) 
  points(1:20,filter(FixedSizeDataset1,site==i)$Nest,col="red",pch=19,cex=1.2) # N mixture
  points(1:20,filter(FixedSizeDataset1,site==i)$avg,col="black",pch=19,cex=1.2) # mean
}
legend(1, 67, legend=c("True value", "Mean of counts","Nmixture estimate approx"),
       fill =c("grey", "black","red"), lty=1:2, cex=0.8)



