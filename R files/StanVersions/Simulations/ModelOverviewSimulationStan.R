### Model fit simulation ###
# The goal here is to simulate very simple count data from a process similar to BTC and RSPB 
# counts and then fit and explore basic N mixture models using stan 

# !Note! we only 'profile' the models here to save computation time, and we, therefore, only avoid divergent transitions 
# as fatal errors for interpretation. For a proper analysis we would need more samples from the posterior
# and to spend more time evaluating chain mixing

library(gstat)
library(tidyverse)
library(rethinking)
library(rstan)

# Model choices:
# global linear trend with site level intercept: BaseNmixStan, BaseNmixStanDensityversion 

## Stage 1: set up naive simulation and run the base model

folderpath = "~/FitNmix"

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

# now sample from this multiple times using poission error 
# and given detection probability 

p = 0.8 # detection prob
#visits = round(runif(1,min=1,max=5))

observdata = data.frame()
for(j in 1:ncol(sitedata)) { # each year
  for(i in 1:nrow(sitedata)) {
    visits = round(runif(1,min=1,max=5)) # between 1 and 5 visits to the site 
    counts = round(rpois(visits,sitedata[i,j]) * p) # counts 
    td=data.frame(site=i,year=j+1990,count=counts,x=xsamp[i],y=ysamp[i])
    observdata=rbind(observdata,td)
    
  }
}

# going to get max count per season as the minimum of the theoretical
# population size that year - this goes into the stan model 
minpops = observdata %>% group_by(site,year) %>%
  summarise(minpop = max(count))

observdata = left_join(observdata,minpops,by=c("site","year"))

# counts across sites - trend somewhat obscured by sampling variation (also apparent population growth and declines coming solely from sampling error!!)
ggplot(observdata,aes(year,count))+
  geom_point(size=2)+
  facet_wrap(~site)+
  theme_classic()+
  theme(axis.text = element_text(size=14),axis.title = element_text(size=18))+
   ylab("Count")+
   xlab("Year")

# wide format for sending to stan 
observdata=observdata %>% group_by(site,year) %>% mutate(visit=row_number())
observdataw=observdata %>% pivot_wider(names_from=visit,values_from= count)
observdataw=replace_na(observdataw,list('1'=-1,'2'=-1,'3'=-1,'4'=-1,'5'=-1)) # pad with -1
# now add visits
observdataw$visits=apply(observdataw[,6:10] >=0,1,sum)

# set up a list for stan 
simlist1 = list(
  site = observdataw$site,
  year = observdataw$year - min(observdataw$year) + 1,
  count = as.matrix(observdataw[,6:10]),
  n = nrow(observdataw),
  nsite = length(unique(observdataw$site)),
  nyear = length(unique(observdataw$year)),
  visits = observdataw$visits,
  maxvisits = 5,
  # stuff to send in to help
  k = 40, # upper feasible pop size - minpop cannot be above this otherwise stan bugs (we need this to integrate over j)
  minpop =observdataw$minpop # minimum possible size
)



# get the file path to the stan model 
filepath = "/Stan files/BaseNmixStan.stan"

# model 1 BaseNmixStan is a random intercept shared slope model
# on the ecological parameters and 
# run the model
set.seed(1234)
nmix1=rstan::stan(
  file = paste0(folderpath,filepath),
  data=simlist1,
  iter = 1250, # only 250 samples for testing 
  warmup = 1000,
  chains=4,
  cores=4,
  init_r = 0.01 
)


traceplot(nmix1) # check chains - also check rhats
WAIC(nmix1) # simulated cross validation score
ab=precis(nmix1,depth=3) # extract posterior estimates in nice table 

# make into a dataframe with rownames
ab=rownames_to_column(data.frame(ab))
# extract the N estimates 
Nestimates = filter(ab, grepl("N",rowname))

# add on the prediction
resultdtw = observdataw # make a copy of original data 
resultdtw$Nest= Nestimates$mean

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
plot(1:20,sitedata[i,],ylim=c(0,35),ylab=paste("pop size site",i),xlab="Year") # real data - site 8
points(1:20,filter(resultdtw,site==i)$Nest,col="red") # N mixture
points(1:20,filter(resultdtw,site==i)$avg,col="blue") # mean
points(1:20,filter(resultdtw,site==i)$sum,col="green") # sum
}

# put this back
par(mfrow=c(1,1))



## Stage 2: Features of the stan model  

# Why does it overestimate abundance?

# We see above that the model is very effective at estimating 
# the trend but overestimates the abundance. This is because the nmixture model is
# set up in a way that we imagine we have a fixed population (not a density - an actual number of individuals) per site and per year 
# that we sample from and which cannot be smaller than the maximum observed value (in the range minpop:k).
# Effectively the lambda in the simulation is really more of a latent density and not a latent population size - causing the mismatch. Or another way of saying
# it, is that the lambda parameter shows what is the likely mean number I would observe in a sample given a fixed population in the local region.

# We can see this feature of the model here in an example with an actual fixed population size: 
# We imagine we have set sized populations, growing by one individual per year, observed while moving around randomly in a 1x1 grid (i.e. a poisson process).
# we observe half the square and so have the possibility of seeing half of these 
# animals on average (lambda is 0.5 * N). I also set the observation prob centered on 80%
# we'll complete the simulation at just two sites 

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


# add the min population size for marginalizing N over- this is 
# the key unusual aspect of the model 
minpop1 = FixedSizeDataset1%>% group_by(site,year) %>%
  summarise(minpop = max(x1,x2,x3))
minpop2 = FixedSizeDataset2%>% group_by(site,year) %>%
  summarise(minpop = max(x1,x2,x3))

# join datasets
FixedSizeDataset1 = left_join(FixedSizeDataset1,minpop1,by=c("site","year"))
FixedSizeDataset2 = left_join(FixedSizeDataset2,minpop2,by=c("site","year"))

# make lists 
simlist2.1 = list(
  site = FixedSizeDataset1$site,
  year = FixedSizeDataset1$year - min(FixedSizeDataset1$year) + 1,
  count = as.matrix(FixedSizeDataset1[,5:7]),
  n = nrow(FixedSizeDataset1),
  nsite = length(unique(FixedSizeDataset1$site)),
  nyear = length(unique(FixedSizeDataset1$year)),
  visits = rep(3,nrow(FixedSizeDataset1)),
  maxvisits =3,
  # stuff to send in to help
  k = 50, # upper feasible pop size - minpop cannot be above this otherwise stan bugs (we need this to integrate over j)
  minpop =FixedSizeDataset1$minpop # minimum possible size
)
simlist2.2 = list(
  site = FixedSizeDataset2$site,
  year = FixedSizeDataset2$year - min(FixedSizeDataset2$year) + 1,
  count = as.matrix(FixedSizeDataset2[,5:7]),
  n = nrow(FixedSizeDataset2),
  nsite = length(unique(FixedSizeDataset2$site)),
  nyear = length(unique(FixedSizeDataset2$year)),
  visits = rep(3,nrow(FixedSizeDataset2)),
  maxvisits =3,
  # stuff to send in to help
  k = 50, # upper feasible pop size - minpop cannot be above this otherwise stan bugs (we need this to integrate over j)
  minpop =FixedSizeDataset2$minpop # minimum possible size
)

# fit models 
set.seed(1234)
nmix2.1=rstan::stan(
  file = paste0(folderpath,filepath),
  data=simlist2.1,
  iter = 1250, 
  warmup = 1000,
  chains=4,
  cores=4,
  init_r = 0.01 
)
set.seed(1234)
nmix2.2=rstan::stan(
  file = paste0(folderpath,filepath),
  data=simlist2.2,
  iter = 1250, 
  warmup = 1000,
  chains=4,
  cores=4,
  init_r = 0.01 
)

# extract data and add data 
ab2.1=precis(nmix2.1,depth=3)
ab2.1=rownames_to_column(data.frame(ab2.1))
Nestimates2.1 = filter(ab2.1, grepl("N",rowname))
FixedSizeDataset1$Nest= Nestimates2.1$mean

ab2.2=precis(nmix2.2,depth=3)  
ab2.2=rownames_to_column(data.frame(ab2.2))
Nestimates2.2 = filter(ab2.2, grepl("N",rowname))
FixedSizeDataset2$Nest= Nestimates2.2$mean

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
stan_dens(nmix2.1,pars=c("p"))

# here is the estimates when sampling the whole square (exactly the same N estimates)
plot(1:20,filter(FixedSizeDataset2,site==1)$wholesquare,ylim=c(0,50),ylab=paste("pop size site",1),xlab="Year") # data from the whole square
points(1:20,filter(FixedSizeDataset2,site==1)$Nest,col="red") # N mixture
points(1:20,filter(FixedSizeDataset2,site==1)$avg,col="blue") # mean

plot(1:20,filter(FixedSizeDataset2,site==2)$wholesquare,ylim=c(0,50),ylab=paste("pop size site",1),xlab="Year") # real data - site 8
points(1:20,filter(FixedSizeDataset2,site==2)$Nest,col="red") # N mixture
points(1:20,filter(FixedSizeDataset2,site==2)$avg,col="blue") # mean

# but now better estimate of the 0.8 observation error 
stan_dens(nmix2.2,pars=c("p"))


# What to do, if we want density?

# An approach is to model the population as a density (not the finite pop size)
# and to have the marginalizing bounded between the observed count per visit and k respectively  
# I implement this in the stan model below 

# we don't need minpop because we only use each count (e.g. count[1:3]) as the baseline - the big penalty though is computation speed 
# as marginalise across every count*site*year rather than every site*year
filepathdens = "/Stan files/BaseNmixStanDensityVersion.stan" # this is called the density version
nmix1dens=rstan::stan(
  file = paste0(folderpath,filepathdens),
  data=simlist2.1,
  iter = 1250, 
  warmup = 1000,
  chains=4,
  cores=4,
  init_r = 0.01
)

ab3=precis(nmix1dens,depth=3)
ab3=rownames_to_column(data.frame(ab3))
Nestimates3 = filter(ab3, grepl("N",rowname))
FixedSizeDataset1$Nestdens= Nestimates3$mean

# now this estimates the lambda as a density correctly 
plot(1:20,filter(FixedSizeDataset1,site==1)$wholesquare,ylim=c(0,50),ylab=paste("pop size site",1),xlab="Year") # data from the whole square
points(1:20,filter(FixedSizeDataset1,site==1)$halfsquare,col="grey") # half square (Expected lambda)
points(1:20,filter(FixedSizeDataset1,site==1)$Nest,col="red") # original N mixture
points(1:20,filter(FixedSizeDataset1,site==1)$avg,col="blue") # mean
points(1:20,filter(FixedSizeDataset1,site==1)$Nestdens,col="green") # dens n mixture

plot(1:20,filter(FixedSizeDataset1,site==2)$wholesquare,ylim=c(0,50),ylab=paste("pop size site",1),xlab="Year") # data from the whole square
points(1:20,filter(FixedSizeDataset1,site==2)$halfsquare,col="grey") # half square (Expected lambda)
points(1:20,filter(FixedSizeDataset1,site==2)$Nest,col="red") # original N mixture
points(1:20,filter(FixedSizeDataset1,site==2)$avg,col="blue") # mean
points(1:20,filter(FixedSizeDataset1,site==2)$Nestdens,col="green") # dens n mixture


# nice output for report 
par(mfrow=c(1,2))
for(i in c(1,2)){
  plot(1:20,filter(FixedSizeDataset1,site==i)$wholesquare,ylim=c(0,70),ylab=paste("Population size site",i),xlab="Year",pch=19,col="grey",cex=1.2,cex.lab=1.2) 
  points(1:20,filter(FixedSizeDataset1,site==i)$Nest,col="red",pch=19,cex=1.2) # N mixture
  points(1:20,filter(FixedSizeDataset1,site==i)$avg,col="black",pch=19,cex=1.2) # mean
}
legend(1, 67, legend=c("True value", "Mean of counts","Nmixture estimate"),
       fill =c("grey", "black","red"), lty=1:2, cex=0.8)

# lots of uncertainty here but peaks over correct value
stan_dens(nmix1dens,pars=c("p"))


# look at computation time comparison 
dentime = system.time(rstan::stan(
  file = paste0(folderpath,filepathdens),
  data=simlist2.1,
  iter = 1250, 
  warmup = 1000,
  chains=4,
  cores=4,
  init_r = 0.01
))

nondentime = system.time(rstan::stan(
  file = paste0(folderpath,filepath),
  data=simlist2.1,
  iter = 1250, 
  warmup = 1000,
  chains=4,
  cores=4,
  init_r = 0.01 
))

# not too bad with this data size but it may scale badly
dentime ; nondentime


## 1 chain version for comparison with nimble speed (note we have to be careful with comparisons because 
# stan samples much more efficiently (i.e higher effective sample size per iteration))
dentime = system.time(rstan::stan(
  file = paste0(folderpath,filepathdens),
  data=simlist2.1,
  iter = 1250, 
  warmup = 1000,
  chains=1,
  cores=1,
  init_r = 0.01
))
nondentime = system.time(rstan::stan(
  file = paste0(folderpath,filepath),
  data=simlist2.1,
  iter = 1250, 
  warmup = 1000,
  chains=1,
  cores=1,
  init_r = 0.01 
))
dentime ; nondentime
