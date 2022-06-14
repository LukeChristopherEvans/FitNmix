library(gstat)
library(tidyverse)
library(rethinking)
library(rstan)

# Here we develop and explore GAM based fits as these are the 
# most realistic approach to modeling changes in population abundance 
# over time.
# 1) Global GAM with site level intercept 
# 2) Site level GAM with shrinkage to global fit and a site level intercept 
# 3) Site level GAM with no shrinkage and site level intercept
# 4) Site level GAM with geographic distance based shrinkage 


# !Note! we only 'profile' the models here to save computation time, and we, therefore, only avoid divergent transitions 
# as fatal errors for interpretation. For a proper analysis we would need more samples from the posterior
# and to spend more time evaluating chain mixing

folderpath = "/home/lukee/Insync/vs917256@reading.ac.uk/OneDrive Biz/FitNmix"

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

sp = bs(1:years,knots=knotsl[-c(1,numknots)],degree=3,intercept = F) # cubic spline
plot(NULL,xlim=c(1,20),ylim=c(0,1),ylab="basis",xlab="year") #basis weights across time 
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



# get the file path to the stan model 
filepath = "/Stan files/GlobalGamNmix.stan"

# 1) Global spline with site level intercept 
GGam1list = list(
  nsplines = ncol(sp),
  splineMat = sp[1:length(unique(observdataWideWiggle$site)),1:ncol(sp)],
  site = observdataWideWiggle$site,
  year = observdataWideWiggle$year - min(observdataWideWiggle$year) + 1,
  count = as.matrix(observdataWideWiggle[,6:10]),
  n = nrow(observdataWideWiggle),
  maxvisits = 5,
  nsite = length(unique(observdataWideWiggle$site)),
  nyear = length(unique(observdataWideWiggle$year)),
  visits = observdataWideWiggle$visits,
  k = 40, 
  minpop =observdataWideWiggle$minpop 
)

nmixspline<-rstan::stan(
  file = paste0(folderpath,filepath),
  data=GGam1list,
  iter = 1250, # only 250 samples for testing 
  warmup = 1000,
  chains=4,
  cores=4,
  init_r = 0.01
)

G1=precis(nmixspline,depth=3) # extract posterior estimates in nice table 
# make into a dataframe with rownames
G1=rownames_to_column(data.frame(G1))
# extract the N estimates 
Nestimates = filter(G1, grepl("N",rowname))

# add on the prediction
resultdtw = observdataWideWiggle # make a copy of original data 
resultdtw$Nest= Nestimates$mean

# quick function for grabbing data from the ragged array and applying a function
finds=function(f,x,y){
  z=c()
  for(i in 1:y) {
    z =c(z, pull(x[1,i]))
  }
  f(z)
}

avg=c()
sums=c()
for(i in 1:400) {
  avg=c(avg,finds(mean,observdataWideWiggle[i,6:10],observdataWideWiggle$visits[i]))
  sums=c(sums,finds(sum,observdataWideWiggle[i,6:10],observdataWideWiggle$visits[i]))
}

resultdtw$avg= avg
resultdtw$sum= sums

# have a look at the fit - remembering caveats from ModelOverviewSimulationStan this is doing pretty well 
par(mfrow=c(2,2))
for (i in c(2,4,8,16)){
  plot(1:20,sitedatawiggle[i,],ylim=c(0,40),ylab=paste("pop size site",i),xlab="Year") # real data - site 8
  points(1:20,filter(resultdtw,site==i)$Nest,col="red") # N mixture
  points(1:20,filter(resultdtw,site==i)$avg,col="blue") # mean
  points(1:20,filter(resultdtw,site==i)$sum,col="green") # sum
}

par(mfrow=c(1,1))

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

GGam2list = list(
  nsplines = ncol(sp),
  splineMat = sp[1:length(unique(observdataWideWiggle$site)),1:ncol(sp)],
  site = observdataWideWigglexy$site,
  year = observdataWideWigglexy$year - min(observdataWideWigglexy$year) + 1,
  count = as.matrix(observdataWideWigglexy[,6:10]),
  n = nrow(observdataWideWigglexy),
  maxvisits = 5,
  nsite = length(unique(observdataWideWigglexy$site)),
  nyear = length(unique(observdataWideWigglexy$year)),
  visits = observdataWideWigglexy$visits,
  k = 40, 
  minpop =observdataWideWigglexy$minpop 
)

filepathshrink = "/Stan files/SiteSpecificShrinkGamNmix.stan"

# got 1 divergent transition - so should probably up the adapt delta
nmixshrinkspline<-rstan::stan(
  file = paste0(folderpath,filepathshrink),
  data=GGam2list ,
  iter = 1250, # only 250 samples for testing 
  warmup = 1000,
  chains=4,
  cores=4,
  init_r = 0.01
)

G2=precis(nmixshrinkspline,depth=3) # extract posterior estimates in nice table 
# make into a dataframe with rownames
G2=rownames_to_column(data.frame(G2))
# extract the N estimates 
Nestimates = filter(G2, grepl("N",rowname))

# add on the prediction
resultdtw2 = observdataWideWigglexy # make a copy of original data 
resultdtw2$Nest= Nestimates$mean

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
  plot(1:20,sitedatawigglexy[i,],ylim=c(0,40),ylab=paste("pop size site",i),xlab="Year") # real data - site 8
  points(1:20,filter(resultdtw2,site==i)$Nest,col="red") # N mixture
  points(1:20,filter(resultdtw2,site==i)$avg,col="blue") # mean
  points(1:20,filter(resultdtw2,site==i)$sum,col="green") # sum
}



# 3) Site level GAM with no shrinkage and site level intercept

# we can use the same dataset for the non-shrinkage version

filepathfree = "/Stan files/SiteSpecificFreeGamNmix.stan"

# got 1 divergent transition - so should probably up the adapt delta
nmixfreespline<-rstan::stan(
  file = paste0(folderpath,filepathfree),
  data=GGam2list ,
  iter = 1250, # only 250 samples for testing 
  warmup = 1000,
  chains=4,
  cores=4,
  init_r = 0.01
)

G3=precis(nmixfreespline,depth=3) # extract posterior estimates in nice table 
# make into a dataframe with rownames
G3=rownames_to_column(data.frame(G3))
# extract the N estimates 
Nestimates = filter(G3, grepl("N",rowname))

# add on the prediction
resultdtw3 = observdataWideWigglexy # make a copy of original data 
resultdtw3$Nest= Nestimates$mean

avg=c()
sums=c()
for(i in 1:400) {
  avg=c(avg,finds(mean,observdataWideWigglexy[i,6:10],observdataWideWigglexy$visits[i]))
  sums=c(sums,finds(sum,observdataWideWigglexy[i,6:10],observdataWideWigglexy$visits[i]))
}

resultdtw3$avg= avg
resultdtw3$sum= sums

# have a look at the fit 
par(mfrow=c(2,2))
for (i in c(2,4,8,16)){
  plot(1:20,sitedatawigglexy[i,],ylim=c(0,40),ylab=paste("pop size site",i),xlab="Year") # real data - site 8
  points(1:20,filter(resultdtw3,site==i)$Nest,col="red") # N mixture
  points(1:20,filter(resultdtw3,site==i)$avg,col="blue") # mean
  points(1:20,filter(resultdtw3,site==i)$sum,col="green") # sum
}


# 4) Site level GAM with geographic distance based shrinkage 
# This final model is a bit of fun - it draws the parameters for the splines from multivariate normal (n = sites)
# where the covariance matrix for the multivariate normal is the distance between the sites. That means that the
# slopes from closer sites should be more similar. The huge downside is estimating a 20X20 covariance matrix i.e. 400 extra parameters  

# add a distance matrix
DMat=as.matrix(dist(distinct(ungroup(observdataWideWigglexy),site,.keep_all = T)[,3:4],diag = T,upper=T))
GGam2list$DMat = DMat

filepathdist = "/Stan files/SiteSpecificDistanceGamNmix.stan"

# got 1 divergent transition - so should probably up the adapt delta
nmixdistspline<-rstan::stan(
  file = paste0(folderpath,filepathdist),
  data=GGam2list ,
  iter = 1250, # only 250 samples for testing 
  warmup = 1000,
  chains=4,
  cores=4,
  init_r = 0.01
)

G4=precis(nmixdistspline,depth=3) # extract posterior estimates in nice table 
# make into a dataframe with rownames
G4=rownames_to_column(data.frame(G4))
# extract the N estimates 
Nestimates = filter(G4, grepl("N",rowname))

# add on the prediction
resultdtw4 = observdataWideWigglexy # make a copy of original data 
resultdtw4$Nest= Nestimates$mean

avg=c()
sums=c()
for(i in 1:400) {
  avg=c(avg,finds(mean,observdataWideWigglexy[i,6:10],observdataWideWigglexy$visits[i]))
  sums=c(sums,finds(sum,observdataWideWigglexy[i,6:10],observdataWideWigglexy$visits[i]))
}

resultdtw4$avg= avg
resultdtw4$sum= sums

# have a look at the fit - they shrink to the change across space 
par(mfrow=c(3,3))
for (i in c(2,4,6,8,10,12,14,16,18)){
  plot(1:20,sitedatawigglexy[i,],ylim=c(0,40),ylab=paste("pop size site",i),xlab="Year") # real data - site 8
  points(1:20,filter(resultdtw4,site==i)$Nest,col="red") # N mixture
  points(1:20,filter(resultdtw4,site==i)$avg,col="blue") # mean
  points(1:20,filter(resultdtw4,site==i)$sum,col="green") # sum
}









#### UNDER DEVELOPMENT - NOT WORKING #################

##### NEW UPDATE THE MEEHAN APPROX IN STAN ####
GGam1list$k = 40 

# get the file path to the stan model 
filepath = "/Stan files/meehan4.stan"

nmixmeehan<-rstan::stan(
  file = paste0(folderpath,filepath),
  data=GGam1list,
  iter = 1250, # only 250 samples for testing 
  warmup = 1000,
  chains=4,
  cores=4,
  init_r = 0.01
)

G1=precis(nmixmeehan,depth=3) # extract posterior estimates in nice table 
# make into a dataframe with rownames
G1=rownames_to_column(data.frame(G1))
# extract the N estimates 
Nestimates = filter(G1, grepl("N",rowname))

# add on the prediction
resultdtw = observdataWideWiggle # make a copy of original data 
resultdtw$Nest= Nestimates$mean

# quick function for grabbing data from the ragged array and applying a function
finds=function(f,x,y){
  z=c()
  for(i in 1:y) {
    z =c(z, pull(x[1,i]))
  }
  f(z)
}

avg=c()
sums=c()
for(i in 1:400) {
  avg=c(avg,finds(mean,observdataWideWiggle[i,6:10],observdataWideWiggle$visits[i]))
  sums=c(sums,finds(sum,observdataWideWiggle[i,6:10],observdataWideWiggle$visits[i]))
}

resultdtw$avg= avg
resultdtw$sum= sums

# have a look at the fit - remembering caveats from ModelOverviewSimulationStan this is doing pretty well 
par(mfrow=c(2,2))
for (i in c(2,4,8,16)){
  plot(1:20,sitedatawiggle[i,],ylim=c(0,40),ylab=paste("pop size site",i),xlab="Year") # real data - site 8
  points(1:20,filter(resultdtw,site==i)$Nest,col="red") # N mixture
  points(1:20,filter(resultdtw,site==i)$avg,col="blue") # mean
  points(1:20,filter(resultdtw,site==i)$sum,col="green") # sum
}

