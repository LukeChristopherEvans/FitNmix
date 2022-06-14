library(tidyverse)
library(splines)
library(rstan)
library(rethinking)
library(coda)

folderpath = "/home/lukee/Insync/vs917256@reading.ac.uk/OneDrive Biz/FitNmix/Data/RSPB/"
modelpath = "/home/lukee/Insync/vs917256@reading.ac.uk/OneDrive Biz/FitNmix/Stan files/Bird/"
filepath1 = "BlueTitWide.txt"

# The plan here is to fit the GAM n mixture models. 
# For model selection, I will first try GlobalGamNmix, SiteSpecificShrinkGamNmix, SiteSpecificFreeGamNmix 
# then for the best model I will try one knot lower and higher than the Fewster recommendation, which is N years * 0.3
# model selection will be based on WAIC 

# !Note! we only 'profile' the models here to save computation time, and we, therefore, only avoid divergent transitions 
# as fatal errors for interpretation. For a proper analysis we would need more samples from the posterior
# and to spend more time evaluating chain mixing


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

#new way of deciding k 
brdt$K = brdt$minpop +  qpois(0.99999,brdt$minpop)+3


# build list
Globallist = list(
  n = nrow(brdt),
  nsite = length(unique(brdt$Gridref)),
  nyear = length(unique(brdt$Year)),
  nsplines = ncol(sp),
  
  k = brdt$K,
  minpop = brdt$minpop,
  count = as.matrix(brdt[,9:16]),
  
  visitmatrix = as.matrix(brdt[,18:25]),
  bandmatrix = as.matrix(brdt[,18:25]),
  
  splineMat = sp,
  
  visits = brdt$visitcount,
  site = as.integer(as.factor(brdt$Gridref)),
  year =  brdt$Year - min(brdt$Year) + 1,

  maxvisits =8
  )

# get the file path to the stan model 
filepath = "GlobalGamNmixMultiProbBird.stan"

# this model will tell you off a lot because p is a ragged array -  key thing 
# to watch out for at this stage is divergent transitions. 
globalgam<-rstan::stan(
  file = paste0(modelpath,filepath),
  data=Globallist,
  iter = 1250, # only 250 samples for testing 
  warmup = 1000,
  chains=2, # runs long
  cores=2,
  init_r = 0.01
)

G1=precis(globalgam,depth=3) # extract posterior estimates in nice table 
# make into a dataframe with rownames
G1=rownames_to_column(data.frame(G1))
# extract the N estimates 
Nestimates = filter(G1, grepl("N",rowname))


# add on the prediction
resultdtw = brdt
resultdtw$Nestglobal= Nestimates$mean

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
  facet_wrap(~Gridref)+
  theme_classic()

# Now I can look at the observational model 
plot(globalgam,pars=c("delta0","delta2","delta1"))

##
# 2) Site level GAM with shrinkage to global fit and a site level intercept - rhats better here
localshrinkpath = "SiteSpecificShrinkGamNmixMultiProbBird.stan"
shrinkgam = rstan::stan(
  file = paste0(modelpath,localshrinkpath),
  data=Globallist,
  iter = 1250,
  warmup = 1000,
  chains=2,
  cores=2,
  init_r = 0.01
)


G2=precis(shrinkgam ,depth=3) 
G2=rownames_to_column(data.frame(G2))
Nestimates = filter(G2, grepl("N",rowname))
resultdtw$Nestshrink= Nestimates$mean


filter(resultdtw, Gridref %in% sample(unique(resultdtw$Gridref),20)) %>%
  ggplot(aes(Year,avg))+
  geom_point()+
  geom_point(aes(Year,globalmean),color="red")+
  geom_point(aes(Year,Nestshrink),color="blue")+
  facet_wrap(~Gridref)+
  theme_classic()

plot(shrinkgam,pars=c("delta0","delta2","delta1"))


# 3) Site level GAM with no shrinkage and site level intercept - high rhats again
localfreepath = "SiteSpecificFreeGamNmixMultiProbBird.stan"
freegam = rstan::stan(
  file = paste0(modelpath,localfreepath),
  data=Globallist,
  iter = 1250, 
  warmup = 1000,
  chains=2,
  cores=2,
  init_r = 0.01
)

G3=precis(freegam ,depth=3)  
G3=rownames_to_column(data.frame(G3))
Nestimates = filter(G3, grepl("N",rowname))
resultdtw$Nestfree= Nestimates$mean

filter(resultdtw, Gridref %in% sample(unique(resultdtw$Gridref),20)) %>%
  ggplot(aes(Year,avg))+
  geom_point()+
  geom_point(aes(Year,globalmean),color="red")+
  geom_point(aes(Year,Nestshrink),color="blue")+
  geom_point(aes(Year,Nestfree),color="green")+
  facet_wrap(~Gridref)+
  theme_classic()

# which has lowest WAIC - approximation of cross validation score 
compare(globalgam,shrinkgam,freegam) # however, note validity

#4) Try other knot numbers - using the shrinkage gam 
## try other knots
numknots = 6
knotsl=quantile(1:years,probs = seq(0,1,length.out=numknots))
knotsl[-c(1,numknots)] # gets rid of first and last vals
sp6 = bs(1:years,knots=knotsl[-c(1,numknots)],degree=3,intercept = F) # cubic
plot(NULL,xlim=c(1,years),ylim=c(0,1),ylab="basis",xlab="year")
for(i in 1:ncol(sp6)){lines(1:years,sp6[,i])}

# one lower
birdlistlow= Globallist
birdlistlow$splineMat = sp6
birdlistlow$nsplines = ncol(sp6) 


localshrinkpath = "SiteSpecificShrinkGamNmixBird.stan"
shrinkgamlow = rstan::stan(
  file = paste0(modelpath,localshrinkpath),
  data=birdlistlow,
  iter = 1250, 
  warmup = 1000,
  chains=4,
  cores=4,
  init_r = 0.01
)


G2low=precis(shrinkgamlow ,depth=3) 
G2low=rownames_to_column(data.frame(G2low))
Nestimates = filter(G2low, grepl("N",rowname))
resultdtw$Nestshrinklow= Nestimates$mean

filter(resultdtw, Gridref %in% sample(unique(resultdtw$Gridref),20)) %>%
  ggplot(aes(Year,avg))+
  geom_point()+
  geom_point(aes(Year,Nestshrinklow),color="lightblue")+
  geom_point(aes(Year,Nestshrink),color="blue")+
  facet_wrap(~Gridref)+
  theme_classic()


# one higher
numknots = 8
knotsl=quantile(1:years,probs = seq(0,1,length.out=numknots))
knotsl[-c(1,numknots)] # gets rid of first and last vals
sp8 = bs(1:years,knots=knotsl[-c(1,numknots)],degree=3,intercept = F) # cubic
plot(NULL,xlim=c(1,years),ylim=c(0,1),ylab="basis",xlab="year")
for(i in 1:ncol(sp8)){lines(1:years,sp8[,i])}

birdlisthigh= Globallist
birdlisthigh$splineMat = sp8
birdlisthigh$nsplines = ncol(sp8) 


shrinkgamhigh = rstan::stan(
  file = paste0(modelpath,localshrinkpath),
  data=birdlisthigh,
  iter = 1250, 
  warmup = 1000,
  chains=4,
  cores=4,
  init_r = 0.01
)


G2high=precis(shrinkgamhigh ,depth=3) 
G2high=rownames_to_column(data.frame(G2high))
Nestimates = filter(G2high, grepl("N",rowname))
resultdtw$Nestshrinkhigh= Nestimates$mean

filter(resultdtw, Gridref %in% sample(unique(resultdtw$Gridref),20)) %>%
  ggplot(aes(Year,avg))+
  geom_point()+
  geom_point(aes(Year,Nestshrinklow),color="lightblue")+
  geom_point(aes(Year,Nestshrink),color="blue")+
  geom_point(aes(Year,Nestshrinkhigh),color="darkblue")+
  facet_wrap(~Gridref)+
  theme_classic()

compare(shrinkgamlow,shrinkgam,shrinkgamhigh)
