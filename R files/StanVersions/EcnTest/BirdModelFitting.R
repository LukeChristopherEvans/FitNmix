library(tidyverse)
library(splines)
library(rstan)
library(rethinking)

folderpath = "/home/lukee/Insync/vs917256@reading.ac.uk/OneDrive Biz/FitNmix/Data/ECN/"
modelpath = "/home/lukee/Insync/vs917256@reading.ac.uk/OneDrive Biz/FitNmix/Stan files/"
filepath1 = "ProcessedBlueTitEcn.txt"

# The plan here is to fit the GAM n mixture models. 
# For model selection, I will first try GlobalGamNmix, SiteSpecificShrinkGamNmix, SiteSpecificFreeGamNmix 
# then for the best model I will try one knot lower and higher than the Fewster recommendation, which is N years * 0.3
# model selection will be based on WAIC 

# !Note! we only 'profile' the models here to save computation time, and we, therefore, only avoid divergent transitions 
# as fatal errors for interpretation. For a proper analysis we would need more samples from the posterior
# and to spend more time evaluating chain mixing


# load data 
brdt = read.table(paste0(folderpath,filepath1))

# make spline matrix
years =length(unique(brdt2$year))
numknots = round(years*0.3) # years * 0.3 Fewster recommendation 
knotsl=quantile(1:years,probs = seq(0,1,length.out=numknots))
knotsl[-c(1,numknots)] # gets rid of first and last vals
sp = bs(1:years,knots=knotsl[-c(1,numknots)],degree=3,intercept = F) # cubic
plot(NULL,xlim=c(1,years),ylim=c(0,1),ylab="basis",xlab="year")
for(i in 1:ncol(sp)){lines(1:years,sp[,i])}


birdlist1 = list(
  n = nrow(brdt),
  nsite = length(unique(brdt$SITECODE)),
  nyear = years,
  nsplines = ncol(sp),
  maxvisits = 6,
  visits = brdt$visits,
  k = round(quantile(brdt$minpop,probs=0.99999))[[1]], # must be above max observed otherwise bugs (note this cuts of p parameter)
  site = as.integer(as.factor(brdt$SITECODE)),
  year = brdt$year - min(brdt$year) + 1,
  count = as.matrix(brdt[,11:16]),
  minpop = brdt$minpop,
  splineMat = sp
)

globalpath = "GlobalGamNmix.stan"
globalgam = rstan::stan(
  file = paste0(modelpath,globalpath),
  data=birdlist1,
  iter = 1250, # only 250 samples for testing
  warmup = 1000,
  chains=4,
  cores=4,
  init_r = 0.01
)

G1=precis(globalgam,depth=3) # extract posterior estimates in nice table 
# make into a dataframe with rownames
G1=rownames_to_column(data.frame(G1))
# extract the N estimates 
Nestimates = filter(G1, grepl("N",rowname))

# add on the prediction
resultdtw = brdt # make a copy of original data 
resultdtw$Nest= Nestimates$mean

ggplot(resultdtw,aes(year,meanind))+
  geom_point()+
  geom_point(aes(year,Nest),color="red")+
  geom_line()+
  facet_wrap(~SITECODE)+
  theme_classic()


##
# 2) Site level GAM with shrinkage to global fit and a site level intercept 
localshrinkpath = "SiteSpecificShrinkGamNmix.stan"
shrinkgam = rstan::stan(
  file = paste0(modelpath,localshrinkpath),
  data=birdlist1,
  iter = 1250, # only 250 samples for testing
  warmup = 1000,
  chains=4,
  cores=4,
  init_r = 0.01
)


G2=precis(shrinkgam ,depth=3) 
G2=rownames_to_column(data.frame(G2))
Nestimates = filter(G2, grepl("N",rowname))
resultdtw$Nestshrink= Nestimates$mean

ggplot(resultdtw,aes(year,meanind))+
  geom_point()+
  geom_point(aes(year,Nest),color="red")+
  geom_point(aes(year,Nestshrink),color="blue")+
  geom_line()+
  facet_wrap(~SITECODE)+
  theme_classic()



# 3) Site level GAM with no shrinkage and site level intercept
localfreepath = "SiteSpecificFreeGamNmix.stan"
freegam = rstan::stan(
  file = paste0(modelpath,localfreepath),
  data=birdlist1,
  iter = 1250, # only 250 samples for testing
  warmup = 1000,
  chains=4,
  cores=4,
  init_r = 0.01
)

G3=precis(freegam ,depth=3) 
G3=rownames_to_column(data.frame(G3))
Nestimates = filter(G3, grepl("N",rowname))
resultdtw$Nestfree= Nestimates$mean

ggplot(resultdtw,aes(year,meanind))+
  geom_point()+
  geom_point(aes(year,Nest),color="red")+
  geom_point(aes(year,Nestshrink),color="blue")+
  geom_point(aes(year,Nestfree),color="green")+
  geom_line()+
  facet_wrap(~SITECODE)+
  theme_classic()


# which has lowest WAIC - approximation of cross validation score 
compare(globalgam,shrinkgam,freegam) # shrink gam wins easily - but need to check all validity 



# 4) Try other knot numbers - using the shrinkage gam 
## try other knots
numknots = 6
knotsl=quantile(1:years,probs = seq(0,1,length.out=numknots))
knotsl[-c(1,numknots)] # gets rid of first and last vals
sp6 = bs(1:years,knots=knotsl[-c(1,numknots)],degree=3,intercept = F) # cubic
plot(NULL,xlim=c(1,years),ylim=c(0,1),ylab="basis",xlab="year")
for(i in 1:ncol(sp6)){lines(1:years,sp6[,i])}

# one lower
birdlistlow= birdlist1
birdlistlow$splineMat = sp6
birdlistlow$nsplines = ncol(sp6) 

localshrinkpath = "SiteSpecificShrinkGamNmix.stan"
shrinkgamlow = rstan::stan(
  file = paste0(modelpath,localshrinkpath),
  data=birdlistlow,
  iter = 1250, # only 250 samples for testing
  warmup = 1000,
  chains=4,
  cores=4,
  init_r = 0.01
)


G2low=precis(shrinkgamlow ,depth=3)
G2low=rownames_to_column(data.frame(G2low))
Nestimates = filter(G2low, grepl("N",rowname))
resultdtw$Nestshrinklow= Nestimates$mean

ggplot(resultdtw,aes(year,meanind))+
  geom_point()+
  geom_point(aes(year,Nestshrinklow),color="lightblue")+
  geom_point(aes(year,Nestshrink),color="blue")+
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

birdlisthigh= birdlist1
birdlisthigh$splineMat = sp8
birdlisthigh$nsplines = ncol(sp8) 

shrinkgamhigh = rstan::stan(
  file = paste0(modelpath,localshrinkpath),
  data=birdlisthigh,
  iter = 1250, # only 250 samples for testing
  warmup = 1000,
  chains=4,
  cores=4,
  init_r = 0.01
)


G2high=precis(shrinkgamhigh ,depth=3)
G2high=rownames_to_column(data.frame(G2high))
Nestimates = filter(G2high, grepl("N",rowname))
resultdtw$Nestshrinkhigh= Nestimates$mean

ggplot(resultdtw,aes(year,meanind))+
  geom_point()+
  geom_point(aes(year,Nestshrinklow),color="lightblue")+
  geom_point(aes(year,Nestshrink),color="blue")+
  geom_point(aes(year,Nestshrinkhigh),color="darkblue")+
  geom_line()+
  facet_wrap(~SITECODE)+
  theme_classic()

compare(shrinkgamlow,shrinkgam,shrinkgamhigh)

