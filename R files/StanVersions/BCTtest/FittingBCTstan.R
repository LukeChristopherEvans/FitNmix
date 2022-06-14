library(tidyverse)
library(splines)
library(rstan)
library(rethinking)
library(coda)

folderpath = "/home/lukee/Insync/vs917256@reading.ac.uk/OneDrive Biz/FitNmix/Data/BCT/"
modelpath = "/home/lukee/Insync/vs917256@reading.ac.uk/OneDrive Biz/FitNmix/Stan files/Bat/"
filepath1 = "FieldCommonPipWide.txt"


# !Note! we only 'profile' the models here to save computation time, and we, therefore, only avoid divergent transitions 
# as fatal errors for interpretation. For a proper analysis we would need more samples from the posterior
# and to spend more time evaluating chain mixing

# The plan here is to fit the model across 100 sites
# A change from the previous model is that we will have detection changed per site 
# and we have parameters on the detector model 

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

#new way of deciding k 
btdt$K = btdt$minpop +  qpois(0.99999,btdt$minpop)+3

Globallist = list(
  # data lengths
  n = nrow(btdt),
  nsite = length(unique(btdt$SiteCode)),
  nyear = length(unique(btdt$CountYear)),
  nsplines = ncol(sp),

  nDetectors = length(1:max(btdt[,13:16])),
  
  maxvisits = 4,
  
  visits =  btdt$visitcount,
  
  # population data
  k = btdt$K,
  minpop = btdt$minpop,
  count = as.matrix(btdt[,8:11]),
  
  site = as.integer(as.factor(btdt$SiteCode)),
  year = btdt$CountYear - min(btdt$CountYear) + 1,
  
  # observation data
  Detector = ifelse(matrix(as.integer(as.factor(as.matrix(btdt[,13:16]))) ,ncol=4)-1==0,-1,matrix(as.integer(as.factor(as.matrix(btdt[,13:16]))) ,ncol=4)-1),
  SpotSections = as.matrix(btdt[,17:20]),
  
  splineMat = sp
)



# get the file path to the stan model 
filepath = "GlobalGamNmixMultiProbBat.stan"

# this model will tell you off a lot because p is a ragged array -  key thing 
# to watch out for at this stage is divergent transitions. Note we also have some high 
# rhats that would need investigating in a proper analysis
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
resultdtw = btdt # make a copy of original data 
resultdtw$Nestglobal= Nestimates$mean

meanind = function(x) mean(x[x!=-1])
sumind = function(x) sum(x[x!=-1])
avg=c()
sums=c()
for(i in 1:nrow(resultdtw)) {
  avg=c(avg,meanind(btdt[i,8:11]))
  sums=c(sums,sumind(btdt[i,8:11]))
}

# add simple
resultdtw$avg= avg
resultdtw$sum= sums

filter(resultdtw, SiteCode %in% sample(unique(resultdtw$SiteCode),20)) %>%
  ggplot(aes(CountYear,avg))+
  geom_point()+
  geom_point(aes(CountYear,Nestglobal),color="red")+
  facet_wrap(~SiteCode)+
  theme_classic()

# Now I can look at the observational model 
plot(globalgam,pars=c("delta0","delta2","delta1"))


##
# 2) Site level GAM with shrinkage to global fit and a site level intercept - rhats better here
localshrinkpath = "SiteSpecificShrinkGamNmixMultiProbBat.stan"
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

filter(resultdtw, SiteCode %in% sample(unique(resultdtw$SiteCode),20)) %>%
  ggplot(aes(CountYear,avg))+
  geom_point()+
  geom_point(aes(CountYear,Nestglobal),color="red")+
  geom_point(aes(CountYear,Nestshrink),color="blue")+
  facet_wrap(~SiteCode)+
  theme_classic()

plot(shrinkgam,pars=c("delta0","delta2","delta1"))

# 3) Site level GAM with no shrinkage and site level intercept - high rhats again
localfreepath = "SiteSpecificFreeGamNmixMultiProbBat.stan"
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

filter(resultdtw, SiteCode %in% sample(unique(resultdtw$SiteCode),20)) %>%
  ggplot(aes(CountYear,avg))+
  geom_point()+
  geom_point(aes(CountYear,Nestglobal),color="red")+
  geom_point(aes(CountYear,Nestshrink),color="blue")+
  geom_point(aes(CountYear,Nestfree),color="green")+
  facet_wrap(~SiteCode)+
  theme_classic()


# which has lowest WAIC - approximation of cross validation score 
compare(globalgam,shrinkgam,freegam) # however, note validity

# 4) Try other knot numbers - using the shrinkage gam 
## try other knots
numknots = 6
knotsl=quantile(1:years,probs = seq(0,1,length.out=numknots))
knotsl[-c(1,numknots)] # gets rid of first and last vals
sp6 = bs(1:years,knots=knotsl[-c(1,numknots)],degree=3,intercept = F) # cubic
plot(NULL,xlim=c(1,years),ylim=c(0,1),ylab="basis",xlab="year")
for(i in 1:ncol(sp6)){lines(1:years,sp6[,i])}

# one lower
batlistlow= Globallist
batlistlow$splineMat = sp6
batlistlow$nsplines = ncol(sp6) 

localshrinkpath = "SiteSpecificShrinkGamNmixBat.stan"
shrinkgamlow = rstan::stan(
  file = paste0(modelpath,localshrinkpath),
  data=batlistlow,
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

batlisthigh= Globallist
batlisthigh$splineMat = sp8
batlisthigh$nsplines = ncol(sp8) 

shrinkgamhigh = rstan::stan(
  file = paste0(modelpath,localshrinkpath),
  data=batlisthigh,
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

ggplot(resultdtw,aes(year,meanind))+
  geom_point()+
  geom_point(aes(year,Nestshrinklow),color="lightblue")+
  geom_point(aes(year,Nestshrink),color="blue")+
  geom_point(aes(year,Nestshrinkhigh),color="darkblue")+
  geom_line()+
  facet_wrap(~SITECODE)+
  theme_classic()

compare(shrinkgamlow,shrinkgam,shrinkgamhigh)
