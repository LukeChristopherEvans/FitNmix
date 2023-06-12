# ECN data
library(readxl)
library(tidyverse)

folderpath = "~FitNmix/Data/ECN/"

filepath1 = "ECN_BAT.csv" #https://catalogue.ceh.ac.uk/documents/2588ee91-6cbd-4888-86fc-81858d1bf085
filepath2 = "batsppcode.xlsx" #see readme from above to find species and site codes
filepath3 = "sitecode.xlsx"

btdt = read.csv(paste0(folderpath,filepath1))
btdt = btdt %>% mutate_all(na_if,"")

# codes 
btsppcode = read_excel(paste0(folderpath, filepath2))
btdt =left_join(btdt,btsppcode,by="FIELDNAME")
rm(btsppcode)

# site codes 
btsitecode = read_excel(paste0(folderpath,filepath3))
btdt =left_join(btdt,btsitecode,by="SITECODE")
rm(btsitecode)

# add on julian day 
library(lubridate)
btdt$SDATE=dmy(btdt$SDATE)
btdt$JDay= yday(btdt$SDATE) # julian day
btdt$year= year(btdt$SDATE) # year

bpip=filter(btdt, FIELDNAME =="Ppl")

ggplot(bpip,aes(JDay,VALUE))+
  geom_point()+
  facet_wrap(~SITECODE)


library(maps)
library(mapdata)
Ukm<-map('worldHires','UK',	xlim=c(-11,3.2), ylim=c(49,60.9))
par(mar=rep(0,4))
points(bpip$Long,bpip$Lat,pch=19,col="red")

# they only record ones but not the zeros, I'm going to add on the zeros
# I'm taking a whole transect as one observation effort
btdt=mutate(btdt,obscode=paste(SITECODE,JDay,year,sep="."))

batinfo = dplyr::select(btdt,FIELDNAME,latin,common)
batinfo = batinfo %>% distinct(FIELDNAME,.keep_all = T)
siteinfo = dplyr::select(btdt,SITECODE,SDATE,ACTS,ACTH,ACTF,`Site name`, Location,Easting,Northing,Long,Lat,JDay,year,obscode)
siteinfo =  siteinfo %>% distinct(obscode,.keep_all = T)

# loop through 
btdt2 = btdt %>% dplyr::select(-TRANSECT,-BATLOC_ID)

# pre allocation should make this a bit quicker
dtlist = vector(mode ="list",length=nrow(siteinfo))

for(i in 1:length(siteinfo$obscode)) {
  siteid = siteinfo$obscode[i]
  td = filter(btdt2,obscode == siteid) # this data
  ts = filter(siteinfo,obscode == siteid) # this site 
  tb = filter(batinfo,!FIELDNAME %in% td$FIELDNAME) # this bat info
  
  
  newd = cbind.data.frame(ts,tb)
  newd$VALUE=0
  newd = rbind.data.frame(td,newd)
  dtlist[[i]] = newd
}

btdt2=do.call(rbind.data.frame,dtlist)
remove(dtlist,newd,td,tb,ts,siteid,i)


# btdt2, the new dataframe, now has the zeros in
# I'm going to generate a version of this where I collapse 
# across the different transects per visit and group together all the counts 
btdtwide=btdt2 %>% group_by(obscode,FIELDNAME) %>% 
  summarise(VALUESUM = sum(VALUE)) %>%
  ungroup()%>%
  right_join(btdt2,by=c("obscode","FIELDNAME")) %>%
  group_by(obscode) %>% 
  distinct(FIELDNAME,.keep_all = T) %>%
  ungroup()%>%
  dplyr::select(-VALUE)


# most of these now more similar in size
oby3 = btdtwide %>% group_by(year,SITECODE) %>%
  summarise(n=n())

# investigate high and low values
t4=filter(btdt2,SITECODE=="T04") # So in particular years, this site just
# get sampled ~ 8 times
t3=filter(btdt2,SITECODE=="T03") # whereas T03 gets visited once in 2000
remove(t4,t3,oby3)

# I'm now going to make a wide version for the pipistrelle 
# and I'm going to group together the counts for both common and pip 
# we would want to be careful about the pips that could be soprano or 
# common in a real analysis! - Maybe could put an uncertainty on just pip 
pipwide=btdtwide %>% filter(FIELDNAME %in% c("Pp","Ppl")) %>%
  group_by(obscode) %>%
  mutate(VALUESUM = sum(VALUESUM)) %>%
  ungroup() %>%
  distinct(obscode,.keep_all = T) %>%
  mutate(FIELDNAME ="ToP") # total pip


# now I will pivot to a wide (er) format for the bayesian models 
pipwide=pipwide %>% group_by(SITECODE,year) %>% mutate(visit=row_number()) %>%ungroup
pipwide=pipwide %>% arrange(SITECODE,year)
# I have to remove the columns with different values across rows, because otherwise I don't get the correct wide format
pipwide2=pipwide %>% dplyr::select(-obscode,-ACTS,-ACTH,-ACTF,-SDATE,-JDay) %>% pivot_wider(names_from=visit,values_from= VALUESUM)
# If I want to send info about detection into the detection model  I can arrange it in the list from pipwide
pipwide2=pipwide2 %>% arrange(SITECODE,year) 
# pad with -1 rather than NA for stan 
pipwide2=replace_na(pipwide2,list('1'=-1,'2'=-1,'3'=-1,'4'=-1,'5'=-1,'6'=-1,'7'=-1,'8'=-1)) 

# lets quick look at raw data - only from first obs 
ggplot(pipwide2,aes(year,`1`,group=SITECODE))+
  geom_point()+
  geom_line()+
  facet_wrap(~SITECODE)

# I should get means -  because this is the simple
# index 
meanm1 = function(x) mean(x[x != -1])
meanm1(c(-1,2,3,4))
mean(c(2,3,4))

pipwide2$meanind = apply(pipwide2[,12:19],MARGIN =1,FUN=meanm1)

ggplot(pipwide2,aes(year,meanind,group=SITECODE))+
  geom_point()+
  geom_line()+
  facet_wrap(~SITECODE)


# add label for number of visits 
numvis = function(x) length(x[x != -1])
numvis(c(-1,2,3,4))
pipwide2$visits = apply(pipwide2[,12:19],MARGIN =1,FUN=numvis)
# add min -used anom function
pipwide2$minpop = apply(pipwide2[,12:19],MARGIN =1,function(x) min(x[x != -1]))
pipwide2$k = apply(pipwide2[,12:19],MARGIN =1,function(x) max(x[x != -1]))


# save this for analysis - we could put refactor all of the above code in a function for cycling 
# through each species and getting a suitable format for the nmixture models
write.table(pipwide2,paste0(folderpath,"ProcessedPipEcn.txt"))
