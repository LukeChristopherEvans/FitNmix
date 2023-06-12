# ECN data
library(readxl)
library(tidyverse)

folderpath = "~/FitNmix/Data/ECN/"

filepath1 = "ECN_BIRD.csv" # https://catalogue.ceh.ac.uk/documents/5886c3ba-1fa5-49c0-8da8-40e69a10d2b5
filepath2 = "birdsppcode.xlsx"
filepath3 = "sitecode.xlsx"

# bats 
brdt = read.csv(paste0(folderpath,filepath1))
brdt = brdt %>% mutate_all(na_if,"")

# codes 
brsppcode = read_excel(paste0(folderpath,filepath2))
brdt =left_join(brdt,brsppcode,by="FIELDNAME")
rm(brsppcode)
# site codes 
brsitecode = read_excel(paste0(folderpath,filepath3))
brdt =left_join(brdt,brsitecode,by="SITECODE")
rm(brsitecode)

# add on julian day 
library(lubridate)
brdt$SDATE=dmy(brdt$SDATE)
brdt$JDay= yday(brdt$SDATE) # julian day
brdt$year= year(brdt$SDATE) # year

bluetit=filter(brdt, common =="BLUE TIT")

ggplot(bluetit,aes(JDay,VALUE))+
  geom_point()+
  facet_wrap(~SITECODE)

min(bluetit$VALUE) # same structure as bats - missing zeros


# add on all transect and species combinations - I'm taking a transect as one observation effort
brdt=mutate(brdt,obscode=paste(SITECODE,JDay,year,sep="."))

birdinfo = dplyr::select(brdt,FIELDNAME,common)
birdinfo = birdinfo %>% distinct(FIELDNAME,.keep_all = T)
siteinfo = dplyr::select(brdt,SITECODE,SDATE,VISIT,`Site name`, Location,Easting,Northing,Long,Lat,JDay,year,obscode)
siteinfo =  siteinfo %>% distinct(obscode,.keep_all = T)


# loop through - lose things I don't account for here 
brdt2 = brdt %>% dplyr::select(-TRANSECT,-LCODE,-DISTANCE)

# pre allocation should make this a bit quicker
dtlist = vector(mode ="list",length=nrow(siteinfo))

for(i in 1:length(siteinfo$obscode)) {
  siteid = siteinfo$obscode[i]
  td = filter(brdt2,obscode == siteid) # this data
  ts = filter(siteinfo,obscode == siteid) # this site 
  tb = filter(birdinfo,!FIELDNAME %in% td$FIELDNAME) # this bat info
  
  
  newd = cbind.data.frame(ts,tb)
  newd$VALUE=0
  newd = rbind.data.frame(td,newd)
  dtlist[[i]] = newd
}

brdt2=do.call(rbind.data.frame,dtlist)
remove(dtlist,newd,td,tb,ts,siteid,i)

# btdt2, the new dataframe, now has all species coverage 
# I'm going to generate a version of this where I collapse 
# the different transects/locations and group together counts 
brdtwide=brdt2 %>% group_by(obscode,FIELDNAME) %>% 
  summarise(VALUESUM = sum(VALUE)) %>%
  ungroup()%>%
  right_join(brdt2,by=c("obscode","FIELDNAME")) %>%
  group_by(obscode) %>% 
  distinct(FIELDNAME,.keep_all = T) %>%
  ungroup()%>%
  dplyr::select(-VALUE)

# most of these now more similar in size
oby3 = brdtwide %>% group_by(year,SITECODE) %>%
  summarise(n=n())

# investigate high and low values
t4=filter(brdt2,SITECODE=="T01") # T01 sampled only once in 2003 - 133 is num spp
t3=filter(brdt2,SITECODE=="T04") # whereas T04 gets visited 6* in 2002
remove(t4,t3,oby3)

# im now going to make bluetit wide
bluewide=brdtwide %>% filter(FIELDNAME %in% c("BT"))

# now I will pivot to a wide (er) format for the stan models 
bluewide=bluewide %>% group_by(SITECODE,year) %>% mutate(visit=row_number()) %>%ungroup
bluewide=bluewide %>% arrange(SITECODE,year)
# I have to remove the columns with different values across rows, because otherwise I don't get the correct wide format
bluewide2=bluewide %>% dplyr::select(-obscode,-SDATE,-JDay,-VISIT) %>% pivot_wider(names_from=visit,values_from= VALUESUM)
# If I want to send info about visit into the detection model  I can arrange it in the list from bluewide
bluewide2=bluewide2 %>% arrange(SITECODE,year) 
# pad with -1 rather than NA for stan 
bluewide2=replace_na(bluewide2,list('1'=-1,'2'=-1,'3'=-1,'4'=-1,'5'=-1,'6'=-1)) 

# lets quick look at raw data - only from first obs 
ggplot(bluewide2,aes(year,`1`,group=SITECODE))+
  geom_point()+
  geom_line()+
  facet_wrap(~SITECODE)


# I should get means -  because this is the simple
# index 
meanm1 = function(x) mean(x[x != -1])
meanm1(c(-1,2,3,4))
mean(c(2,3,4))

bluewide2$meanind = apply(bluewide2[,11:16],MARGIN =1,FUN=meanm1)

ggplot(bluewide2,aes(year,meanind,group=SITECODE))+
  geom_point()+
  geom_line()+
  facet_wrap(~SITECODE)


# add label for number of visits 
numvis = function(x) length(x[x != -1])
numvis(c(-1,2,3,4))
bluewide2$visits = apply(bluewide2[,11:16],MARGIN =1,FUN=numvis)
# add min -used anom function
bluewide2$minpop = apply(bluewide2[,11:16],MARGIN =1,function(x) min(x[x != -1]))
bluewide2$k = apply(bluewide2[,11:16],MARGIN =1,function(x) max(x[x != -1]))

# save this for analysis - we could put refactor all of the above code in a function for cycling 
# through each species and getting a suitable format for the nmixture models
write.table(bluewide2,paste0(folderpath,"ProcessedBlueTitEcn.txt"))
