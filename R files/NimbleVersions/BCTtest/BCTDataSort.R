library(tidyverse)

folderpath = "/home/lukee/Insync/vs917256@reading.ac.uk/OneDrive Biz/FitNmix/Data/BCT/"
filepath1 = "fieldsurvey.txt"

FieldDat = read.delim(paste0(folderpath,filepath1))

## make wide 
FieldDat =FieldDat %>% mutate(GridEstimateLabel = paste0(GridReference,CountYear))

# get rid of empties
FieldDat = FieldDat %>% mutate_all(na_if,"")
# make detector an index
FieldDat$DetectorInd = as.integer(as.factor(FieldDat$Detector)) # bat detector 

# add an upper and lower bound for the population  
# need to split by species 
colnames(FieldDat)
FNoctule = FieldDat[,c(1:20,25,29,30)]
FSerotine = FieldDat[,c(1:19,21,25,29,30)]
FLeisler = FieldDat[,c(1:19,23,25,29,30)]
FCommonPip = FieldDat[,c(1:19,26,25,29,30)]
FSopranoPip = FieldDat[,c(1:19,27,25,29,30)]

# I need to replace the empties with zero in the count columns 
FNoctule=FNoctule %>% dplyr::mutate(SumOfNoctuleCount = replace_na(SumOfNoctuleCount, 0))
FSerotine=FSerotine  %>% dplyr::mutate(SumOfSerotineCount = replace_na(SumOfSerotineCount, 0))
FLeisler=FLeisler %>% dplyr::mutate(SumOfLeislerCount = replace_na(SumOfLeislerCount, 0))
FCommonPip=FCommonPip %>% dplyr::mutate(SumOfCommonPipCount= replace_na(SumOfCommonPipCount, 0))
FSopranoPip=FSopranoPip %>% dplyr::mutate(SumOfSopranoPipCount = replace_na(SumOfSopranoPipCount, 0))


FNoctule=FNoctule %>% group_by(GridEstimateLabel) %>%
  mutate( minpop = max(SumOfNoctuleCount)) 
FSerotine=FSerotine %>% group_by(GridEstimateLabel) %>%
  mutate(minpop = max(SumOfSerotineCount)) 
FLeisler=FLeisler %>% group_by(GridEstimateLabel) %>%
  mutate( minpop = max(SumOfLeislerCount)) 
FCommonPip=FCommonPip %>% group_by(GridEstimateLabel) %>%
  mutate(minpop = max(SumOfCommonPipCount)) 
FSopranoPip=FSopranoPip %>% group_by(GridEstimateLabel) %>%
  mutate(minpop = max(SumOfSopranoPipCount)) 


FNoctule=FNoctule%>% group_by(GridEstimateLabel) %>% mutate(visitNUM=row_number()) %>%ungroup
FSerotine=FSerotine%>% group_by(GridEstimateLabel) %>% mutate(visitNUM=row_number()) %>%ungroup
FLeisler=FLeisler%>% group_by(GridEstimateLabel) %>% mutate(visitNUM=row_number()) %>%ungroup
FCommonPip=FCommonPip%>% group_by(GridEstimateLabel) %>% mutate(visitNUM=row_number()) %>%ungroup
FSopranoPip=FSopranoPip%>% group_by(GridEstimateLabel) %>% mutate(visitNUM=row_number()) %>%ungroup

FNoctule=FNoctule%>% arrange(GridEstimateLabel)
FSerotine=FSerotine%>% arrange(GridEstimateLabel)
FLeisler=FLeisler%>% arrange(GridEstimateLabel)
FCommonPip=FCommonPip%>% arrange(GridEstimateLabel)
FSopranoPip=FSopranoPip%>% arrange(GridEstimateLabel)


## wide version of counts - the ecological component 
WFNoctule=FNoctule[,-c(6:19,21,23)]   %>% pivot_wider(names_from=visitNUM,values_from= SumOfNoctuleCount)
WFSerotine=FSerotine[,-c(6:19,21,23)]   %>% pivot_wider(names_from=visitNUM,values_from= SumOfSerotineCount)
WFLeisler=FLeisler[,-c(6:19,21,23)]   %>% pivot_wider(names_from=visitNUM,values_from= SumOfLeislerCount)
WFCommonPip=FCommonPip[,-c(6:19,21,23)]   %>% pivot_wider(names_from=visitNUM,values_from= SumOfCommonPipCount)
WFSopranoPip=FSopranoPip[,-c(6:19,21,23)]   %>% pivot_wider(names_from=visitNUM,values_from= SumOfSopranoPipCount)

# detection components 
## wide detector model - species doesn't matter because the same visit 
WDetector=FNoctule[,-c(6:21)]   %>% pivot_wider(names_from=visitNUM,values_from= DetectorInd)
## wide number of spots
WSpots=FNoctule[,-c(6:20,23)]   %>% pivot_wider(names_from=visitNUM,values_from= NumberOfSpotsSurveyed)
## wide number of walks 
WWalks=FNoctule[,-c(6:18,20,21,23)]   %>% pivot_wider(names_from=visitNUM,values_from= NumberOfWalksSurveyed)


# remove nas for stan
WFNoctule=replace_na(WFNoctule,list('1'=-1,'2'=-1,'3'=-1,'4'=-1)) 
WFSerotine=replace_na(WFSerotine,list('1'=-1,'2'=-1,'3'=-1,'4'=-1)) 
WFLeisler=replace_na(WFLeisler,list('1'=-1,'2'=-1,'3'=-1,'4'=-1)) 
WFCommonPip=replace_na(WFCommonPip,list('1'=-1,'2'=-1,'3'=-1,'4'=-1)) 
WFSopranoPip=replace_na(WFSopranoPip,list('1'=-1,'2'=-1,'3'=-1,'4'=-1)) 

WDetector=replace_na(WDetector,list('1'=-1,'2'=-1,'3'=-1,'4'=-1)) 
WSpots=replace_na(WSpots,list('1'=-1,'2'=-1,'3'=-1,'4'=-1)) 
WWalks=replace_na(WWalks,list('1'=-1,'2'=-1,'3'=-1,'4'=-1)) 

numvis = function(x) length(x[x != -1])
WFNoctule$visitcount = apply(WFNoctule[,8:11],MARGIN =1,FUN=numvis)
WFSerotine$visitcount = apply(WFSerotine[,8:11],MARGIN =1,FUN=numvis)
WFLeisler$visitcount = apply(WFLeisler[,8:11],MARGIN =1,FUN=numvis)
WFCommonPip$visitcount = apply(WFCommonPip[,8:11],MARGIN =1,FUN=numvis)
WFSopranoPip$visitcount = apply(WFSopranoPip[,8:11],MARGIN =1,FUN=numvis)


# add different column names 
colnames(WDetector)[8:11] =c("DetectorIndex1","DetectorIndex2","DetectorIndex3","DetectorIndex4")
colnames(WSpots)[8:11] =c("NumSpotsIndex1","NumSpotsIndex2","NumSpotsIndex3","NumSpotsIndex4")
colnames(WWalks)[8:11] =c("NumWalksIndex1","NumWalksIndex2","NumWalksIndex3","NumWalksIndex4")


# join
WFNoctule=cbind.data.frame(WFNoctule,WDetector[,8:11],WSpots[,8:11],WWalks[,8:11])
WFSerotine=cbind.data.frame(WFSerotine,WDetector[,8:11],WSpots[,8:11],WWalks[,8:11])
WFLeisler=cbind.data.frame(WFLeisler,WDetector[,8:11],WSpots[,8:11],WWalks[,8:11])
WFCommonPip=cbind.data.frame(WFCommonPip,WDetector[,8:11],WSpots[,8:11],WWalks[,8:11])
WFSopranoPip=cbind.data.frame(WFSopranoPip,WDetector[,8:11],WSpots[,8:11],WWalks[,8:11])


WFNoctule=arrange(WFNoctule,SiteCode,CountYear)
WFSerotine=arrange(WFSerotine,SiteCode,CountYear)
WFLeisler=arrange(WFLeisler,SiteCode,CountYear)
WFCommonPip=arrange(WFCommonPip,SiteCode,CountYear)
WFSopranoPip=arrange(WFSopranoPip,SiteCode,CountYear)

# add detector label back on 

detectorcode = distinct(FieldDat[,c(10,30)])

# I remove this site because its got mixed up - 
WFCommonPip = filter(WFCommonPip,SiteCode !="120614")
WFNoctule = filter(WFNoctule,SiteCode !="120614")
WFSerotine = filter(WFSerotine,SiteCode !="120614")
WFLeisler = filter(WFLeisler,SiteCode !="120614")
WFSopranoPip= filter(WFSopranoPip,SiteCode !="120614")

# save
write.table(WFCommonPip,paste0(folderpath,"FieldCommonPipWide.txt"))
write.table(detectorcode ,paste0(folderpath,"Detectorcode.txt"))
