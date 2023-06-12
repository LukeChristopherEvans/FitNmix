library(tidyverse)

folderpath = "~/FitNmix/Data/RSPB/"
filepath1 = "BBS_MalcolmBurgess_March2022_CSV.csv"

brdt = read.csv(paste0(folderpath,filepath1))

# Column info:	
# ObserverID	Observer unique ID
# Visit	E = Early, L = Late. BBS is a two visit survey, with visits at least 2 weeks apart and both within a set time window
# Species	BT= Blue tit, GT = Great tit, CB = Corn bunting, SL = Barn swallow, GP = Golden plover, S. = Skylark, P. = Grey partridge
# DistanceBand	F = Flying over, 1 = within 25m of transect, 2 = between 25-100m distant from transect, 3 = over 100m distant from transect
# S1-S10	These are the 10 x 200m long sections of the BBS transect

BlueTit=filter(brdt, Species == "BT")

# quick look at data 
BlueTit %>% filter(Gridref =="SK5251") %>%
  ggplot(aes(Year,S1,group=Visit,color=Visit))+
  geom_point()+
  geom_line()

# sum across transects - Could model these separately? Makes sense to sum them for the population at the site 
BlueTit$PopCount = rowSums(BlueTit[,9:18])

BlueTit = BlueTit[,-c(9:18)] # remove s1-10

# add label for a year*site id - called GridEstimate in stan
BlueTit =BlueTit %>% mutate(GridEstimateLabel = paste0(Gridref,Year))

# add an upper and lower bound for the population  # top out at 150 or 3* the observed sample size
BlueTit=BlueTit %>% group_by(GridEstimateLabel) %>%
  mutate(minpop = max(PopCount))


# now I will pivot to a wide(er) format 
BlueTit=BlueTit%>% group_by(GridEstimateLabel) %>% mutate(visitNUM=row_number()) %>%ungroup
BlueTit=BlueTit %>% arrange(GridEstimateLabel)

# change to factor levels 
BlueTit$DistanceBand =  as.integer(as.factor(BlueTit$DistanceBand))
BlueTit$Visit =  as.integer(as.factor(BlueTit$Visit))

# I have to remove the columns with different values across rows, because otherwise I don't get the correct wide format
Bcounts=BlueTit%>% dplyr::select(-ObserverId,-Date,-DistanceBand,-Visit ) %>% pivot_wider(names_from=visitNUM,values_from= PopCount)

Bbands=BlueTit%>% dplyr::select(-ObserverId,-Date,-PopCount,-Visit ) %>% pivot_wider(names_from=visitNUM,values_from= DistanceBand)

Bvisits=BlueTit%>% dplyr::select(-ObserverId,-Date,-PopCount,-DistanceBand ) %>% pivot_wider(names_from=visitNUM,values_from= Visit)

# pad with -1 rather than NA for stan 
Bcounts=replace_na(Bcounts,list('1'=-1,'2'=-1,'3'=-1,'4'=-1,'5'=-1,'6'=-1,'7'=-1,'8'=-1)) # counts
Bbands=replace_na(Bbands,list('1'=-1,'2'=-1,'3'=-1,'4'=-1,'5'=-1,'6'=-1,'7'=-1,'8'=-1)) # bands 
Bvisits=replace_na(Bvisits,list('1'=-1,'2'=-1,'3'=-1,'4'=-1,'5'=-1,'6'=-1,'7'=-1,'8'=-1)) # visits 

numvis = function(x) length(x[x != -1])
Bcounts$visitcount = apply(Bcounts[,9:16],MARGIN =1,FUN=numvis)


colnames(Bbands)[9:16] =c("BandsIndex1","BandsIndex2","BandsIndex3","BandsIndex4","BandsIndex5","BandsIndex6","BandsIndex7","BandsIndex8")
colnames(Bvisits)[9:16] =c("NumVisitsIndex1","NumVisitsIndex2","NumVisitsIndex3","NumVisitsIndex4","NumVisitsIndex5","NumVisitsIndex6","NumVisitsIndex7","NumVisitsIndex8")


BluetitWideComb = cbind.data.frame(Bcounts,Bbands[,9:16],Bvisits[,9:16])

write.table(BluetitWideComb ,paste0(folderpath,"BlueTitWide.txt"))

