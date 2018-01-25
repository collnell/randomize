########################
## combination, randomization, premutation - refined
########################

##############################
library(gtools)
library(dplyr)
library(ggplot2)
library(reshape2)
library(mosaic)

#1. use randomization test to see whether group means differ than null hypothesis (random)
#2. for all permutations or treatment assignments that do not differ from null, calculate LRR and DD - obtain distribution of these values
#3. for each iteratio of LRR and DD, run regression and recrod P, R2, sign

## run same procedure for all species - write as function of species
## apply across list of dataframes for each species
all.plant<-read.csv('/Users/colleennell/Dropbox/Projects/CSS exclusion/data/2018/CSS_plant_groups.csv')%>%dplyr::select(-T_group)
t.plant<-all.plant%>%filter(treat=='T')
c.plant<-all.plant%>%filter(treat=='C')%>%mutate(calc = 'ID_C')

ggplot(all.plant, aes(sample, herb_mg_dens))+geom_text(aes(label=sample, color=treat))+facet_wrap(~species)
#lengths by order - sizes<-read.csv('data/2018/CSS_arth_size.csv')
# plant g - methods<-read.csv('data/2017/CSS18_plants.csv')
methods<-read.csv('data/2017/CSS18_plants.csv')
str(methods)
View(methods)

##############################
## find_combos 
# for a plant species (sps, 4-letter code) found in df$species
# find all combinations that samples in treat 'T' can be split into 2 groups, 'DD', and 'ID_T'
# with n_dd samples assigned to the 'DD' group
# Value - a list of dataframes - subset of species' data, all possible combos, combos in melted format, combos with all data

find_combos<-function(sps, df=t.plant, n_dd=3){
  sp.df<-df%>%filter(species == sps)
  #find al possibl combinations without replacemnt
  comb<-t(combn(as.character(sp.df$sample), n_dd, simplify=TRUE))
  #generate useful dataframes
  sp.comb<-as.data.frame(as.table(comb))%>%dcast(Var1~Var2)%>%mutate(species = paste(sps), calc='DD')
  sp.melt<-sp.comb%>%melt(id.vars=c('species','Var1','calc'))
  sp.cast<-sp.melt%>%dcast(species+value~Var1, value.var='calc', fill = 'ID_T')
  full.df<-sp.df%>%left_join(sp.cast, by=c('species','sample'='value'))
  out=list(t.data=sp.df, combos=sp.comb, combo_melt = sp.melt, combo_cast = sp.cast, combo_all = full.df)
  return(out)
}

# for a single species:
ARCA<-find_combos('ARCA',t.plant)
str(ARCA) 
View(ARCA$combo_cast)
View(ARCA$t.data)
View(ARCA$combo_all)

# apply to all species
sps.list<-as.character(unique(all.plant$species))
str(sps.list)
all.combos<-sapply(sps.list, find_combos, USE.NAMES=TRUE, simplify=FALSE) # lapply but retains names from input into list hierarchy
str(all.combos, max.level=2)
## a function that for a given response variable
#calculate LRR
#makes figures - all possible DD/ID combos, 
#with test significance indicated on the fig
#posthoc contrasts sp*treat - which species differed between treatments?



##############################
#is the difference between the two group means more than expected by random?
ARCA['combo_cast']
13:length(colnames(ARCA$combo_all))
colnames(ARCA$combo_all)[13:length(colnames(ARCA$combo_all))]

#function that tests null for all combo cols
combo_legit<-function(combo_all, var=herb_mg_dens,combos =13:length(colnames(combo_all)), permutations = 1000){
  pval.df<-NULL
  combo_cols<-colnames(combo_all)[13:length(colnames(combo_all))]
                                  
  for(i in combo_cols){
    obs <- diff(mean(var~i, data=combo_all))
    h0<-do(permutations)*diff(mean(var~shuffle(calc), data=combo_all)) 
    res.df<-data.frame(species = paste(unique(combo_all$species)),
                       iter = paste(i),
                       sum = sum(h0$ID_T>=obs),
                       n = length(h0$ID_T),
                       p = sum/length,
                       observed = obs)
    pval.df<-rbind(pval.df, res.df)
  }
  return(pval.df)
}

ARCA_p<-combo_legit(ARCA$combo_all)


  for (j in 4:length(i.iter)){
    i.rand<-i.df%>%left_join(i.iter%>%dplyr::select(calc,sample=j), by='sample')
    obs <- diff(mean(herb_mg_dens~calc, data=i.rand))
    i.h0<-do(permutations)*diff(mean(herb_mg_dens~shuffle(calc), data=i.rand))
    pval<-sum(i.h0$ID_T>=obs)/length(i.h0$ID_T)
    res.df<-data.frame(species = paste(i),
                       iter = paste(colnames(i.iter)[j]),
                       p = pval,
                       sum = sum(i.h0$ID_T>=obs),
                       n = length(i.h0$ID_T),
                       observed = obs)
    i.out<-rbind(res.df, i.out)
  }
}

##match by species and value

for (i in levels(t.plant$species)){ #for each species
  filt.temp<-t.plant%>%filter(species == i)
  temp.cast<-sp.cast%>%filter(species == i)
  temp.cast2<-temp.cast[,colSums(is.na(temp.cast)) != nrow(temp.cast)]
  
  for (j in 4:length(colnames(temp.cast2))){
    id.samps<-setdiff(filt.temp$sample, temp.cast[,j])#which samples are missing from within treat 
    id.df<-data.frame(species = paste(i),
                      Var1 = rep(paste(colnames(temp.cast2)[j]), length(id.samps)),
                      calc = rep('ID_T',length(id.samps)),
                      variable = LETTERS[seq(from=length(temp.cast2[,j])+1, to = length(filt.temp$sample))],
                      value=id.samps)
    sp.melt<-rbind(sp.melt, id.df)#combine with id groupings
  }
}
str(sp.melt)
write.csv(sp.melt, 'rand_end_of_day.csv', row.names=FALSE)


