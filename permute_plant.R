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

##############################
#find all combinations can be made by splitting 'T' samples into 2 groups - ID_T (n=5), DD (n=3)
find_combos<-function(sps, df, n_dd=3){
  sp.df<-df%>%filter(species == sps)
  #find al possibl combinations without replacemnt
  comb<-t(combn(as.character(sp.df$sample), n_dd, simplify=TRUE))
  #generate useful dataframes
  sp.comb<-as.data.frame(as.table(comb))%>%dcast(Var1~Var2)%>%mutate(species = paste(sps), calc='DD')
  sp.melt<-sp.comb%>%melt(id.vars=c('species','Var1','calc'))
  sp.cast<-sp.melt%>%dcast(species+value~Var1, value.var='calc', fill = 'ID_T')
  full.df<-sp.df%>%left_join(sp.cast, by=c('species','sample'='value'))
  return(list(t.data=sp.df, combos=sp.comb, combo_melt = sp.melt, combo_cast = sp.cast, combo_all = full.df))
}
ARCA<-find_combos('ARCA',t.plant)
str(ARCA) # returns list with 2 df - the raw t.data, possible combos, melted combos, and cast cmbos (col each combo)
View(ARCA$combo_cast)
View(ARCA$t.data)
View(ARCA$combo_all)

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


