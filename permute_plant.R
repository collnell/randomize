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
setwd('/Users/colleennell/Dropbox/Projects/CSS exclusion/')
all.plant<-read.csv('data/2018/CSS_plant_data.csv')
#t.plant<-all.plant%>%filter(treat=='T')
#c.plant<-all.plant%>%filter(treat=='C')%>%mutate(calc = 'ID_C')

ggplot(all.plant, aes(sample, herb_mg_dens))+geom_text(aes(label=sample, color=treat))+facet_wrap(~species)
#lengths by order - sizes<-read.csv('data/2018/CSS_arth_size.csv')
# plant g - methods<-read.csv('data/2017/CSS18_plants.csv')
#methods<-read.csv('data/2017/CSS18_plants.csv')
#str(methods)
#View(methods)

##############################
## find_combos 
# for a plant species (sps, 4-letter code) found in df$species
# find all combinations that samples in treat 'T' can be split into 2 groups, 'DD', and 'ID_T'
# with n_dd samples assigned to the 'DD' group
# Value - a list of dataframes - subset of species' data, all possible combos, combos in melted format, combos with all data

find_combos<-function(sps, df=all.plant, n_dd=3){
  c.plant<-df%>%filter(treat == 'C', species==sps)%>%mutate(calc = 'ID_C') #control data
  t.plant<-df%>%filter(treat=='T', species == sps) #bird exclusion data
  #find al possibl combinations without replacemnt
  comb<-t(combn(as.character(t.plant$sample), n_dd, simplify=TRUE))
  #generate useful dataframes
  sp.comb<-as.data.frame(as.table(comb))%>%dcast(Var1~Var2)%>%mutate(species = paste(sps), calc='DD')
  sp.melt<-sp.comb%>%melt(id.vars=c('species','Var1','calc'))
  sp.cast<-sp.melt%>%dcast(species+value~Var1, value.var='calc', fill = 'ID_T')
  full.df<-t.plant%>%left_join(sp.cast, by=c('species','sample'='value'))
  out=list(t.data=t.plant, c.data=c.plant, combos=sp.comb, combo_melt = sp.melt, combo_cast = sp.cast, combo_all = full.df)
  return(out)
}

# for a single species:
#ARCA<-find_combos('ARCA',all.plant)

# apply to all species
sps.list<-as.character(unique(all.plant$species))
all.combos<-sapply(sps.list, find_combos, USE.NAMES=TRUE, simplify=FALSE) # lapply but retains names from input into list hierarchy
#str(all.combos, max.level=2)

##############################
## combo_vars
# calculate all possible species-level means using combos for any given response var
# LRR bird effect, DD, control - 3 groups
# response - herb density, total biomass, herb/pred, pred density

combo_vars<-function(sps, master=all.combos){
  require(metafor) #for LRR
  combos<-master[[sps]]
  #mean values in control treatment
  c.plant.mean<-combos$c.data%>% #for each species the mean in control treat - used for ID
    group_by(species)%>%
    summarize_at(vars(herb_mg_dens), funs(ID_C_mean = mean(., na.rm=TRUE), ID_C_sd = sd(., na.rm=TRUE), ID_C_se = se, ID_C_n=length))
  
  lrr.df<-combos[['combo_all']]%>%
    melt(id.vars=colnames(combos[['t.data']]), variable.name='combo',value.name='calc')%>% #melt combos
    group_by(species, combo, calc)%>%
    summarize_at(vars(herb_mg_dens), funs(mean(., na.rm=TRUE), sd(., na.rm=TRUE), se, n=length))%>%
    melt(id.vars = c('species','combo','calc'))%>%
    dcast(species+combo~calc+variable)%>% #cast so there are columnsfor each group, row per combo
    left_join(c.plant.mean, by='species') ## add control
  
  #calculate lrr
  lrr.out<-summary(escalc('ROM', m1i=ID_C_mean, m2i=ID_T_mean, sd1i=ID_C_sd, sd2i=ID_T_sd, n1i=ID_C_n, n2i=ID_T_n, 
                          var.names=c('lrr','lrr_var'), digits=3, data=lrr.df))
  
  return(lrr.out)
}  
#ARCA.lrr<-combo_vars(sps='ARCA')
#View(ARCA.lrr$lrr.df)

##all possible lrr & DD for all species
all.lrr<-sapply(sps.list, combo_vars, USE.NAMES=TRUE, simplify=FALSE)
#str(all.lrr, max.level=2)

# combine all lrr.df into single df
lrr.df<-bind_rows(all.lrr)
View(lrr.df)

write.csv(lrr.df, 'data/2018/CSS_all_lrrs.csv', row.names=FALSE)
##############################
## get mean and 95% CI for each species 

ggplot(lrr.df, aes(DD_mean, lrr))+
  geom_errorbar(data=lrr.df, aes(ymin=lrr-sei, ymax=lrr+sei),color='grey')+
  geom_errorbarh(aes(xmin=DD_mean-DD_se, xmax=DD_mean+DD_se), color='grey')+
  geom_point(aes(color=species))

##find mean and 95ci based on distribution of possible values
##95%CI
lrr_mean<-function(lrr.df){
  ##find mean and 95ci based on distribution of possible values
  lrr.avg<-lrr.df%>%group_by(species)%>%
    summarize_at(vars(lrr, DD_mean, ID_C_mean, ID_T_mean), funs(mean, se, n=length))
  lrr.fit<-lm(lrr~1+species, data=lrr.df)
  dd.fit<-lm(DD_mean~1+species, data=lrr.df)
  ci95<-lrr.avg%>%
    left_join(
      as.data.frame(confint(lrr.fit))%>%
        mutate(species=ifelse(rownames(.)=='(Intercept)', 'ARCA', 
                              rownames(.)%>%gsub(pattern='species',replacement=''))))%>%
    mutate(lrr.lb=ifelse(species == 'ARCA',`2.5 %`,`2.5 %`+lrr.fit$coefficients[1]), 
           lrr.ub=ifelse(species == 'ARCA',`97.5 %`, `97.5 %`+lrr.fit$coefficients[1]))%>%
    dplyr::select(-`2.5 %`, -`97.5 %`)%>%
    left_join(
      as.data.frame(confint(dd.fit))%>%
        mutate(species=ifelse(rownames(.)=='(Intercept)', 'ARCA', 
                              rownames(.)%>%gsub(pattern='species',replacement=''))))%>%
    mutate(dd.lb=ifelse(species == 'ARCA',`2.5 %`,`2.5 %`+dd.fit$coefficients[1]), 
           dd.ub=ifelse(species == 'ARCA',`97.5 %`, `97.5 %`+dd.fit$coefficients[1]))
  return(ci95)
}

lrr.summary<-lrr_mean(lrr.df)
#View(lrr.summary)

ggplot(lrr.summary, aes(DD_mean_mean, lrr_mean))+
  geom_errorbar(aes(ymin=lrr.lb, ymax=lrr.ub),color='grey')+
  geom_errorbarh(aes(xmin=dd.lb, xmax=dd.ub), color='grey', height=0)+
  geom_point(aes(color=species))+
  geom_hline(yintercept=0, lty='dashed')

summary(lm(lrr_mean~log(1+DD_mean_mean), data=lrr.summary))

##############################
## plot_lrr
## a function that for a given response variable
#makes figures - all possible DD/ID combos, mean+sd, se, n
#with test significance indicated on the fig
#posthoc contrasts sp*treat - which species differed between treatments?


lrr.hpq<-lrr.df%>%left_join()
ggplot(lrr.df, aes(DD_mean, lrr))+
  geom_errorbar(data=lrr.df, aes(ymin=lrr-sei, ymax=lrr+sei),color='grey')+
  geom_errorbarh(aes(xmin=DD_mean-DD_se, xmax=DD_mean+DD_se), color='grey')+
  geom_point(aes(color=species))



plot_lrr<-function(lrr.df){
  dd_id<-ggplot(all.lrr)
  
}

##bind all lrr.df together
View(all.lrr$ARCA$lrr.df)


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


