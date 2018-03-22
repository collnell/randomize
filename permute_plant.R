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
#3. for each iteratio of LRR and DD

## run same procedure for all species - write as function of species
## apply across list of dataframes for each species
setwd('/Users/colleennell/Dropbox/Projects/CSS exclusion/')

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
##############################
## get mean and 95% CI for each species 
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
# for a single species:
#ARCA<-find_combos('ARCA',all.plant)

##############################
all.plant<-read.csv('data/2018/CSS_plant_data.csv')
str(all.plant)
# do same for predators
#all.plant$herb_mg_dens<-all.plant$pred_mg_dens
ggplot(all.plant, aes(sample, herb_mg_dens))+geom_text(aes(label=sample, color=treat))+facet_wrap(~species)

# apply to all species
sps.list<-as.character(unique(all.plant$species))
all.combos<-sapply(sps.list, find_combos, USE.NAMES=TRUE, simplify=FALSE) # lapply but retains names from input into list hierarchy
#str(all.combos, max.level=2)

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


##all possible lrr & DD for all species
all.lrr<-sapply(sps.list, combo_vars, USE.NAMES=TRUE, simplify=FALSE)
#str(all.lrr, max.level=2)

# combine all lrr.df into single df
lrr.df<-bind_rows(all.lrr)
#View(lrr.df)

lrr.summary<-lrr_mean(lrr.df)
#View(lrr.summary)

ggplot(lrr.df, aes(DD_mean, lrr))+
  geom_errorbar(data=lrr.df, aes(ymin=lrr-sei, ymax=lrr+sei),color='grey')+
  geom_errorbarh(aes(xmin=DD_mean-DD_se, xmax=DD_mean+DD_se), color='grey')+
  geom_point(aes(color=species))+
  geom_hline(yintercept=0, lty='dotted')+
  labs(x='Herbivore density in exclusion',y='LRR Bird effect')

#write.csv(lrr.df, 'data/2018/CSS_all_lrrs.csv', row.names=FALSE)

ggplot(lrr.summary, aes(DD_mean_mean, lrr_mean))+
  geom_errorbar(aes(ymin=lrr.lb, ymax=lrr.ub),color='grey')+
  geom_errorbarh(aes(xmin=dd.lb, xmax=dd.ub), color='grey', height=0)+
  geom_point(aes(color=species))+
  geom_smooth(method='lm', se=F, lty='dashed', color='grey')+
  geom_hline(yintercept=0, lty='dotted')+
  labs(x='Herbivore density in exclusion',y='LRR Bird effect')

summary(lm(lrr_mean~DD_mean_mean, data=lrr.summary))
# P = 0.1127
##############################
# test for proportionality
log.x<-lrr.summary$lrr_mean
log.y<-log(lrr.summary$DD_mean_mean)
log.xy<-log.x+log.y
log.x.on.y<-log.x-log.y
tstat<-cor(log.xy, log.x.on.y)^2
permutations=10000
p.value<-sum(replicate(permutations, cor(log.xy, sample(log.x.on.y))^2>tstat))/permutations
p.value ##0.4345
#do not reject the null hypothesis that the slope is 1 - DD and LRR are proportiona;

library(propr)

# phi and rho
# phit, perb, phis
?perb
perb(lrr.df, ivar=0)

##############################
## is this different than just the raw means without groups?
raw.mean<-all.plant%>%
  group_by(species, treat)%>%
  summarize_at(vars(herb_mg_dens), funs(mean(., na.rm=TRUE), sd(., na.rm=TRUE), se, length))%>%
  melt(id.vars=c('species','treat'))%>%
  dcast(species~treat+variable)

raw.lrr<-summary(escalc('ROM', m1i=C_mean, m2i=T_mean, sd1i=C_sd, sd2i=T_sd, n1i=C_length, n2i=T_length, data=raw.mean))

ggplot(raw.lrr, aes(T_mean, yi))+
  geom_point()+
  geom_smooth(method='lm', se=FALSE)+
  geom_hline(yintercept=0, lty='dotted')

summary(lm(yi~log(1+T_mean), data=raw.lrr))

##############################
# select 1 combo from each species
# test species' correlations
samp_res<-function(df){
  # select 1 combo from each species
  samp_combo<-df%>%group_by(species)%>%sample_n(size=1)
  # test correlation
  samp.lm<-lm(samp_combo$lrr~log(1+samp_combo$DD_mean))
  # p, r2, est, tval - make a df with results
  samp.df<-data.frame(rsq = summary(samp.lm)$r.squared, pval = summary(samp.lm)$coefficients[2,4], 
                      tval = summary(samp.lm)$coefficients[2,3],est = summary(samp.lm)$coefficients[2,1])
}

##takes a while..s
samp_reps<-do(10000)*samp_res(lrr.df)
#n/total = p prop
samp_reps$sig<-ifelse(samp_reps$pval<=0.05, 1, 0)
sum(samp_reps$sig)/length(samp_reps$sig)
## p =0.1189 
mean(samp_reps$rsq) ##0.2189

histogram(samp_reps$rsq)
mean(samp_reps$rsq)
sum(samp_reps$sig) #over 20% of them are significant
histogram(samp_reps$est)
# if the correlation is significant less than 5% of the time...

##OR - get null distribution of correlation when species is random - shuffle 
###########################
## hpq
ggplot(raw.lrr, aes(C_mean, T_mean))+
  geom_point()+
  geom_smooth(method='lm', se=F)

hpq<-read.csv('data/2017/CSS_means.csv')%>%dplyr::select(sp, comp, comp_se, hpq, hpq_se)

hpt<-raw.lrr%>%
  left_join(hpq, by=c('species'='sp'))

summary(aov(lm(yi~log(1+hpq), data=hpt)))
summary(aov(lm(T_mean~log(1+hpq), data=hpt)))
## hpq and herb dens in exclusion is correlated

ggplot(hpt, aes(T_mean, hpq))+
  geom_point()+
  geom_smooth(method='lm', se=F)+
  labs(x='Host plant quality', y='Herbivore density in exclusion')


rlong<-raw.lrr%>%dplyr::select(species, C_mean, T_mean)%>%
  melt(id.vars=c('species'))%>%
  left_join(hpq, by=c('species'='sp'))

#KM - usedlsmean for bird predation for each species controlling for density
#LRR~species+density
#density controls for association between density and bird predation as well as nonindependence
Anova(lm(yi~T_mean*log(1+hpq), data=hpt), type='III')
#only T_mean is sig

levels(rlong$variable)<-c('birds','no birds')
ggplot(rlong, aes(hpq, value))+geom_point(aes(color=variable))+
  geom_smooth(method='lm', se=F, aes(color=variable))+
  scale_x_log10()+
  labs(x='Host plant quality\n(caterpillar weight gain)', y='Herbivore density')

Anova(lm(value~log(1+hpq)*variable, data=rlong), type='III')

## hpq correlated with herb dens


##############################
## account for phylogenetic relatedness (phylogenetic non-independence) in trait-trait or trait-habitat correlations


lrr.hpq<-lrr.df%>%left_join()

ggplot(lrr.df, aes(DD_mean, lrr))+
  geom_errorbar(data=lrr.df, aes(ymin=lrr-sei, ymax=lrr+sei),color='grey')+
  geom_errorbarh(aes(xmin=DD_mean-DD_se, xmax=DD_mean+DD_se), color='grey')+
  geom_point(aes(color=species))


##############################
## dealing with spurrious correlation
library(propr)#https://github.com/tpq/propr/blob/master/vignettes/a_introduction.Rmd
# phi
# rho
#phi_s

# centered log-ratio

# non-linear regression for ratio variables



# test null hypthesis - https://www.nature.com/articles/srep23247
# appropriate null
# distribution of correlation coefficient between metrics
# standard null - coefficient of correlation is 0 (rx,x-y) - erroneous (range would be between 1 and 01)
# pearson correlation between change (x-y) and pretreatment value (x) 
# if changes in variables are relateed to baseline values, consequential changes in variances of these vars
# 1. determine range of the correlation coefficient between LRR and DD, conditional on the correlation between x any y

# cor coef
arca.lrr<-lrr.df%>%filter(species =='ARCA')
cor(arca.lrr$DD_mean, arca.lrr$lrr) #.99

# mean spurious coefficient of determination r2
# sample correlation coefficient varies a lot depending on what data subset is used

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


