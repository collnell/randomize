########################
## combination, randomization, premutation
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
#4. then the results are....?

##all possible combinations of gorup assignments within arca..to start
all.plant<-read.csv('/Users/colleennell/Dropbox/Projects/CSS exclusion/data/2018/CSS_plant_groups.csv')
t.plant<-all.plant%>%filter(treat=='T')%>%mutate(group = ifelse(T_group == 'ID_T', 'T_ID', ifelse(T_group == 'DD', 'T_DD', NA)))%>%transform(group=as.factor(group))
arca<-t.plant%>%filter(species=='ARCA')

##possible combinations of group such that there is 4 in each group
cob<-t(combn(as.character(arca$sample), 4, simplify=TRUE))
cob

arca.comb<-as.data.frame(as.table(cob))%>%dcast(Var1~Var2)%>%mutate(calc='DD')##all possible combos within a single treatment, var1 is unique combo
arca.melt<-arca.comb%>%melt(id.vars=c('Var1','calc'))
arca.cast<-arca.melt%>%dcast(variable+calc~Var1)
View(arca.cast)

##so now each column isa combination, 35 columns for all possible combos
#calc is the assigned calculation group
##needs the ID group now for each possible combo
for (i in 3:length(colnames(arca.cast))){#for each combo
  id.samps<-setdiff(arca$sample, arca.cast[,i])#which samples are missing from within treat- arca$sample
  ##new df with samples 'value' paired to the combination (Var1), assigned calc var (ID), varible is simply unique for each new value
  id.df<-data.frame(Var1 = rep(paste(colnames(arca.cast)[i]), length(id.samps)),
                    calc = rep('ID_T',length(id.samps)),
                    variable = LETTERS[seq(from=length(arca.cast[,i])+1, to = length(arca$sample))],
                    value=id.samps)
  arca.melt<-rbind(arca.melt, id.df)#combine with id groupings
}
arca.new<-arca.melt
View(arca.new)

##now cast arca.melt so each colum is a combination with the treatments assigned
arca.iter<-arca.new%>%dcast(variable+calc~Var1)
View(arca.iter)

###################################

##for each column in arca.iter
#bind with arca (density) 
arca.ran<-arca%>%left_join(arca.iter%>%dplyr::select(calc,sample=A), by='sample')
View(arca.ran)
#use randomization test - shuffle 'calc' - are differences in means due to chance or not?
#distribution of the test statistic uneder the null hypothesis iobtained by calculating all possible values of tstat under randomized groupings
#observed difference in means
obs <- diff(mean(herb_mg_dens~calc, data=arca.ran))
obs
arca.h0<-do(1000)*diff(mean(herb_mg_dens~shuffle(calc), data=arca.ran))
arca.h0
confint(arca.h0, level=0.975, method="quantile")#confidence intervals..of the difference?
histogram(~ ID_T, groups=(ID_T >= obs), data=arca.h0, width=0.001,xlab="Distribution of difference in means\nunder the null hypothesis")

#legit?
sum(arca.h0$ID_T>=obs)/length(arca.h0$ID_T)##P = 0.951...so this grouping combination does not reject null - should be more reliable

##do this FOR ALL COMBOS - can use same null
#is the difference between the two group means more than expected by random?
arca.out<-NULL
for (i in 3:length(arca.iter)){
  arca.ran<-arca%>%left_join(arca.iter%>%dplyr::select(calc,sample=i), by='sample')
  obs <- diff(mean(herb_mg_dens~calc, data=arca.ran))
  arca.h0<-do(1000)*diff(mean(herb_mg_dens~shuffle(calc), data=arca.ran))
  pval<-sum(arca.h0$ID_T>=obs)/length(arca.h0$ID_T)
  res.df<-data.frame(iter = paste(colnames(arca.iter)[i]),
                     p = pval,
                     sum = sum(arca.h0$ID_T>=obs),
                     n = length(arca.h0$ID_T),
                     observed = obs)
  arca.out<-rbind(res.df, arca.out)
}
# takes a bit (3 min)...code is super clunky
View(arca.out)
arca.tested<-arca.out%>%filter(p < 0.05)
arca.tested$iter
#drop iter T from 
arca.iter.filt<-arca.iter[,-which(names(arca.iter) %in% arca.tested$iter)]

#MCPP - monte carlo permutation procedure
# calcualte tstat, compare to true, repeat
#if tstat is greater than 95% of random values, can reject null at P<0.025 (two.tailed test)
##########################################

##now using combinations that do not reject null
#iterate through each column to calculate the LRR and mean DD +sd, se, n

#control treatment plants to combine with other df
arca.c<-all.plant%>%filter(treat=='C' & species =='ARCA')%>%dplyr::select(-T_group)%>%mutate(calc = 'ID_C')
#use by= to iterate through columns, and combine with the cacl col?

##now iterate that through all possible treatment combinations
iter.means<-NULL
for (i in 3:length(colnames(arca.iter.filt))){
  arca.one<-arca%>%
    dplyr::select(-T_group, -group)%>%
    left_join(arca.iter%>%
                dplyr::select(calc,sample=i),
              by=c('sample'))%>%
    rbind(arca.c)%>%
    group_by(species, calc)%>%
    summarize(mean =mean(herb_mg_dens, na.rm=TRUE), sd=sd(herb_mg_dens, na.rm=TRUE), se=se(herb_mg_dens), n=length(unique(sample)))%>%
    melt(id.vars=c('species','calc'))%>%
    dcast(species~calc+variable)%>%
    mutate(iter = paste(colnames(arca.iter.filt)[i]))#result in group means for single iteration
  ##combine with all other iterations
  iter.means<-rbind(iter.means, arca.one)
}
str(iter.means)

##bird effect - LRR - ln(control/treatment)
arca.lrrs.all<-summary(escalc('ROM', m1i=ID_C_mean, m2i=ID_T_mean, sd1i=ID_C_sd, sd2i=ID_T_sd, n1i=ID_C_n, n2i=ID_T_n, data=iter.means))
View(arca.lrrs.all) ##all possible LRR & DD combos from given data for ARCA

##add var for which iteration it was
##all possible LRR for ARCA
ggplot(arca.lrrs.all, aes(reorder(iter, yi), yi))+
  geom_point()+
  geom_errorbar(aes(ymin=yi-sei, ymax=yi+sei), width=0)+
  geom_hline(yintercept=0, lty='dashed')+
  ylim(-1.5,1.5)+
  labs(x='Combination', y='Bird effect - ln(birds/no birds)')

##LRR and DD
ggplot(arca.lrrs.all, aes(DD_mean, yi))+
  geom_point()+
  geom_errorbar(aes(ymin=yi-sei, ymax=yi+sei), width=0)+
  geom_errorbarh(aes(xmin=DD_mean-DD_se, xmax=DD_mean+DD_se))+
  labs(x='Herbivore density in bird exclusion', y='Bird effect - ln(birds/no birds)')
#there is a linear relationship between DD and LRR

summary(lm(yi~DD_mean, data=arca.lrrs.all))
##can use the coef? sig?

##for all of the possible combinations - only use those in which the group means within the T group do not rejet the null hypothesis that the means are differet

########################
## run same procedure for all species
all.plant<-read.csv('/Users/colleennell/Dropbox/Projects/CSS exclusion/data/2018/CSS_plant_groups.csv')
t.plant<-all.plant%>%filter(treat=='T')%>%mutate(group = ifelse(T_group == 'ID_T', 'T_ID', ifelse(T_group == 'DD', 'T_DD', NA)))%>%transform(group=as.factor(group))
c.plant<-all.plant%>%filter(treat=='C')%>%dplyr::select(-T_group)%>%mutate(calc = 'ID_C')

##combinations
out.df<-NULL
for (i in levels(t.plant$species)){
  sp.df<-t.plant%>%filter(species == i)
  cob<-t(combn(as.character(sp.df$sample), 4, simplify=TRUE))
  sp.comb<-as.data.frame(as.table(cob))%>%dcast(Var1~Var2)%>%mutate(species = paste(i), calc='DD')
  out.df<-rbind(out.df, sp.comb)
}
str(out.df)#490 total combinations across all species

#organize
sp.out<-out.df
sp.melt<-sp.out%>%melt(id.vars=c('species','Var1','calc'))
str(sp.melt)
sp.cast<-sp.melt%>%dcast(species+variable+calc~Var1)

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

##now cast so each colum is a combination with the treatments assigned
sp.iter<-sp.melt%>%dcast(species+variable+calc~Var1)
View(sp.iter)

sp.ind.iter<-split(sp.melt, sp.melt$species, drop=TRUE)##individual dfs for each species
str(sp.ind.iter)
i.melt<-sp.ind.iter[['SAME']]#%>%dcast(species+variable+calc~Var1)
str(i.melt)
i.iter<-i.melt%>%dcast(species+variable+calc~Var1)
str(i.iter)
str(i.melt[[1]])
str(arca.melt)
unique(sp.melt$species)

#is the difference between the two group means more than expected by random?
##could take a while
#start at: 8:34
i.out<-NULL
for (i in unique(sp.melt$species)){
  i.melt<-sp.ind.iter[[i]]
  i.df<-t.plant%>%filter(species == i)
  i.iter<-i.melt%>%dcast(species+variable+calc~Var1)
  
  for (j in 4:length(i.iter)){
    i.rand<-i.df%>%left_join(i.iter%>%dplyr::select(calc,sample=j), by='sample')
    obs <- diff(mean(herb_mg_dens~calc, data=i.rand))
    i.h0<-do(1000)*diff(mean(herb_mg_dens~shuffle(calc), data=i.rand))
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
str(i.out)
write.csv(i.out, '/Users/colleennell/Dropbox/Projects/CSS exclusion/data/2018/CSS_all_possible_lrr.csv',row.names=FALSE)
#p values for all possible groupings
##filter to only those with non significant differences
View(i.out)
i.out.1<-i.out%>%filter(p>=0.1)
get.rid<-i.out%>%filter(p<=0.1)
anti.rows<-get.rid%>%dplyr::select(species, iter)
##41 to drop
dim(i.out.1)
dim(i.out)
str(sp.melt)
str(anti.rows)

sp.melt.filt<-sp.melt%>%anti_join(anti.rows, by=c('species','Var1'='iter'))
str(sp.melt.filt)
write.csv(sp.melt.filt, '/Users/colleennell/Dropbox/Projects/CSS exclusion/data/2018/CSS_filtered_combos_1.csv',row.names=FALSE)

########################
##now using combinations that do not reject null
#iterate through each column to calculate the LRR and mean DD +sd, se, n
sp.means.all<-NULL
for (i in unique(sp.melt.filt$species)){
  i.melt.filt<-sp.melt.filt%>%filter(species == i)
  i.iter.filt<-i.melt.filt%>%dcast(species+variable+calc~Var1)
  i.df<-t.plant%>%filter(species == i)
  i.c<-c.plant%>%filter(species == i)
  
  for (j in 4:length(colnames(i.iter.filt))){
    i.one<-i.df%>%
      dplyr::select(-T_group, -group)%>%
      left_join(i.iter.filt%>%
                  dplyr::select(calc,sample=j),
                by=c('sample'))%>%
      rbind(i.c)%>%
      group_by(species, calc)%>%
      summarize(mean =mean(herb_mg_dens, na.rm=TRUE), sd=sd(herb_mg_dens, na.rm=TRUE), se=se(herb_mg_dens), n=length(unique(sample)))%>%
      melt(id.vars=c('species','calc'))%>%
      dcast(species~calc+variable)%>%
      mutate(iter = paste(colnames(i.iter.filt)[j]))#result in group means for single iteration
    ##combine with all other iterations
    sp.means.all<-rbind(sp.means.all, i.one)
    
  }
}

str(sp.means.all)#all possible group means for each combination for each species
write.csv(sp.means.all, '/Users/colleennell/Dropbox/Projects/CSS exclusion/data/2018/CSS_all_sp_means.csv',row.names=FALSE)

##bird effect - LRR - ln(control/treatment)
lrrs.all<-summary(escalc('ROM', m1i=ID_C_mean, m2i=ID_T_mean, sd1i=ID_C_sd, sd2i=ID_T_sd, n1i=ID_C_n, n2i=ID_T_n, data=sp.means.all))
View(lrrs.all) ##all possible LRR & DD combos from given data for ARCA
write.csv(lrrs.all, '/Users/colleennell/Dropbox/Projects/CSS exclusion/data/2018/CSS_all_lrrs.csv',row.names=FALSE)

ggplot(lrrs.all, aes(reorder(species, yi), yi))+
  geom_errorbar(aes(ymin=yi-sei, ymax=yi+sei), width=0, color='grey')+
  geom_hline(yintercept=0, lty='dashed')+
  geom_point(aes(color=species), alpha=.7)+
  ylim(-1.5,1.5)+
  labs(x='Species', y='Bird effect - ln(birds/no birds)')

##LRR and DD
ggplot(lrrs.all, aes(DD_mean, yi))+
  geom_errorbar(aes(ymin=yi-sei, ymax=yi+sei), width=0, color='grey')+
  geom_errorbarh(aes(xmin=DD_mean-DD_se, xmax=DD_mean+DD_se), color='grey')+
  geom_point(aes(color=species))+
  geom_hline(yintercept=0, lty='dashed')+
  labs(x='Herbivore density in bird exclusion', y='Bird effect - ln(birds/no birds)')
#there is a linear relationship between DD and LRR

#### next - permutate regression
# for each of the 9 species - all possible combos
##combinations - 9 with ~40 within each
choose()
##choose 9 from 9 groups of ~40
#what are the possible combinations of those 9?

##restricted permutation without replacement
library(permute)
set.seed(4)


lrr.avg<-lrrs.all%>%group_by(species)%>%
  summarize_at(vars(yi, DD_mean, ID_C_mean, ID_T_mean), funs(mean, se, n=length))##confint??
View(lrr.avg)

fit<-lm(yi~species, data=lrrs.all)
ci95
ci95<-as.data.frame(confint(fit))
ci95$species<-rownames(confint(fit))
ci95$species<-gsub('species','',ci95$species)
#ci95$species<-gsub("(Intercept)", 'ARCA', ci95$species, fixed=TRUE)
ci95

lrr.95ci<-lrr.avg%>%left_join(ci95, by='species')%>%
  mutate(species = ifelse(species == '(Intercept)', 'ARCA', species))%>%
  mutate(ci.ub =ifelse(species != 'ARCA',`97.5 %` -0.19848134, -0.19848134),
         ci.lb=ifelse(species != 'ARCA',`2.5 %`-0.37975306, -0.37975306))
lrr.95ci

#fit dd
fit.dd<-lm(DD_mean~species, data=lrrs.all)
ci95.dd<-as.data.frame(confint(fit.dd))
ci95.dd$species<-rownames(confint(fit.dd))
ci95.dd$species<-gsub('species','',ci95.dd$species)
str(ci95.dd)
colnames(ci95.dd)

#add direct defense w/95ci also
dd.cis<-ci95.dd%>%dplyr::select(species,dd.lb =`2.5 %`, dd.ub=`97.5 %`)%>%
  mutate(species = ifelse(species == '(Intercept)', 'ARCA', species))%>%
  mutate(dd.ci.ub =ifelse(species != 'ARCA',dd.ub+0.020831628, 0.020831628),
         dd.ci.lb=ifelse(species != 'ARCA', dd.lb+0.0171131696,0.0171131696))
dd.95ci<-lrr.95ci%>%left_join(dd.cis,by='species')

##mean and 95ci predicted from model
ggplot(lrr.95ci, aes(reorder(species, yi_mean), yi_mean))+
  geom_errorbar(aes(ymin=ci.lb, ymax=ci.ub), width=0, color='black')+
  geom_hline(yintercept=0, lty='dashed')+
  geom_point(aes(color=species),size=2, alpha=1)+
  labs(x='Species', y='Bird effect - ln(birds/no birds)')

##LRR and DD
str(dd.95ci)
ggplot(dd.95ci, aes(DD_mean_mean, yi_mean))+
  geom_errorbarh(aes(xmin=dd.ci.lb, xmax=dd.ci.ub), color='grey')+
  geom_errorbar(aes(ymin=ci.ub, ymax=ci.ub), size-4,color='black')+
  geom_point(aes(color=species))+
  geom_hline(yintercept=0, lty='dashed')+
  labs(x='Herbivore density in bird exclusion', y='Bird effect - ln(birds/no birds)')+
  geom_smooth(method='lm',se=FALSE, lty='dashed',color='grey')

summary(lm(yi_mean~log(1/DD_mean_mean), data=dd.95ci))
##p=0.23

write.csv(dd.95ci, '/Users/colleennell/Dropbox/Projects/CSS exclusion/data/2018/CSS_lrr_modci.csv', row.names=FALSE)
