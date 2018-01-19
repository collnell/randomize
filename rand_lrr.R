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

#####

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

sum(arca.h0$ID_T>=obs)/length(arca.h0$ID_T)##P = 0.951...so this grouping combination does not reject null - should be more reliable

##do this FOR ALL COMBOS - can use same null
#the difference between the two group means is more than expected by random
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

######

#control treatment plants to combine with other df
arca.c<-all.plant%>%filter(treat=='C' & species =='ARCA')%>%dplyr::select(-T_group)%>%mutate(calc = 'ID_C')
#use by= to iterate through columns, and combine with the cacl col?

##now using combinations that do not reject null
#iterate through each column to calculate the LRR and mean DD +sd, se, n

##now iterate that through all possible treatment combinations
iter.means<-NULL
for (i in 3:length(colnames(arca.iter))){
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
    mutate(iter = paste(colnames(arca.iter)[i]))#result in group means for single iteration
  ##combine with all other iterations
  iter.means<-rbind(iter.means, arca.one)
  
}
str(iter.means)
##now take LRR
arca.lrrs.all<-summary(escalc('ROM', m1i=ID_C_mean, m2i=ID_T_mean, sd1i=ID_C_sd, sd2i=ID_T_sd, n1i=ID_C_n, n2i=ID_T_n, data=iter.means))
View(arca.lrrs.all) ##all possible LRR & DD combos from given data for ARCA

##add var for which iteration it was
##all possible LRR for ARCA
ggplot(arca.lrrs.all, aes(iter, yi))+geom_point()+geom_errorbar(aes(ymin=yi-sei, ymax=yi+sei), width=0)+geom_hline(yintercept=0, lty='dashed')+ylim(-1.5,1.5)

##LRR and DD
ggplot(arca.lrrs.all, aes(DD_mean, yi))+geom_point()+geom_errorbar(aes(ymin=yi-sei, ymax=yi+sei), width=0)+
  geom_errorbarh(aes(xmin=DD_mean-DD_se, xmax=DD_mean+DD_se))
#there is a linear relationship between DD and LRR
summary(lm(yi~DD_mean, data=arca.lrrs.all))
##can use the coef? sig?

##for all of the possible combinations - only use those in which the group means within the T group do not rejet the null hypothesis that the means are differet

