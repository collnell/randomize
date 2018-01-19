## randomiztion tests
# https://cran.r-project.org/web/packages/randomizeR/vignettes/article.pdf

#randomization - critical goals
# 1. balance knwn and unknown covariate - structural equality of treatment groups
# 2. blinding of treatment allocation to prevent bias in selection
# 3. internal validity - basis for statistical inference

# randomization procedures - Rosenberger and Lachin 2016
# 

## data in - arthropod density across 114 plants from 2 treatments (T and C)
plants<-read.csv('/Users/colleennell/Dropbox/Projects/CSS exclusion/data/2018/CSS_plant_data.csv')

## need to randomly subset within treatment T to create 2 groups T_ID & T_DD
## small sample size and high variance within groups - 
## in theory the mean density within two subsetted groups is comparable
## purpose: ln(C/T_ID) ~ T_DD  - doing this to avoid spurrious correlation

## however, seems the results are highly unstable
## outcome is variable depending on randomization
## need to reduce vulnerability to chance differences between the treatment groups
## is there a way to do all possible combinations to allocate groups
## run test for each combo
## record each outcome - composite is true result
## or ways to assess the randomization procedure?

##assessing randomization procedures
# assumptions - normally distributed responses
# treatments have equal expectation and variances


# bootstrapping an dpermutation tests
#estimation and testing
# permutation test - two group comparison - permute labels for grouping variable
# then calc sample stat (difference between 2 groups using new assignments) to empirically construct the null distribution

install.packages('mosaic')
library(mosaic)
require(mosaicData)
options(digits=3)

# selling prices for used cars, age, miles for 25 cars
Mustangs <- read.file("http://www.mosaic-web.org/go/datasets/MustangPrice.csv")

#construct 90% CI for the mean price
histogram(~Price, data=Mustangs)
mean(~Price, data=Mustangs)

simple=c(1, 2, 3, 4, 5) #with replacement
resample(simple)
resample(simple)

resample(Mustangs)
#resample mean
mean(~Price, data=resample(Mustangs))
mean(~Price, data=resample(Mustangs))
mean(~Price, data=resample(Mustangs))

#5 times
do(5) * mean(~Price, data=resample(Mustangs))

#using parallel package
set.rseed(19)

#1000 times
Mustangs.Price.boot <- do(1000) * mean(~Price, data=resample(Mustangs))
histogram(~ mean, data=Mustangs.Price.boot, xlab="Mean Mustang Price (in thousand dollars)")
#get CI 
confint(Mustangs.Price.boot, level=0.90, method="quantile")#distribution quantiles
confint(Mustangs.Price.boot, level=0.90, method="stderr") #normal theory/standard errr

############################
##simulation tests

# simulate whether the binomial outcome randomly sampled from 428 trials is even (.5)
#what proportion of trials give a result as or more extreme than observed(240)
#
prop(~ rbinom(1000, prob=0.5, size=428) >= 240)
#unlikely if null is true (equally liekly outcomes), that single outcome would occur 240 times or more
#variation between simualtions decreases with more trials

#deterministic probability calculation 
pbinom(239, prob=0.5, size=428)
binom.test(240, 248)

#explicitly simulate binomial outcome
do(1)*rflip(428)
do(100)*rflip(428)
#what fraction of the trials is the outcome 240/428 attempts?
NFL.null <- do(1000) * rflip(428)
prop(~ heads >= 240, data=NFL.null)
histogram(~ heads, groups=(heads >= 240), data=NFL.null)#shading - 
#the observed 240/420 is unliekly under the null hypothesis (outcomes occur evenly (50%))

## permutation test of means from two groups

#assign randomly to 2 groups for treatment, test for difference in mean Words between treatments
Sleep <- read.file("http://www.mosaic-web.org/go/datasets/SleepCaffeine.csv")
mean(Words ~ Group, data=Sleep)
obs <- diff(mean(Words ~ Group, data=Sleep))
obs
bwplot(Group ~ Words, data=Sleep)
#implement nul hypothesis
#scramble group with respect to outcome, Words
diff(mean(Words ~ shuffle(Group), data=Sleep))#single trial - difference in means between groups
diff(mean(Words ~ shuffle(Group), data=Sleep))       
#get distribution under null using 1000 trials
Sleep.null <- do(1000) * diff(mean(Words ~ shuffle(Group), data=Sleep))
histogram(~ Sleep, groups=(Sleep >= obs), data=Sleep.null, width=0.4,
          xlab="Distribution of difference in means\nunder the null hypothesis")
#one-sided p is the proportion of tials which yielded as extreme or more extreme as observeddifference
#calculate p by summing number of trials with extreme
#here 30/1000 trials were extreme, so P = 0.03 - and thus unliekly that the two groups have the same mean word in respective populations
z<-Sleep.null$Sleep>=obs
table(z)['TRUE']
length(z[z==TRUE])
sum(z)

sum(Sleep.null$Sleep>=obs)/length(Sleep.null$Sleep)

mean(Words~1, data=Sleep)
mean( Words ~ Group, data=Sleep )
lm( Words ~ Group, data=Sleep )
diffmean( Words ~ Group, data=Sleep )

#test whether group assignment is representative vs null hypothesis - should have equal mean density
#scramble/shuffle ID/DD assignemnt, do many trials to calculate the mean of each group, difference
all.plant<-read.csv('/Users/colleennell/Dropbox/Projects/CSS exclusion/data/2018/CSS_plant_groups.csv')
#just bird exclusion treatment, 2 levels
t.plant<-all.plant%>%filter(treat == 'T')%>%mutate(group = ifelse(T_group == 'ID_T', 'T_ID', ifelse(T_group == 'DD', 'T_DD', NA)))%>%transform(group=as.factor(group))

mean(herb_mg_dens~group+species, data=t.plant) #mean for each species in the group assignments
obs <- diff(mean(herb_mg_dens~group|species, data=t.plant))
obs #difference between groups
.02031-.01609
0.06961-.02411

#this isnt working because it is not doing it within each species
#for each species, shuffle group assignment and claculate means and idfference between groups

#test wheterh groups assignment is representative - asess groupings
#for combinations in which means for 2 groups are within the normal distribution

#ARCA
arca<-t.plant%>%filter(species=='ARCA')
mean(herb_mg_dens~group, data=t.plant)
obs <- diff(mean(herb_mg_dens~group+species, data=arca))
obs #observed difference in means
bwplot(herb_mg_dens~group, data=arca)
#null hypothesis - groups randomly assigned
arca.null<-do(1000)*diff(mean(herb_mg_dens~shuffle(group, groups=species), data=t.plant))
arca.null
histogram(~T_ID, groups=(T_ID >= obs), data=arca.null)
#likely they do have same mean? - does not reject nul hypothesis, so grouping vars are legit?

diffmean(herb_mg_dens~species|group, data=t.plant, only.2=FALSE)
?diffmean

#randomizationw ith linear models
#more generalized - multiple explanatory variables
#permutation test on the difference between densities 
diffprop( homeless ~ sex, data=HELPrct)

Mustangs.lm.boot <- do(1000) * lm(Price ~ Miles, data=resample(Mustangs))
confint(Mustangs.lm.boot)
resample(Mustangs)

HELPrct.null <- do(1000) * lm(homeless=="homeless" ~ shuffle(sex), data=HELPrct)
prop(~ (abs(sexmale) > 0.1146), data=HELPrct.null)

Mustangs.boot3 <- do(1000) * lm( Price ~ Miles + Age, data=resample(Mustangs))
confint(Mustangs.boot3)
#shuffle Age to compare the change in R to what would be expected under the null
Mustangs.Age.boot <- do(1000) * lm( Price ~ Miles + shuffle(Age), data=Mustangs )

#notes from https://cran.r-project.org/web/packages/mosaic/vignettes/Resampling.pdf

##############################
#1. use randomization test to see whether group means differ than null hypothesis (random)
#2. for all permutations or treatment assignments that do not differ from null, calculate LRR and DD - obtain distribution of these values
#3. for each iteratio of LRR and DD, run regression and recrod P, R2, sign
#4. then the results are....?

##all possible permutations of gorup assignments within arca
arca
##possible combinations of group such that there is 4 in each group
library(gtools)
cob<-t(combn(as.character(arca$sample), 4, simplify=TRUE))
cob

arca.comb<-as.data.frame(as.table(cob))%>%dcast(Var1~Var2)%>%mutate(calc='DD')##all possible combos within a single treatment, var1 is unique combo
str(arca.comb)
View(arca.comb)
arca.comb[1,-1]

arca.melt<-arca.comb%>%melt(id.vars=c('Var1','calc'))
arca.cast<-arca.melt%>%dcast(variable+calc~Var1)
#View(arca.melt)
#View(arca.cast)

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
View(arca.melt)

##now cast arca.melt so each colum is a combination with the treatments assigned
arca.iter<-arca.melt%>%dcast(variable+calc~Var1)
View(arca.iter)

##now iterate through each column to calculate the LRR and mean DD +sd, se, n

#control treatment plants to combine with other df
arca.c<-all.plant%>%filter(treat=='C' & species =='ARCA')%>%dplyr::select(-T_group)%>%mutate(calc = 'ID_C')
#use by= to iterate through columns, and combine with the cacl col?

### a single iteration
#combine by that column, selecting only the 2 cols
arca.all<-arca%>%dplyr::select(-T_group, -group)%>%
  left_join(arca.iter[,c(2,3)], by=c('sample'='A'))%>%
  rbind(arca.c) ##should be all samples with their groups assigned
#hen use that to calculate LRR and DD
arca.mean<-arca.all%>%group_by(species, calc)%>%
  summarize(mean =mean(herb_mg_dens, na.rm=TRUE), sd=sd(herb_mg_dens, na.rm=TRUE), se=se(herb_mg_dens), n=length(unique(sample)))%>%
  melt(id.vars=c('species','calc'))%>%
  dcast(species~calc+variable)
View(arca.mean)##means for each group for single iteration

#LRR 
arca.lrr<-summary(escalc('ROM', m1i=ID_C_mean, m2i=ID_T_mean, sd1i=ID_C_sd, sd2i=ID_T_sd, n1i=ID_C_n, n2i=ID_T_n, data=arca.mean))
View(arca.lrr)

ggplot(arca.lrr, aes(1, yi))+geom_point()+geom_errorbar(aes(ymin=yi-sei, ymax=yi+sei), width=0)+geom_hline(yintercept=0, lty='dashed')+ylim(-1.5,1.5)


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

