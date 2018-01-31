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

#install.packages('mosaic')
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



summary(aov(lm(lrr~DD_mean, data=resample(lrr.df, replace=TRUE, groups=lrr.df$species))))


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

