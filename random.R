## randomiztion tests
# https://cran.r-project.org/web/packages/randomizeR/vignettes/article.pdf

#randomization - critical goals
# 1. balance knwn and unknown covariate - structural equality of treatment groups
# 2. blinding of treatment allocation to prevent bias in selection
# 3. internal validity - basis for statistical inference

# randomization procedures - Rosenberger and Lachin 2016
# 

## data in - arthropod density across 114 plants from 2 treatments (T and C)
plants<-read.csv('data/2018/CSS_plant_data.csv')

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



##reample and calculate mean x times

