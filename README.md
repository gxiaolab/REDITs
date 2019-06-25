# Statistical inference of differential RNA editing sites from RNA-sequencing data 

## Description
This repository holds two functions in R to run REDITs (RNA editing tests) for calling differential RNA editing sites:
1. RNA editing sites that are significantly different between case-control (or condition1 vs condition2) cohorts
   - Handled by the REDIT-LLR (RNA editing test - log likelihood ratio)
2. RNA editing sites that significantly correlate with categorical (e.g. sex, gender) and/or quantitative (expression of ADAR, age) variables
   - Handled by REDIT-Regression (RNA editing test - regression)

These tests consider both variance in editing from biological variance and intrinsic inaccuracy from calculating editing from count data such as RNA-seq. Thus they wield greater power and lower false positives, at and below the 5% false positive threshold, than commonly used alterantives such as the t-test, Wilcoxon's rank-sum test, or pooled Fisher's Exact test.

## Related publication
This work is from: <link to publication to be added>

## Table of Contents
1. [Installation](https://github.com/tranStephen/beta_binomial_tests/blob/master/README.md#installation)
2. [Usage](https://github.com/tranStephen/beta_binomial_tests/blob/master/README.md#usage)
   - [Beta binomial log likelihood](https://github.com/tranStephen/beta_binomial_tests/blob/master/README.md#beta-binomial-log-likelihood-test)
   - [Beta binomial regression](https://github.com/tranStephen/beta_binomial_tests/blob/master/README.md#beta-binomial-regression)
   - [Parallelization](https://github.com/tranStephen/beta_binomial_tests/blob/master/README.md#example-of-running-beta-binomial-tests-with-parallelization)
3. [Credits](https://github.com/tranStephen/beta_binomial_tests/blob/master/README.md#credits)
4. [License](https://github.com/tranStephen/beta_binomial_tests/blob/master/README.md#license)

## Installation
**You only need base** [R](https://www.r-project.org)

Both beta-binomial tests use the optim function from the stats package which is a base package that is automatically included and loaded in R.

If you want to use the below example code with parallelization, then you will need to install 
1. doParallel
2. foreach 

But the actual beta-binomial tests, which alone are quite fast, do not use these packages. Nor do they need parallelization, unless you plan on running on a dataset with millions of editing sites.

## Usage
## REDIT-LLR
```
source("beta_binomial_log_likelihood_test.R")
```
### Function documentation
**REDIT_LLR**(*data*, *groups*)

*data*: a 2xn numeric matrix. The first row holds the number of counts (e.g. RNA-seq reads) supporting editing. The second row holds the number of counts supporting non-editing. Each column corresponds to data from one sample.

*groups*: a character vector corresponding to the condition/cohort/diease-control membership of corresponding column in the 2xn data matrix. It must have exactly 2 unique strings. `length(unique(groups)) == 2`

*Returns*: A list
- *data and groups*: The original data and groups parameters
- *mle.for.group.disease, mle.for.group.control, and mle.for.null.model*: the maximum likelihood estimates for the underlying beta distribution of the disease and controls groups and null model respectively. You can obtain maximum likelihood estimates of the average editing level per condition by dividing their respective alpha / (alpha+beta). Note, however, that the estimates are NOT stable enough to obtain a maximum likelihood variance of editing level.
- *log.likelihood.for.group.disease, log.likelihood.for.group.control, and log.likelihood.for.null*: the log likelihoods for the corresponding maximum likelihood estimates
- *p.value*: the p-value testing whether this editing site is statistically significant between disease-control groups.

The actual names of the list elements will vary based on the two unique strings in the *groups* parameter.

### Code example 
```
>source("beta_binomial_log_likelihood_test.R")
>the_data = matrix( c(1,9, 2,9, 8,1,10,0),nrow=2)
>the_groups = c('disease','disease','control','control');
>beta_binomial_log_likelihood_test(data=the_data, groups=the_groups)
>> 
$data
     [,1] [,2] [,3] [,4]
[1,]    1    2    8   10
[2,]    9    9    1    0

$groups
[1] "disease" "disease" "control" "control"

$mle.for.group.disease
   alpha     beta 
165118.1 990660.1 

$mle.for.group.control
    alpha      beta 
977596.39  54293.39 

$mle.for.null.model
   alpha     beta 
1.025074 1.000000 

$log.likelihood.for.group.disease
[1] -8.61245

$log.likelihood.for.group.control
[1] -3.917653

$log.likelihood.for.null
[1] -18.0895

$p.value
[1] 0.003851092
```
## REDIT-Regression
```
source("beta_binomial_regression.R")
```
### Function documentation
**REDIT_regression**(*data*, *covariates*, *covariates_to_get_p_values*)

*data*: a 2xn numeric matrix. The first row holds the number of counts (e.g. RNA-seq reads) supporting editing. The second row holds the number of counts supporting non-editing. Each column corresponds to data from one sample.

*covariates*: a dataframe where each column holds a categorical (chacter/factor) variable or quantitative (numeric, integer) variable. Each row corresponds 1 to 1 with each respective column in the 2xn data matrix.

*covariates_to_get_p_values*: a character vector specifying for which variables to get p-values. The strings must exactly match to colnames from the covariates dataframe. If you omit this argument, then the function will calculate p values for all variables. Note that all variables in the covariates dataframe will be included in the regression regardless if the covariates_to_get_p_values argument is used or not.

*Returns* a list
- *parameter_estimates*: the maximum likelihood estimates of all regression terms
- *converged*: 0 if maximum likelihood converged. 1 if not. We found that statistical inference is generally robust to whether or not maximum likelihood converged.
- *age.p.value*: p-value of the effect of age on editing
- *sex.M.p.value and sex.unknown.p.value*: p-value of the effect of sex on editing. Note that the beta_binomial_regression function automatically converts categorical variables into dummy representations.

The actual names of the list elements will vary based on the colnames in the *covariates* parameter.


### Code example
```
>source("beta_binomial_regression.R")
>the_data = matrix( c(1,9, 2,9, 8,1,10,0),nrow=2)
>the_groups = c('disease','disease','control','control');
>beta_binomial_log_likelihood_test(data=the_data, groups=the_groups)
> the_data = matrix(c(0,1,2,10,10,10,10,9,8,0,0,0),nrow=2,byrow=TRUE) 
>the_covariates = data.frame( age=c(1,1,1,8,8,8), sex=c("unknown",'M','F','M','F','M'))
>beta_binomial_regression(data=the_data, covariates=the_covariates, covariates_to_get_p_values=c('age','sex')) #you can omit the covariates_to_get_p_values argument
>>
$parameter_estimates
          age         sex.M   sex.unknown     intercept         sigma 
 1.183318e-01  4.154607e-07 -6.952360e-03  5.334515e-02  4.342248e-01 

$converged
[1] 0

$age.p.value
[1] 0.007673521

$sex.M.p.value
[1] 1

$sex.unknown.p.value
[1] 0.6962553
```

## Example of running beta-binomial tests with parallelization
Parallelization doesn't actually use any code implemented by us. You are free to use any packages suitable for parallelization. We provide below an example using doParallel and foreach

You need to install 
1. doParallel 
2. foreach 

to run the below example code
### Example code
```
>library(doParallel)
>library(foreach)
>source("beta_binomial_regression.R")

#initiate a cluster
>noCores = detectCores() -1
>cl = makeCluster(noCores,outfile="")
>registerDoParallel(cl,cores=noCores)

#create some data for illustration
>number_of_editing_sites_to_test = 2
>G_reads = matrix( c(0,    1,    2,   10,   10,   10,
  2, 5, 2, 3, 8, 1), nrow=2, byrow=TRUE)
>A_reads = matrix( c(10,    9,    8 ,   0,    0,    0,
  8, 5, 8, 7, 1, 10), nrow=2, byrow=TRUE)
>the_covariates = data.frame( age=c(1,1,1,8,8,8), sex=c("unknown",'M','F','M','F','M'))

#run the parallelization
>output_matrix = foreach(i=1:number_of_editing_sites_to_test,.combine='rbind') %dopar%{
  G_reads_of_editing_site = G_reads[i,]
  A_reads_of_editing_site = A_reads[i,]
  bb_regression_info = beta_binomial_regression(data=rbind(G_reads_of_editing_site,A_reads_of_editing_site), covariates=the_covariates)
  return( as.matrix( data.frame(p_value_age= bb_regression_info$age.p.value )) )
}
>output_matrix
>>
     p_value_age
[1,] 0.007673521
[2,] 0.708027338

>stopCluster(cl) #to stop using the cluster you created
```
## Credits
[Grace Xiao webpage](https://www.ibp.ucla.edu/research/xiao/Home.html)

[Grace Xiao github](https://github.com/gxiaolab)

[Qing Zhou webpage](http://www.stat.ucla.edu/~zhou/)

## License
The code here is freely distributed with no restrictions.


