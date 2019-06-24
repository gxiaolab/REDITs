#Holds functions to run REDIT-Regression

#function to validate user input data
validate_input_data = function(data,covariates,covariates_to_get_p_values){
	#make sure no data is NA
	if(any(is.na(data))){ #try with data= matrix(c(NA,2,3,4),nrow=2)
		stop("data cannot have NA values. You can probably substitute zero for missing (NA) values")
	}
	if(any(is.na(covariates))){ #try with covariates = data.frame(hello=c(NA,2,3))
		stop("covariates matrix cannot have NA values")
	}
	#make sure data is a matrix
	if(!is.matrix(data)){
		stop("data argument must be a matrix") #try on data = data.frame(apple=c(1,2),pears=c(2,3))
	}
	#make sure data is a 2xn matrix 
	if(nrow(data) != 2){
		stop("data argument must be a 2xn matrix") #try on data = matrix(c(1,2,3,2,3,4),nrow=3)
	}
	#make sure data is numeric, without decimals
	if(!all(is.numeric(data))){ 
		stop("elements in data must be integers") #try on data = matrix(c('1',2,3,4,5,6),nrow=2)
	}
	if(!all(data %% 1 == 0)){
		stop("elements in data must be integers") #try on data = matrix(c(0.2,2,1,7),nrow=2)
	}
	#make sure covariates is a data.frame
	if(!is.data.frame(covariates)){
		stop("covariates must be a data.frame") #test on covariates = matrix(c(1,2,3,4))
	}
	#make sure number of covariate entries and columns in data match 
	if(ncol(data) != nrow(covariates)){
		stop("nrow(covariates) must equal ncol(data)") 
	}
	#make sure each column in covariates is either integer, numeric, factor, or character
	if(any(! sapply(covariates,class) %in% c('integer','numeric','factor','character'))){
		stop("each covariates in covariates data.frame must be of type integer, numeric, factor, or character")	
	}
	if(any(! covariates_to_get_p_values %in% names(covariates))){
		stop("values in covariates_to_get_p_values must correspond to colnames in covariates")
	}
}

#function to refactor categorical variables into dummy variables. Returns a matrix
convert_categorical_variables_into_dummy = function(covariates){
	covariate_types = sapply(covariates,class)
	non_categorical_covariates = covariates[,covariate_types %in% c('integer','numeric'),drop=FALSE]
	categorical_covariates = covariates[,covariate_types %in% c('character','factor'),drop=FALSE]
	if(ncol(categorical_covariates)==0){ #if have no categorical variables
		return(as.matrix(covariates))			
	}
	names(categorical_covariates) = paste0(names(categorical_covariates),'.')
	dummy_matrix = model.matrix(~.,data=categorical_covariates)
	dummy_matrix = as.data.frame(dummy_matrix)
	dummy_matrix[['(Intercept)']] = NULL
	new_covariates = as.matrix(cbind(non_categorical_covariates,dummy_matrix))
	return(new_covariates)
}

#fit beta binomial regression with identity link (no logistic fuction) using maximum likelihood estimation
get_maximum_likelihood_beta_binomial_regression_linear_mapping = function(data,covariates=data.frame()){ 
	#data is 2xn matrix of n samples. first row are number of G reads. Second row are number of A reads
	#covariates is n x k   data.frame of k covariates for n samples. The sample order must correspond with the data argument. The colnames will be the names of the covariates that this method will test are all significant
	covariate_names = colnames(covariates)
	all_covariate_names = c(covariate_names,'intercept','sigma')
	output_object = list() #will have p values for each covariate
	G_reads = data[1,]
	A_reads = data[2,]
	N = G_reads+A_reads
	y = G_reads
	linear_covariates_comma = paste(covariate_names,collapse=',')
	if(length(covariate_names) > 0){
		mu_regression_term = paste0("as.numeric( covariates"," %*% ",'matrix(c(',linear_covariates_comma,'),ncol=1)','+intercept)')
		#also make an expression for abstracting parameter names to and from the optim function
		parameter_terms = rep(NA,length(all_covariate_names))
		for(i in 1:length(all_covariate_names)){
			the_covariate = all_covariate_names[i]
			parameter_terms[i] = paste0(the_covariate,'=par[',i,']')
		}
		parameter_term = paste(parameter_terms,collapse=';')
	}else{
		mu_regression_term = paste0('intercept')
	}
	if(length(covariate_names) > 0){ #for maximum likelihood convert back to alpha, beta convention and set LL to 0 if either are <= 0
		likelihood_function = paste0('LL = function(par){',
			parameter_term,'
			mu = ',mu_regression_term,'
			if(any(mu > 1 | mu < 0)){
				return(NA)
			}
			if(sigma <= 0){ #sigma must > 0
				return(NA)
			}
			mu = ifelse(mu==1, 0.9999999999, ifelse(mu==0, 0.00000000001, mu))
			alpha = mu / sigma
			beta = (1-mu) / sigma
			log_likelihoods = -sum( lgamma(alpha+beta) - lgamma(alpha) - lgamma(beta) - lgamma(alpha+beta+N) + lgamma(alpha+y) + lgamma(N-y+beta) )
			return(log_likelihoods)
		}')
	}else{
		likelihood_function = paste0('LL = function(par){
			intercept=par[1]; sigma=par[2]
			mu = intercept
			#forbid parameters that make mu estimates illegal (outside 0 or 1)
			if(any(mu > 1 | mu < 0)){
				return(NA)
			}
			if(sigma <= 0){ #sigma must > 0
				return(NA)
			}
			mu = ifelse(mu==1, 0.9999999999, ifelse(mu==0, 0.00000000001, mu))
			alpha = mu / sigma
			beta = (1-mu) / sigma
			log_likelihoods = -sum( lgamma(alpha+beta) - lgamma(alpha) - lgamma(beta) - lgamma(alpha+beta+N) + lgamma(alpha+y) + lgamma(N-y+beta) )
			return(log_likelihoods)
		}')
	}		
	eval(parse(text=likelihood_function))
	if(length(covariate_names) > 0){ 
		covariates_initial_points = rep(NA,length(all_covariate_names))
		for(i in 1:length(all_covariate_names)){
			the_covariate = all_covariate_names[i]
			if(the_covariate != 'sigma'){
				covariates_initial_points[i] = 0
			}else{
				covariates_initial_points[i] = 0.5
			}
		}
	}else{
		covariates_initial_points = c(0,0.5) #c(intercept, sigma)
	}
	parameters = optim(par=covariates_initial_points,LL)
	log_likelihood = -1 * parameters$value
	all_parameters = parameters$par; names(all_parameters) = all_covariate_names
	converged = parameters$convergence
	output_object = list()
	output_object[['log_likelihood']] = log_likelihood
	output_object[['parameters_estimates']]=all_parameters
	output_object[['converged']]=converged #0 if converged. 1 otherwise
	return(output_object)
}


#main function to run beta-binomial regression
REDIT_regression = function(data,covariates=data.frame(),covariates_to_get_p_values=c()){
	#data is 2xn matrix of n samples. first row are number of G reads. Second row are number of A reads
	#coviarates is n x k   matrix of k coviarates for n samples. The sample order must correspond with the data argument. The colnames will be the names of the covariates that this method will test are all significant. If user omits this argument, then assumed to run beta binomial regression without any covariates
	#covariates_to_get_p_values: vector of names of the covariates from which you want to get p values for. If omitted then algorith will find p values for all covariates in the covariate matrix
	validate_input_data(data,covariates,covariates_to_get_p_values)
	covariates = convert_categorical_variables_into_dummy(covariates)
	full_model = get_maximum_likelihood_beta_binomial_regression_linear_mapping(data,covariates)
	covariate_names = colnames(covariates)
	if(length(covariates_to_get_p_values) == 0){ #then get all covariates
		covariates_to_get_p_values = covariate_names
	}
	#determine the p value of each covariate desired
	output_object = list()
	output_object[['parameter_estimates']] = full_model$parameters_estimates
	output_object[['converged']] = full_model$converged
	for(covariate in covariates_to_get_p_values){
		mini_covariates = covariates[,setdiff(covariate_names,covariate),drop=FALSE]
		null_model = get_maximum_likelihood_beta_binomial_regression_linear_mapping(data,mini_covariates)
		D0 <- null_model$log_likelihood
		D1 <- full_model$log_likelihood
		df = 1
		chi = -2 * (D0 - D1)
		p.val = pchisq(chi,df,lower.tail=FALSE)
		output_object[[paste0(covariate,'.p.value')]] <- p.val
	}
	return(output_object)
}

#Developmental variables
#covariates = data.frame( age=c(1,1,1,8,8,8),weight=c(2,9,2,8,3,5),sex=c("unknown",'M','F','M','F','M'), dynamite=c('field','mine','sourcery','field','mine','sourcery' ),useless=c('a','a','a','a','a','a') )
#covariates = data.frame( age=c(1,1,1,8,8,8),weight=c(2,9,2,8,3,5),sex=c("unknown",'M','F','M','F','M'))
#covariates = data.frame( sex=c("unknown",'M','F','M','F','M'))
#data = matrix(c(0,1,2,10,10,10,10,9,8,0,0,0),nrow=2,byrow=TRUE)
