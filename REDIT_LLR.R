#Holds functions to run REDIT-LLR

#Fits data to a beta-binomial distribution restricted to uniform and unimodal shapes. 
get_maximm_likelihood_parameters_beta_binomial_unimodal = function(data){
	#data is matrix. First row are G reads. Second row are A reads
	N = data[1,]+data[2,]
	y = data[1,]
	LL = function(par){
		alpha = par[1]
		beta = par[2]
		if(alpha < 1 | beta < 1 | alpha > 1e6 | beta > 1e6){ #boundaries of the optimization
			return(NA)
		}
		log_likelihood = -sum( lgamma(alpha+beta) - lgamma(alpha) - lgamma(beta) - lgamma(alpha+beta+N) + lgamma(alpha+y) + lgamma(beta+N-y) )
		return(log_likelihood)
	}
	parameters = optim(par=c(1,1), LL)
	log_likelihood = -1 * parameters$value
	all_parameters = parameters$par; names(all_parameters) = c('alpha','beta')
	converged = parameters$convergence
	output_object = list()
	output_object$output_parameters = all_parameters
	output_object$converged = converged
	output_object$log_likelihood = log_likelihood
	return(output_object)
}

#a generic function for log likelihood ratio test
run_likelihood_ratio_test = function(L0,La1, La2){
	#L0 is maximum log likelihood of null hypothesis
	#La1 is maximum log likelihood of alternative hypothesis for group1
	#La2 is max log likelihood of alternative hypothesis for group2
	chi_square_stat = -2*(L0 - (La1 + La2))
	df = 2
	p_value = pchisq(chi_square_stat,df=df,lower.tail=FALSE)
	return(p_value)
}

validate_input_data = function(data,groups){ #validates user input data
	#make sure no data is NA
	if(any(is.na(data))){
		stop("data cannot have NA values. You can probably substitute zero for missing (NA) values")
	}
	if(any(is.na(groups))){
		stop("groups cannot have NA values")
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
	#make sure groups is a vector
	if(!is.vector(groups)){
		stop("groups must be a vector") #test on groups = matrix(c(1,2,3,4))
	}
	if(!is.character(groups)){
		stop("groups must be a chracter vector") #test on groups = as.factor(c('a','b','a','b')). groups=c(1,2,3,4)
	}
	#make sure number of groups and columns in data match 
	if(ncol(data) != length(groups)){
		stop("length(groups) must equal ncol(data)") 
	}
	#make sure there are only two groups
	if(length(unique(groups)) != 2){
		stop("must have exactly two groups") #test on groups=c('a','b','c','a'); data= matrix(c(1,2,3,4,5,6,7,8),nrow=2)
	}
}

#splits data by experimental condition
split_data_by_group = function(data,groups){
	type_groups = unique(groups)
	data1 = data[,groups==type_groups[1],drop=FALSE]
	data2 = data[,groups==type_groups[2],drop=FALSE]
	group1 = type_groups[1]
	group2 = type_groups[2]
	return(list(data1=data1,data2=data2,group1 = group1,group2=group2))
}

#the actual function for running the beta binomial log likelihood test
beta_binomial_log_likelihood_test = function(data,groups){
	validate_input_data(data,groups)
	split_data = split_data_by_group(data,groups)
	data1 = split_data$data1                                                                                               
	data2 = split_data$data2
	group1 = split_data$group1
	group2 = split_data$group2
	model_fit0 = get_maximm_likelihood_parameters_beta_binomial_unimodal(data=data)
	model_fit1 = get_maximm_likelihood_parameters_beta_binomial_unimodal(data=data1)
	model_fit2 = get_maximm_likelihood_parameters_beta_binomial_unimodal(data=data2)
	parameters0 = model_fit0$output_parameters
	parameters1 = model_fit1$output_parameters
	parameters2 = model_fit2$output_parameters
	L0 = model_fit0$log_likelihood
	La1 = model_fit1$log_likelihood
	La2 = model_fit2$log_likelihood
	p_value = run_likelihood_ratio_test(L0=L0,La1=La1,La2=La2)
	output_object = list()
	output_object$data = data
	output_object$groups = groups
	output_object[[paste0('mle.for.group.',group1)]] = parameters1
	output_object[[paste0('mle.for.group.',group2)]] = parameters2
	output_object[[paste0('mle.for.null.model')]] = parameters0
	output_object[[paste0('log.likelihood.for.group.',group1)]] = La1
	output_object[[paste0('log.likelihood.for.group.',group2)]] = La2
	output_object[[paste0('log.likelihood.for.null')]] = L0
	output_object[['p.value']] = p_value
	return(output_object)
}


#example:
#data = matrix(c(2,2,2,2,5,0,5,0),nrow=2)
#data = matrix(c(0,5,0,5,5,0,5,0),nrow=2)
#groups = c('a','a','b','b')
#data =matrix(c(10,0,0,10),nrow=2); groups = c('a','b')

