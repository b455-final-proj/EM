list.of.packages = c('knitr', 'mvtnorm')
new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
rm(list.of.packages, new.packages)

# Generates MOG data given # of observations, means and standard deviatins of each Normal Distribution

# n_data: integer number of data points to generate
# means: vector of means for each class
# sds: vector of standard deviations of each class
# priors: vector of prior probabilities for each class

# e.g.: MOG for 100 observations from N(0,1) and 500 observations from N(10,2) with equal priors would be:
# generate_mog(c(100,500),c(0,10),c(1,2))

generate_mog <- function(n_data,means,sds,priors=NULL) {
  n_classes = length(means)
  if(is.null(priors)){ priors = rep(1/n_classes, n_classes) }
  classes = sample(1:n_classes, size=n_data, replace=TRUE, prob=priors)
  data = rnorm(n_data, mean = means[classes], sd = sds[classes])
  return (data)
}

EM_setup = function(n_classes, prior_choice='equal_priors_unif_params', data){
	n_vars = ncol(data)
	if(prior_choice=='equal_priors_random_params'){
		prior = matrix(rep(1/n_classes, n_classes), nrow=n_classes) #All classes are equally as likely.
		mu = matrix(runif(n_classes*n_vars, 0, 1), nrow=n_classes) #All means are initially ~Unif[0,1].
		sigma = replicate(n_classes, matrix(runif(n_vars^2, 0, 1), nrow=n_vars))
	} else {
		stop('invalid prior choice.')
	}
	return(list(prior, mu, sigma))
}

expectation = function(data, prior, mu, sigma){
	#Get the probabilities of each data point belonging to each class.
	n_classes = ncol(prior)
	n_vars = ncol(data)
	n_data = nrow(data)
	P.Xj.and.Ci = matrix(rep(NA, length=n_classes*n_data), ncol=n_classes)
	for(i in 1:n_classes){
		P.Ci.and.Xj[i] = prior[i] * dmvnorm(data, mu=mu[i,], sigma=sigma[,,i])
	}
	P.Ci.given.Xj = P.Xj.and.Ci / apply(P.Xj.and.Ci, 1, sum)
	return(P.Ci.given.Xj)
}

maximization = function(P.Ci.given.Xj, data, prior){
  n_classes = ncol(prior)
  n_vars = ncol(data)
  n_data = nrow(data)
	N = matrix(apply(P.Ci.given.Xj, 1, sum), ncol=1)
	mu = matrix(rep(NA, length=n_classes*n_vars, 0, 1), nrow=n_classes)
	sigma = replicate(n_classes, matrix(rep(NA, n_vars^2), nrow=n_vars))
	for(i in 1:n_classes){
		mu[i,] = colMeans(data * P.Ci.given.Xj[i] / N[i])
		diff = data - mu[rep(i, nrow(data)),]
		sigma[,,i] = diff %&% t(diff) / n_data
	}
	prior = N / sum(N)
}


#prior_and_params = EM_setup(3)
#full_vector = expec(dat,prior,params)
#[params,prior] = maxim(full_vector,dat,prior)

#for i = 1:15
#    full_vector = expec(dat,prior,params)
#    [params,prior] = maxim(full_vector,dat,prior)
#end
