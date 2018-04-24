list.of.packages = c('knitr', 'mvtnorm')
new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
rm(list.of.packages, new.packages)

# Generates MOG data given # of observations, means and standard deviatins of each Normal Distribution

# n.data: integer number of data points to generate
# means: vector of means for each class
# sds: vector of standard deviations of each class
# priors: vector of prior probabilities for each class

# e.g.: MOG for 100 observations from N(0,1) and 500 observations from N(10,2) with equal priors would be:
# generate_mog(c(100,500),c(0,10),c(1,2))

generate_mog <- function(n.data, means, sds, priors=NULL) {
  n.classes = length(means)
  if(is.null(priors)){ priors = rep(1/n.classes, n.classes) }
  if((length(means) != length(sds)) | (length(sds) != length(priors))){
    stop('Each class must have a mean, sd, and prior.')
  }
  classes = sample(1:n.classes, size=n.data, replace=TRUE, prob=priors)
  data = rnorm(n.data, mean = means[classes], sd = sds[classes])
  return (data)
}

EM_setup = function(data, n.classes, prior.choice='equal_priors_unif_params', ...){
	n_vars = ncol(data)
	if(prior.choice=='equal_priors_random_params'){
		prior = matrix(rep(1/n.classes, n.classes), nrow=n.classes) #All classes are equally as likely.
		mu = matrix(runif(n.classes*n_vars, 0, 1), nrow=n.classes) #All means are initially ~Unif[0,1].
		sigma = replicate(n.classes, matrix(runif(n_vars^2, 0, 1), nrow=n_vars))
	} else if(prior.choice=='means'){
	  #do something.
	} else {
		stop('Invalid prior choice.')
	}
	out = list(mu, sigma, prior)
	names(out) = c('mu','sigma','prior')
	return(out)
}

expectation = function(data, mu, sigma, prior){
	#Get the probabilities of each data point belonging to each class.
	n.classes = ncol(prior)
	n_vars = ncol(data)
	n.data = nrow(data)
	P.Xj.and.Ci = matrix(rep(NA, length=n.classes*n.data), ncol=n.classes)
	for(i in 1:n.classes){
		P.Ci.and.Xj[i] = prior[i] * dmvnorm(data, mu=mu[i,], sigma=sigma[,,i])
	}
	P.Ci.given.Xj = P.Xj.and.Ci / apply(P.Xj.and.Ci, 1, sum)
	
	most.likely.class = apply(P.Xj.and.Ci, 1, which.max)
	log.like = sum(log(diag(P.Ci.given.Xj[most.likely.class, 1:n.classes])))
	out = list(P.Ci.given.Xj, log.like)
	names(out) = c('P.Ci.given.Xj', 'log.like')
	return(out)
}

maximization = function(P.Ci.given.Xj, data, prior){
  n.classes = ncol(prior)
  n_vars = ncol(data)
  n.data = nrow(data)
	N = matrix(apply(P.Ci.given.Xj, 1, sum), ncol=1)
	mu = matrix(rep(NA, length=n.classes*n_vars, 0, 1), nrow=n.classes)
	sigma = replicate(n.classes, matrix(rep(NA, n_vars^2), nrow=n_vars))
	for(i in 1:n.classes){
		mu[i,] = colMeans(data * P.Ci.given.Xj[i] / N[i])
		diff = data - mu[rep(i, nrow(data)),]
		sigma[,,i] = diff %&% t(diff) / n.data
	}
	out = list(mu, sigma, prior)
	names(out) = c('mu','sigma','prior')
	return(out)
}

EM(data, n.classes, prior.choice, tol=.001, max_iters=1000, ...){
  init = EM_setup(data, n.classes, prior.choice, ...)
  mu = init$mu
  sigma = init$sigma
  prior = init$prior
  
  log.like = rep(-Inf, max_iters)
  for(i in 1:max_iters){
    exp = expectation(data, mu, sigma, prior)
    P.Ci.given.Xj = exp$P.Ci.given.Xj
    log.like[i] = exp$log.like
    
    max = maximization(P.Ci.given.Xj, data, prior)
    mu = max$mu
    sigma = max$sigma
    prior = max$prior
    
    if(abs(log.like) < tol){ break }
  }
  return(max)
}