list.of.packages = c('knitr', 'mvtnorm','MASS')
new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
rm(list.of.packages, new.packages)

# Generates MOG data given # of observations, means and standard deviatins of each Normal Distribution

# n.data: integer number of data points to generate
# means: vector of means for each class
# sigma: covariance matrix between the classes
# priors: vector of prior probabilities for each class

# e.g.: MOG for 100 observations from N(0,1) and 500 observations from N(10,2) with equal priors would be:
# generate.mog(c(100,500),c(0,10),c(1,2))

generate.mog <- function(n.data, mu, sigma, prior) {
  if(is.vector(mu)){ 
    n.classes = length(mu)
    set.seed(0)
    classes = sample(1:n.classes, size=n.data, replace=TRUE, prob=prior)
    data = rnorm(n=n.data, mean=mu[classes], sd=sigma[classes])
  } else { 
    n.classes = nrow(mu) 
    set.seed(0)
    classes = sample(1:n.classes, size=n.data, replace=TRUE, prob=prior)
    data = mvrnorm(n=n.data, mu=mu[classes,], Sigma=sigma[,,classes])
  }
  return(data)
}

EM.setup = function(data, n.classes, prior.choice='equal.priors.unif.params', ...){
  if(is.vector(data)){
    data = as.matrix(data, ncol=1)
  }
	n.vars = ncol(data)
	if(prior.choice=='equal.priors.random.params'){
		prior = rep(1/n.classes, n.classes) #All classes are equally as likely.
		mu = matrix(runif(n.classes*n.vars, 0, 1), nrow=n.classes) #All means are initially ~Unif[0,1].
		sigma = replicate(n.classes, matrix(runif(n.vars^2, 0, 1), nrow=n.vars))
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
	n.classes = length(prior)
	n.vars = ncol(data)
	n.data = nrow(data)
	P.Xj.and.Ci = matrix(rep(NA, length=n.classes*n.data), ncol=n.classes)
	if(n.vars==1){
	  for(i in 1:n.classes){
	    #P(Xj.and.Ci) = P(Ci) * P(Xj | Ci)
	    P.Xj.and.Ci[,i] = prior[i] * dnorm(data, mean=mu[i,],
	                                       sd=sigma[i])
	  }
	} else {
  	for(i in 1:n.classes){
  	  #P(Xj.and.Ci) = P(Ci) * P(Xj | Ci)
	    P.Xj.and.Ci[,i] = prior[i] * dmvnorm(data, mean=mu[i,], sigma=sigma[,,i])
	  }
	}
	P.Ci.given.Xj = P.Xj.and.Ci / apply(P.Xj.and.Ci, 1, sum)
	most.likely.class = apply(P.Xj.and.Ci, 1, which.max)
	log.like = sum(log(diag(P.Ci.given.Xj[1:n.data, most.likely.class])))
	out = list(P.Ci.given.Xj, log.like)
	names(out) = c('P.Ci.given.Xj', 'log.like')
	return(out)
}

maximization = function(P.Ci.given.Xj, data, prior){
  n.classes = length(prior)
  n.vars = ncol(data)
  n.data = nrow(data)
	N = apply(P.Ci.given.Xj, 2, sum)
	mu = matrix(rep(NA, length=n.classes*n.vars, 0, 1), nrow=n.classes)
	sigma = replicate(n.classes, matrix(rep(NA, n.vars^2), nrow=n.vars))
	for(i in 1:n.classes){
	  mu[i,] = sum(data * P.Ci.given.Xj[,i]) / N[i]
	  diff = data - mu[rep(i, nrow(data)),]
	  if(n.vars==1){
	    sigma[i] = sqrt(t(diff) %*% (diff * P.Ci.given.Xj[,i]) / N[i])
	  } else {
	    sigma[,,i] = sqrt(t(diff) %*% (diff * P.Ci.given.Xj[,i]) / N[i])
	  }
	}
	prior = N/sum(N)
	out = list(mu, sigma, prior)
	names(out) = c('mu','sigma','prior')
	return(out)
}

EM = function(data, n.classes, prior.choice, tol=1e-3, max.iters=500, ...){
  init = EM.setup(data, n.classes, prior.choice, ...)
  mu = init$mu
  sigma = init$sigma
  prior = init$prior
  
  if(is.vector(prior)){
    mu = as.matrix(mu, nrow=1)
    data = as.matrix(data, ncol=1)
  }
  
  log.like = rep(-Inf, max.iters)
  for(i in 1:max.iters){
    exp = expectation(data, mu, sigma, prior)
    P.Ci.given.Xj = exp$P.Ci.given.Xj
    log.like[i] = exp$log.like
    
    max = maximization(P.Ci.given.Xj, data, prior)
    mu = max$mu
    sigma = max$sigma
    prior = max$prior
    
    #For debugging
    print(c('iter: ', as.character(i)))
    print(mu)
    print(sigma)
    print(prior)
    
    if(abs(log.like[i]) < tol){ break }
  }
  log.like = log.like[1:i]
  out = list(mu, sigma, prior, log.like)
  names(out) = c('mu','sigma','prior', 'log.like')
  return(out)
}