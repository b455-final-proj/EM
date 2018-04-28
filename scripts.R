list.of.packages = c('knitr', 'mixtools','MASS', 'Matrix', 'stats')
new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
rm(list.of.packages, new.packages)

generate.mog <- function(n.data, mu, sigma, prior) {
  if(is.vector(mu)){ 
    n.classes = length(mu)
    set.seed(0)
    classes = sample(1:n.classes, size=n.data, replace=TRUE, prob=prior)
    data = rnorm(n=n.data, mean=mu[classes], sd=sigma[classes])
  } else { 
    n.classes = nrow(mu) 
    n.vars = ncol(mu)
    set.seed(0)
    classes = sample(1:n.classes, size=n.data, replace=TRUE, prob=prior)
    classes = c(0, table(classes))
    data = matrix(NA, nrow=n.data, ncol=n.vars)
    for(i in 1:n.classes){
      data[(sum(classes[1:i])+1):(sum(classes[1:i])+classes[i+1]),] = mvrnorm(
        n=classes[i+1], mu = mu[i,], Sigma=sigma[,,i])
    }
  }
  return(data)
}

EM.setup = function(data, n.classes, prior.choice='equal.priors.unif.params', ...){
  if(is.vector(data)){
    data = as.matrix(data, ncol=1)
  }
	n.vars = ncol(data)
	
	if(prior.choice=='equal.priors.random.params'){
	  generate.sigma = function(nrow=n.vars){
	    A = matrix(runif(n.vars^2, 0, 1), nrow=n.vars)
	    return(t(A) %*% A)
	  }
		prior = rep(1/n.classes, n.classes) #All classes are equally as likely.
		mu = matrix(runif(n.classes*n.vars, 0, 1), nrow=n.classes) #All means are initially ~Unif[0,1].
		sigma = replicate(n.classes, generate.sigma(nrow=n.vars))
		
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
	    P.Xj.and.Ci[,i] = prior[i] * dmvnorm(data, mu=mu[i,], sigma=sigma[,,i])
	  }
	}
	#Infinite values are returned at very low densities. Replace them with a low value.
	P.Xj.and.Ci[which(is.infinite(P.Xj.and.Ci), arr.ind=TRUE)] = min(P.Xj.and.Ci) * 1e-10
	P.Ci.given.Xj = P.Xj.and.Ci / apply(P.Xj.and.Ci, 1, sum)
	N = apply(P.Ci.given.Xj, 2, sum)
	log.like = -sum(log(apply(P.Xj.and.Ci, 1, sum)))
	out = list(P.Ci.given.Xj, N, log.like)
	names(out) = c('P.Ci.given.Xj', 'N', 'log.like')
	return(out)
}

maximization = function(P.Ci.given.Xj, N, data, prior){
  n.classes = length(prior)
  n.vars = ncol(data)
  n.data = nrow(data)
	mu = matrix(rep(NA, length=n.classes*n.vars, 0, 1), nrow=n.classes)
	sigma = replicate(n.classes, matrix(rep(NA, n.vars^2), nrow=n.vars))
	for(i in 1:n.classes){
	  mu[i,] = t(data) %*% P.Ci.given.Xj[,i] / N[i]
	  if(n.vars==1){
	    sigma[i] = sqrt(t(data^2) %*% P.Ci.given.Xj[,i] / N[i] - mu[i,]^2)
	  } else {
	    diff = data - apply(data, 2, mean)
	    sigma[,,i] = sqrt(t(diff) %*% (diff * P.Ci.given.Xj[,i]) / N[i])
	  }
	}
	prior = N/sum(N)
	out = list(mu, sigma, prior)
	names(out) = c('mu','sigma','prior')
	return(out)
}

EM = function(data, n.classes, prior.choice, tol=1e-5, max.iters=500, ...){
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
    N = exp$N
    
    max = maximization(P.Ci.given.Xj, N, data, prior)
    mu = max$mu
    sigma = max$sigma
    prior = max$prior
    
    if(i!=1){ if(log.like[i - 1] - log.like[i] < tol){ break } }
  }
  log.like = log.like[1:i]
  out = list(mu, sigma, prior, log.like)
  names(out) = c('mu','sigma','prior', 'log.like')
  return(out)
}