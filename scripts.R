list.of.packages = c('knitr', 'mvtnorm')
new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
rm(list.of.packages, new.packages)

initial_setup = function(n_classes, prior_choice='equal_priors_unif_params', data){
	n_vars = ncol(data)
	if(prior_choice=='equal_priors_random_params'){
		prior = matrix(rep(1/n_classes, n_classes), nrow=n_classes) #All classes are equally as likely.
		mu = matrix(runif(n_classes*n_vars, 0, 1), nrow=n_classes) #All means are initially ~Unif[0,1].
		sigma = replicate(n_classes, matrix(runif(n_vars^2, 0, 1), nrow=n_vars))
	} else {
		stop('invalid prior choice.')
	}
	return(list(prior, params)))
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
	N = matrix(apply(P.Ci.given.Xj, 1, sum), ncol=1)
	mu = matrix(rep(NA, length=n_classes*n_vars, 0, 1), nrow=n_classes)
	for(i in 1:n_classes){
		mu[i,] = colMeans(data * P.Ci.given.Xj[i] / N[i])
		sigma[,,i] = 
	}
	prior = N / sum(N)
	mu = P.Ci.given.Xj * data

}
        sigma <- (Y %*% t(Y) + diag(extra_var)) / n

function [new_params,new_priors] = maxim(full_vector,data,priors)
    new_params = []
    new_priors = []
    for i = 1:length(priors)
        total=0
        total_squares=0
        points_distribution = sum(full_vector)
        points = points_distribution(i)
        for j = 1:length(data)
            total = total + full_vector(j,i)*data(j)
            total_squares = total_squares + full_vector(j,i)*data(j)^2
        end
        mean = total/points
        sd = sqrt((total_squares/points)-(mean^2))
        prior = points/sum(points_distribution)
        new_params = [new_params,mean,sd]
        new_priors = [new_priors,prior]



prior_and_params = initial_setup(3)
full_vector = expec(dat,prior,params)
[params,prior] = maxim(full_vector,dat,prior)

for i = 1:15
    full_vector = expec(dat,prior,params)
    [params,prior] = maxim(full_vector,dat,prior)
end