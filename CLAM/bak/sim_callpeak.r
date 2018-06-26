## simulate read counts in bins and
## perform LRT test as peak calling
## *-- prototyping --*
## Zijun Zhang
## 9.1.2017


sim_bin_counts = function(mu_vec, beta, delta)
{
	others = c(
		rpois(n=1, lambda=exp(mu_vec[1])),
		rpois(n=1, lambda=exp(mu_vec[2]))
		)
	this = c(
		rpois(n=1, lambda=exp(mu_vec[1]+beta+delta)),
		rpois(n=1, lambda=exp(mu_vec[2]+beta))
		)
	res = matrix(c(this,others), nrow=2, byrow=T)
	rownames(res) = c('this','others')
	colnames(res) = c('IP','Input')
	return(as.data.frame(res))
}


loglik_constrain = function(par, data)
{
	ll = 0
	mu1=par[1]; mu0=par[2]; beta = par[3]
	lamb1.this = exp(mu1 + beta)
	lamb1.others = exp(mu1)
	lamb0.this = exp(mu0+beta)
	lamb0.others = exp(mu0)
	ll = ll + dpois(data['this','IP'], lamb1.this, log=T)
	ll = ll + dpois(data['this','Input'], lamb0.this, log=T)
	ll = ll + dpois(data['others','IP'], lamb1.others, log=T)
	ll = ll + dpois(data['others','Input'], lamb0.others, log=T)
	return(ll)
}

loglik_unconstrain = function(par, data)
{
	ll = 0
	mu1=par[1]; mu0=par[2]; beta = par[3]; delta=par[4]
	lamb1.this = exp(mu1 + beta + delta)
	lamb1.others = exp(mu1)
	lamb0.this = exp(mu0+beta)
	lamb0.others = exp(mu0)
	ll = ll + dpois(data['this','IP'], lamb1.this, log=T)
	ll = ll + dpois(data['this','Input'], lamb0.this, log=T)
	ll = ll + dpois(data['others','IP'], lamb1.others, log=T)
	ll = ll + dpois(data['others','Input'], lamb0.others, log=T)
	return(ll)
}


callpeak_LRT = function(data)
{
	ll0 = optim(rep(1,3), loglik_constrain, control=list(fnscale=-1), data=data)
	ll1 = optim(rep(1,4), loglik_unconstrain, control=list(fnscale=-1), data=data)
	pval = 1-pchisq(2*(ll1$value-ll0$value),1)
	pval
}



B=200
res = matrix(NA, nrow=B, ncol=2)
colnames(res) = c('fisher', 'lrt')
for(b in 1:B) {
	data = sim_bin_counts(c(2.5,2), -0.5, 1)
	p1=fisher.test(data)$p.value
	p2=callpeak_LRT(data)
	res[b,] = c(p1, p2)
}

plot(res[,'fisher'], res[,'lrt'])
abline(0,1)
mean(res[,'fisher']<0.05)
mean(res[,'lrt']<0.05)
