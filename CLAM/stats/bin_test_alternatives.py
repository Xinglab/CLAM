""" A script for storing alternative peak calling
statistical models / bin tests other than negative binomial
Zijun Zhang
Last revisited: 9.6.2017
"""


def test_bin_poisson(intv_bin_ip, intv_bin_con, correction_method='fdr_bh'):
	"""DOCSTRING
	Args
	Returns
	"""
	def _par_to_vec(par, data, is_constrained):
		if is_constrained:
			beta = par[0]
			mu_vec = par[1::]
			delta = 0
		else:
			beta, delta = par[0], par[1]
			mu_vec = par[2::]
		ip_counter = data['this_ip'].shape[0]
		con_counter = data['this_con'].shape[0]
		mu0 = np.asarray(mu_vec[0:con_counter])
		mu1 = np.asarray(mu_vec[con_counter::])
		lamb1_this = np.exp(mu1 + beta + delta)
		lamb1_others = np.exp(mu1)
		lamb0_this = np.exp(mu0 + beta)
		lamb0_others = np.exp(mu0)
		return (lamb1_this, lamb1_others, lamb0_this, lamb0_others)
		
	def _neg_loglik_unconstrain(par, data):
		(l1, l2, l3, l4) = _par_to_vec(par, data, False)
		ll = np.sum(poisson.logpmf(data['this_ip'], mu=l1)) + \
			np.sum(poisson.logpmf(data['others_ip'], mu=l2)) + \
			np.sum(poisson.logpmf(data['this_con'], mu=l3)) + \
			np.sum(poisson.logpmf(data['others_con'], mu=l4))
		return -ll
	
	def _neg_loglik_constrain(par, data):
		(l1, l2, l3, l4) = _par_to_vec(par, data, True)
		ll = np.sum(poisson.logpmf(data['this_ip'], mu=l1)) + \
			np.sum(poisson.logpmf(data['others_ip'], mu=l2)) + \
			np.sum(poisson.logpmf(data['this_con'], mu=l3)) + \
			np.sum(poisson.logpmf(data['others_con'], mu=l4))
		return -ll
		
	intv_counter = intv_bin_ip.shape[1]
	assert intv_counter == intv_bin_con.shape[1]
	binscore = np.empty(intv_counter)
	binsignal = np.empty(intv_counter)
	ip_sum = np.apply_along_axis(np.sum, 1, intv_bin_ip)
	con_sum = np.apply_along_axis(np.sum, 1, intv_bin_con)
	for i in range(intv_counter):
		this_ip = intv_bin_ip[:, i]
		others_ip = ip_sum - this_ip
		this_con = intv_bin_con[:, i]
		others_con = con_sum - this_con
		if this_ip == 0:
			binsignal[i], binscore[i] = np.nan, 1.0
			continue
		## because Poisson (and other count-based methods) only
		## takes integers, here we take the floor of the fractional
		## multi-reads as a conservative approach
		data = {
				'this_ip':np.floor(this_ip),
				'others_ip':np.floor(others_ip),
				'this_con':np.floor(this_con),
				'others_con':np.floor(others_con)
			}
		
		res_constrain = optimize.minimize(
				x0=np.ones(1+this_ip.shape[0]+others_ip.shape[0]), 
				fun=_neg_loglik_constrain,
				args=(data),
				method='Nelder-Mead',
				options={'disp':False}
			)
		
		res_unconstrain = optimize.minimize(
				x0=np.ones(2+this_ip.shape[0]+others_ip.shape[0]), 
				fun=_neg_loglik_unconstrain,
				args=(data),
				method='Nelder-Mead',
				options={'disp':False}
			)
		
		delta_mle = res_unconstrain.x[1]
		pval = 1 - chi2.cdf(2*(res_constrain.fun - res_unconstrain.fun), 1)
		binscore[i] = pval
		binsignal[i] = delta_mle
	adj = multipletests(binscore, alpha=0.05, method=correction_method)
	binscore_adj = adj[1]
	return binsignal, binscore_adj


def test_bin_fisher(intv_bin_ip, intv_bin_con, with_control=True, correction_method='fdr_bh'):
	"""DOCSTRING
	Args
	Returns
	"""
	if intv_bin_ip.shape[0] != 1:
		raise Exception('Fisher exact test does not deal with replicates.')
	intv_counter = intv_bin_ip.shape[1]
	assert intv_counter == intv_bin_con.shape[1]
	binscore = np.empty(intv_counter)
	binsignal = np.empty(intv_counter)
	ip_sum = np.sum(intv_bin_ip[0,])
	con_sum = np.sum(intv_bin_con[0,])
	for i in range(intv_counter):
		this_ip = intv_bin_ip[0, i]
		others_ip = ip_sum - this_ip
		this_con = intv_bin_con[0, i]
		others_con = con_sum - this_con
		if this_ip == 0:
			binsignal[i], binscore[i] = np.nan, 1.0
			continue
		_, binscore[i] = fisher_exact([[this_ip, others_ip], [this_con, others_con]], alternative='greater')
		if with_control:
			binsignal[i] = this_ip/others_ip / this_con*others_con
		else:
			binsignal[i] = this_ip
		
	adj = multipletests(binscore, alpha=0.05, method=correction_method)
	binscore_adj = adj[1]
	return binsignal, binscore_adj
