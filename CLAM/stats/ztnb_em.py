import scipy.special as special
import numpy as np
from numpy.random import negative_binomial
import mpmath
from collections import defaultdict


##############################
## distribution characterizing functions
##############################

def trunc_logLik(data, mu, alpha):
	log_1_plus_a_mu = np.log(1 + alpha*mu)
	log_1_minus_prob_zero = np.log(1.0 - np.exp(-np.log(1.0+alpha*mu)/alpha))
	alpha_inv = 1.0/alpha
	lim = int(np.max(data.keys()))
	holding_val=0.0
	log_L=0.0
	for i in range(1, lim+1):
		holding_val += np.log(1+alpha*(i-1))
		log_L += data[i]* (holding_val - special.gammaln(i) + i*np.log(mu)-(i+alpha_inv)*log_1_plus_a_mu - log_1_minus_prob_zero)
	return log_L

def ztnb_pmf(y, mu, alpha):
	r = 1.0 / alpha
	if y <= 0:
		raise Exception('y must be larger than 0.')
	p = mu/(mu+r+0.0)
	ztnbin_mpmath = lambda y, p, r: mpmath.gamma(y + r)/(mpmath.gamma(y+1)*mpmath.gamma(r))*np.power(1-p, r)*np.power(p, y)/(1-np.power(1-p, r))
	ztnbin = np.frompyfunc(ztnbin_mpmath, 3, 1)
	return float(ztnbin(y, p, r))

def ztnb_cdf(y, mu, alpha):
	r = 1.0/alpha
	if y <= 0:
		raise Exception('y must be larger than 0.')
	p = mu/(mu+r+0.0)
	F_ztnb = ( 1 - special.btdtr(y+1, r, p) - np.power(1-p, r) ) / (1-np.power(1-p,r))
	return F_ztnb

def ztnb_pval(y, mu, alpha):
	pval = 1 - ztnb_cdf(y, mu, alpha) + ztnb_pmf(y, mu, alpha)
	if pval <= 10**-5:
		return 0
	else:
		return pval

def rztnb(mu=3, alpha=0.5, size=100):
	r = 1.0/alpha
	p = mu/(mu+r+0.0)
	ztnb=[]
	while(len(ztnb)<size):
		x = negative_binomial(n=r,p=1-p)
		if x>0:
			ztnb.append(x)
	return ztnb

def collapse_data(data):
	col_data = defaultdict(int)
	for i in data:
		col_data[i] += 1
	return col_data	

##############################
## parameter estimation functions	
##############################

def EM_estim_params(height, tol = 10**-4, max_iter = 1000, verbose = False, mu = None, alpha = None):
	tot_size = np.sum([height[x] for x in height if x>0])
	error = 10000
	prev_score = 10000
	score = 0.0
	if mu is None or alpha is None:
		mu = np.sum([x*height[x] for x in height])/(tot_size+0.0)
		var = np.sum([ height[x]*(x-mu)**2 for x in height]) / tot_size
		alpha = (var - mu) / (mu * mu)
	#mu = 0.01
	#alpha = 100
	for h in range(1, int(max(height.keys()))+1):
		if not h in height:
			height[h]=0
	
	for i in range(1, (max_iter+1)):
		height[0] = expected_zeros(tot_size, mu, alpha)
		mu, alpha = estim_params(height, tol)
		height[0] = 0
		score = trunc_logLik(height, mu, alpha)
		if score == 0:
			raise ZeroDivisionError('invalid loglik function value')
		error = abs((score - prev_score)/score)
		if verbose:
			print('Iter ' + str(i) + ': eps = ' + str(error) + '; mu = ' + str(mu) + '; alpha = ' + str(alpha))
		if(error < tol):
			break
		prev_score = score
	return (trunc_logLik(height, mu, alpha), mu, alpha)

def expected_zeros(pseudo_size, mu, alpha):
	min_allowed_alpha=10**-4
	max_allowed_prob_zero=0.99
	if alpha < min_allowed_alpha:
		prob_zero = max_allowed_prob_zero
	else:
		prob_zero = np.min([np.power(1.0+alpha*mu, -1.0/alpha), 0.99])
	expected_zeros = int(pseudo_size*(prob_zero/(1-prob_zero)))
	return expected_zeros


def estim_params(pseudo_hist, tolerance = 10**-4):
	min_allowed_alpha = 10**-3
	max_allowed_alpha = 1000
	
	mu = compute_mean(pseudo_hist)  
	pseudo_size = np.sum(pseudo_hist.values())
  
	a_low = min_allowed_alpha
	a_high = max_allowed_alpha
  
	diff = 10000
	prev_val = 10000
  
	while diff > tolerance and movement(a_high, a_low) > tolerance:
		a_mid = (a_low + a_high)/2
		mid_val = alpha_score_function(pseudo_hist, mu, a_mid, pseudo_size)
		#print str(a_mid) + '; ' + str(mid_val) + '; ' + str(trunc_logLik(pseudo_hist, mu, a_mid))
		if (mid_val < 0): 
			a_high = a_mid
		else:
			a_low = a_mid
		diff = np.abs((prev_val - mid_val)/prev_val)
		prev_val = mid_val
	
	alpha = a_mid
	return mu, alpha

def alpha_score_function(vals_hist, mean, a_mid, vals_count):
	one_plus_alpha_mu = 1.0 + a_mid*mean
	return (score_fun_first_term(vals_hist, a_mid)/(vals_count+0.0) + (np.log(one_plus_alpha_mu)/(a_mid+0.0) - mean)/(a_mid+0.0))

def score_fun_first_term(vals_hist,a_mid):
	sum = 0.0
	lim = int(np.max(vals_hist.keys()))
	for i in range(0, lim+1):
		if (vals_hist[i] > 0):
			inner_sum = 0.0
			for j in range(0, i):
				inner_sum += j/(1.0 + a_mid*j)
			sum += vals_hist[i]*inner_sum
    
	return sum 

	
##############################	
## in-line functions
##############################

def compute_mean(height):
	tot_size = np.sum(height.values())
	mean = np.sum([x*height[x] for x in height])/(tot_size + 0.0)
	return(mean)

def movement(a, b):
	return abs(a - b)/max(a, b)

##############################
## testing function
##############################

def test(size=10**3, mu=0.01, alpha=50, max_iter=100):
	data=rztnb(mu, alpha, size)
	height=collapse_data(data)
	return EM_estim_params(height, max_iter=max_iter, verbose=True)