import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import multivariate_normal
from sklearn import preprocessing 
import math
import time
import concurrent.futures
import multiprocessing as mp
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import sys
from os import makedirs
from errno import EEXIST

pool_number = int(sys.argv[1])

### SET RUN HYPERPARAMETERS ###

MOUSE = sys.argv[2].lower() == 'true'

# Set constants
#L = 0.8
#c = 10
G= 1

# settings
num_params = 17
hits_needed = 1000
cutoff1 = 5 # accepts like 40%
target_thresh = 0.2

prop = 0.60

prior_min = -3
prior_max = 3

gene_names=['IRX3', 'PAX6', 'OLIG2', 'NKX2-2']

if MOUSE:
    # MOUSE INITIAL CONDITIONS - if mouse == True
    X0 = [0.9673932616381912, 0.9905698391107712, 0.22032352020457502, 0.006276499439069383]
    t_range = [0.5, 7]
    timepts=[0.5,1, 1.5, 2, 2.5, 3, 3.5, 4,5,6,7]
    spec='mouse'

else:
    # HUMAN INITIAL CONDITIONS - if mouse == False
    X0 = [0.9978271229639294, 0.9976673744053118, 0.029296758259304223, 0.03117422425335577]
    t_range = [2,15]
    timepts=[2,4,5,6,7,8,9,10,11,12,13,14,15]
    spec='human'

# ODE functions
paramnames = ['r_s', 'gamma_s', 'K_NI', 'K_OI', 'K_NP', 'K_OP', 'K_IO', 'K_NO', 'K_GO', 'K_IN', 'K_PN', 'K_ON', 'K_GN', 'psiI', 'psiP', 'psiO', 'psiN']

# I P O N

def ODEs(t, X, r_s, gamma_s, K_NI, K_OI, K_NP, K_OP, K_IO, K_NO, K_GO, K_IN, K_PN, K_ON, K_GN, psiI, psiP, psiO, psiN):

	phi_I = (1/((1+((10**K_OI)*X[2])**2)*(1+((10**K_NI)*X[3])**2))+ (10**psiI))
	phi_P = (1/((1+((10**K_OP)*X[2])**2)*(1+((10**K_NP)*X[3])**2))+ (10**psiP))
	phi_O = ((((10**K_GO*G))**2)/((1+((10**K_IO)*X[0])**2)*(1+((10**K_NO)*X[3])**2)*(1+((10**K_GO)*G)**2)) + (10**psiO))
	phi_N = ((((10**K_GN*G))**2)/((1+((10**K_IN)*X[0])**2)*(1+((10**K_PN)*X[1])**2)*(1+((10**K_ON)*X[2])**2)*(1+((10**K_GN)*G)**2)) + (10**psiN))

	dI_dt = (10**r_s)*phi_I - (10**gamma_s)*X[0]
	dP_dt = (10**r_s)*phi_P - (10**gamma_s)*X[1]
	dO_dt = (10**r_s)*phi_O- (10**gamma_s)*X[2]
	dN_dt = (10**r_s)*phi_N- (10**gamma_s)*X[3]


	return np.array([dI_dt, dP_dt, dO_dt, dN_dt])


  # simulator function - only going to vary alpha and beta for now 
def simulate_data(parameter_list):
	run= solve_ivp(ODEs, y0=X0, t_span = t_range, args =tuple(parameter_list), t_eval = timepts, dense_output=True) 
	run_array = np.array(run.y).T
	return run_array 


### load data functions ###
def load_data(csvname, gene_names, mouse=True):
    df = pd.read_csv(csvname) # load csv
    df = pd.read_csv(csvname) # load csv
    df = df.loc[df['Gene'].isin(gene_names)] # select genes of interest only
    min_max_scaler = preprocessing.MinMaxScaler(feature_range=(0,1))

    # initialise empty variables
    start=0
    end=0
    timepts=[]

    if mouse:
        df = df[df.columns.drop(list(df.filter(regex='h_')))]
        df = df[df.columns.drop(list(df.filter(regex='m_0_')))]
        end=33

    else:
        df = df[df.columns.drop(list(df.filter(regex='m_')))]
        df = df[df.columns.drop(list(df.filter(regex='h_0_')))]
        df = df[df.columns.drop(list(df.filter(regex='h_1_')))]
        end=39
            
    means=[]

    for gene_id in gene_names:
        genearr = np.array(df[df.Gene == gene_id].drop(labels='Gene', axis =1))
        genes_scaled = min_max_scaler.fit_transform(genearr.reshape(len(genearr[0]),1)) # normalise data from 0 to 1
        genes_scaled = genes_scaled.reshape(1,len(genes_scaled))[0]
        rep1 = genes_scaled[[i for i in range(0,end,3)]]
        rep2 = genes_scaled[[i for i in range(1,end,3)]]
        rep3 = genes_scaled[[i for i in range(2,end,3)]]

        repmean=(rep1+rep2+rep3)/3
        means.append(repmean)
     
    run_array = np.array(means).T

    return run_array


###### ABC #######
def test_parset(parset, obs_data):
    go = simulate_data(parset)
    score = np.max((go-obs_data)**2)
    if np.sqrt(score) < cutoff1:
        return False
    else:
        return True

def get_successful_parset(obs_data):
	prior_space = np.random.uniform(prior_min, prior_max, 10000)
	miss = True 
	while miss:
		np.random.seed()
		param_list = [] 
		for _ in range(num_params):
			rand_prior = prior_space[np.random.randint(10000)]
			param_list.append(rand_prior)

		miss = test_parset(param_list, obs_data) # keeps rurnning until miss is False
#	print('hit')
	return param_list

def ABC_run(obs_data):
	print('ABC_run')
	good_params = []
	good_runs=[]
	good_scores = []

	def add_success(parset):
		go = simulate_data(parset)
		good_runs.append(go)
		score = np.max((go-observed)**2)
		good_scores.append(np.sqrt(score))
		good_params.append(parset) # now is a list of particles, each particle being 11 paramaters (or num_parameters)
		# print('parset', parset)

	pool = mp.Pool(processes = pool_number)

	for i in range(hits_needed):
		pool.apply_async(get_successful_parset, args=(obs_data,), callback=add_success) # call back function decides what to do with the return of the pooled function

	pool.close()
	pool.join()
	print('First ABC complete')
	return good_params, good_runs, good_scores

###### SMC #######

# test_parest function needs to have adjustable threshold
def test_SMCparset(parset, thresh, obs_data):
	observed = obs_data
	# print('parset', parset[0])
	go = simulate_data(parset[0])
#	score = np.sum((go-observed)**2, axis=0)
	score = np.max((go-observed)**2)
	if np.sqrt(score) < thresh:
		return False
	else:
		return True

# adjust get_successful_parset so that param varied based on previous set 
def get_successfulSMC_parset(old_params, thresh, prev_weig, prev_covar, obs_data):

	miss = True
	# np.random.seed(myidentity)
	while miss:
		np.random.seed()
		# print('simulating SMC...')
		redlight = 1
		while redlight != 0:
			redlight = 0
			num = np.random.choice(hits_needed, 1, p =prev_weig)[0] # so parameter set treated as one particle - try vectorise ?
			#print('random num & process:', num,mp.current_process().name, mp.current_process()._identity[0])
			old = old_params[num] # randomly selects a PARTICLE
			# multivatiate normal pertubation kernel
			mean = np.zeros(num_params)
		
			pertubation = np.random.multivariate_normal(mean=mean, cov = 2*prev_covar/10, size=1)		
			# print('ran num & proc:', num, mp.current_process().name, mp.current_process()._identity[0], pertubation)
			new = old + pertubation
			# print('new', new)
			for par in new[0]:
				if (par < prior_min) or (par > prior_max):
					redlight += 1

		miss = test_SMCparset(new, thresh, obs_data)

	### calculate weight for the chosen set
	denom = []
	for h in range(hits_needed):

		denom.append(prev_weig[h]*multivariate_normal.pdf(new[0], mean=old_params[h], cov=prev_covar, allow_singular=True))
	
	weight = (1/np.sum(denom))
	# print('weight', weight)

	return (list(new[0]), weight)

	
def SMCABC_run(good_params, good_runs, obs_data, good_scores):
    print('running SMC...')

    # initialise lists etc to store variables
    runs=good_runs
    good_scores.sort()
    final_good_scores=good_scores
    all_good_params={0: good_params} # key represents the layer, for each layer there will be hits_needed lists of num_params parameters (representing hits_needed particles)
    weig = {}
    weig[0]=[(1.0/hits_needed) for _ in range(hits_needed)]
#    max_good_scores ={0:good_scores}
    all_thresh = [cutoff1,final_good_scores[int(prop*hits_needed)]]
    threshold = all_thresh[1]

    #### need to set initial covar based on previous run
    covar = np.cov(all_good_params[0], rowvar=False)

    def add_SMCsuccess(param_tuple):
#        print('hitq
        param_list = param_tuple[0]
        weight = param_tuple[1]
        go = simulate_data(param_list)

        # calculate and add    new score
        score = np.max((go-observed)**2)
        final_good_scores.append(np.sqrt(score))
#        max_good_scores[l+1].append(np.max(np.sqrt(score)))

        # store run to plot
        runs.append(go)

        # store good particles
        all_good_params[l+1].append(param_list)

        # weights added
        weig[l+1].append(weight)

    l=0
    red_counter=0
    smallthresh = [0.3,0.25,0.2,0.0]
    proceed = True
    while proceed:
        if threshold == target_thresh:
            proceed=False
        final_good_scores=[]

        all_good_params[l+1] = []
#        max_good_scores[l+1]=[]

        #### adaptive weight ####
        for n in range(hits_needed):
            xkernel = multivariate_normal.pdf(observed.flatten(), runs[n].flatten()) # accounts for how close a run was to the data
            weig[l][n] = weig[l][n] * xkernel

        weig[l] = weig[l]/np.sum(weig[l]) # normalise again

        runs=[]
        weig[l+1] =[]

        pool = mp.Pool(processes = pool_number)

        for i in range(hits_needed):
            pool.apply_async(get_successfulSMC_parset, args=(all_good_params[l],threshold, weig[l], covar, obs_data), callback=add_SMCsuccess)

        pool.close()
        pool.join()

        print("Done Stage", l+2)

        # normalise weights
        weig[l+1] = weig[l+1]/np.sum(weig[l+1])
        covar = np.cov(all_good_params[l+1], rowvar=False)
        plot_trajectories(runs, threshold)
        save_csv(all_good_params[l+1], final_good_scores, threshold)

        final_good_scores.sort()

        threshold = final_good_scores[int(prop*hits_needed)]
        if (all_thresh[-1]-threshold) < 0.1: # to prevent getting stuck in local minima
            print('Thresholds converging, force reduce by 0.1')
            threshold = all_thresh[-1]-0.1

        if threshold < target_thresh:
            threshold = target_thresh

        elif threshold < target_thresh+0.1:
             threshold = smallthresh[red_counter]
             red_counter+=1

        all_thresh.append(threshold)
        print('Threshold schedule', all_thresh)
        print('\n')

        l+=1

    return runs, all_good_params[len(all_good_params)-1], final_good_scores, weig, threshold


###### Results processing ######

def plot_trajectories(run_output, thresh):
	run = run_output
	_, ax = plt.subplots()
	for i in range(len(run)):
	    ax.plot(timepts, run[i][:, 0], "--", color="peachpuff")
	    ax.plot(timepts, run[i][:, 1], "--", color="palegreen")
	    ax.plot(timepts, run[i][:, 2], "--", color="lightcoral")
	    ax.plot(timepts, run[i][:, 3], "--", color="cornflowerblue")

	ax.plot(timepts, observed[:, 0], "-", label="[I]", color = "orange")
	ax.plot(timepts, observed[:, 1], "-", label="[P]", color = "green")
	ax.plot(timepts, observed[:, 2], "-", label="[O]", color = "crimson")
	ax.plot(timepts, observed[:, 3], "-", label="[N]", color = "blue")

	ax.set_xlabel("time steps")
	ax.set_ylabel("concentration")
	ax.set_title("finding " + str(num_params) +" parameters")
	ax.legend(loc="upper right");
	plt.savefig("simulated_trajectories_"+spec+"/"+time.strftime("%Y_%m_%d_%H%M") + "_trajectories_at_thresh" + str(thresh) + ".png")


def save_csv(params, good_scores, thresh):
	# # load into pandas df
	param_df = pd.DataFrame(params)
	param_df.columns = ['param_' + str(n) for n in range(1, num_params+1)] # think about logging
	param_df['scores'] = good_scores
	param_df.to_csv("csvs_"+ spec +"/_"+time.strftime("%Y_%m_%d_%H%M") + "_final_params_thresh" + str(thresh) + ".csv")
	return param_df

### plot best parameters and get scores 

def get_posteriors(observed, params, good_scores, thresh):
	# # load into pandas df
	param_df = save_csv(params, good_scores, thresh)

	# # grid plot - separate out diff parameters
	print("Making gridplot...")
	g = sns.PairGrid(param_df, hue='scores', vars = ['param_' + str(n) for n in range(1,num_params+1)])
	g.map_diag(sns.kdeplot, hue= None)
	g.map_upper(sns.scatterplot, palette="icefire", s=10)
	g.map_lower(sns.scatterplot, palette="icefire", s=10)
	g.set(xlim=(prior_min, prior_max))
	g.set(ylim=(prior_min, prior_max))
	plt.savefig(spec+"_"+time.strftime("%Y_%m_%d_%H%M")+"_SMCgridplot_thresh" + str(thresh) + ".png")

	# # get best scores
	param_df.drop(labels='scores', inplace = True, axis =1)
	peaks =[]
	for col in param_df:
		ax = sns.kdeplot(param_df[col], label=str(col))
		points = ax.get_lines()[0].get_data()
		max_y = np.max(points[1])
		index = np.where(points[1] == max_y)
		max_x = points[0][index]
		peaks.append(max_x[0])
		ax.clear()

	unlogged = list(10**np.array(peaks))
	resultsdf = pd.DataFrame(list(zip(paramnames, peaks, unlogged)), columns = ['param', 'logged posterior', 'posterior'])
	print('\nposterior values:\n')
	print(resultsdf)

	simulateposterior = simulate_data(peaks)

	def get_cohen_score(go, obs):
	    diff = abs(go-obs)
	    scores = np.array([0,0,0,0])
	    for point in diff:
	        for gene in range(4):
	            if point[gene] > 0.2:
	                scores[gene] +=1
	                
	    return scores/len(diff)

	euclideandist_score = np.sqrt(np.sum((simulateposterior-observed)**2, axis=0))
	cohen_score = get_cohen_score(simulateposterior, observed)

	scoredf = pd.DataFrame(list(zip(gene_names, euclideandist_score, cohen_score)), columns = ['genes', 'euclidean distance', 'Cohen score'])
	print('\n posterior scores:\n')
	print(scoredf)

	### plot simulation from posterior against observed data ###

	sns.color_palette('colorblind')
	_, ax = plt.subplots()
	ax.plot(timepts, simulateposterior[:, 0], "--", color="orange")
	ax.plot(timepts, simulateposterior[:, 1], "--", color="green")
	ax.plot(timepts, simulateposterior[:, 2], "--", color="red")
	ax.plot(timepts, simulateposterior[:, 3], "--", color="cornflowerblue")

	ax.plot(timepts, observed[:, 0], "-", label="[I]", color = "orange")
	ax.plot(timepts, observed[:, 1], "-", label="[P]", color = "green")
	ax.plot(timepts, observed[:, 2], "-", label="[O]", color = "crimson")
	ax.plot(timepts, observed[:, 3], "-", label="[N]", color = "blue")

	ax.set_xlabel("time steps")
	ax.set_ylabel("concentration")
	ax.set_title(" varying K parameters")
	ax.legend(loc ='upper right');

	plt.savefig(spec +"_"+ time.strftime("%Y_%m_%d_%H%M") + "_posterior_trajectory"+ "_thresh" + str(thresh) + ".png")


if __name__ == '__main__':
    start_time = time.perf_counter()

    try:
        makedirs("simulated_trajectories_"+ spec)
    except OSError as exc:
        if exc.errno == EEXIST:
            pass

    try:
        makedirs("csvs_"+ spec)
    except OSError as exc:
        if exc.errno == EEXIST:
            pass

    observed = load_data(csvname='vst_transformed.csv', gene_names=gene_names, mouse = MOUSE)
    first_run = ABC_run(observed)
    SMC_run = SMCABC_run(first_run[0], first_run[1], observed, first_run[2])
    finish_time = time.perf_counter()
    time_taken = finish_time-start_time
    print("Time taken: %d s" %time_taken)
    print("Cores:", pool_number)
    plot_trajectories(SMC_run[0], SMC_run[4])
    get_posteriors(observed, SMC_run[1], SMC_run[2], SMC_run[4])




