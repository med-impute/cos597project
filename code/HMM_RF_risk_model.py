import sklearn as sk
import cPickle as pickle
from sklearn.hmm import GaussianHMM
from sklearn.decomposition import PCA
from sklearn.metrics import roc_curve, auc
from sklearn.ensemble import RandomForestClassifier
import numpy as np
import pandas as pd
import sys, os, re
from math import *
#from med_data import *
#from medutils_new import *
#from model_utils import *
import matplotlib.pyplot as plt


model_template = { 
	'keep_ratio_half' : True,
	'ntest_alive' : 500,
	'ntest_dead' : 500,
	'Npatients' : 1000,
	'min_duration' : 120.,
	'max_time_before_end' : 120.,
	'min_time_before_end_for_testing' : 0.,
	'min_nobs' : 6,
	'min_npats' : 10,
	'nstates' : 3,
	'covariance_type' : 'full',
	'n_estimators' : 50
}
models_to_do = {
	'deltat_0' : { 'min_time_before_end_for_testing' : 0. },
	'deltat_6' : { 'min_time_before_end_for_testing' : 6. },
	'deltat_12' :  { 'min_time_before_end_for_testing' : 12. },
	#'deltat_18' : { 'min_time_before_end_for_testing' : 18. },
	'deltat_24' : { 'min_time_before_end_for_testing' : 24. },
	#'deltat_30' : { 'min_time_before_end_for_testing' : 30. },
	#'deltat_36' : { 'min_time_before_end_for_testing' : 36. },
	#'deltat_42' : { 'min_time_before_end_for_testing' : 42. },
	'deltat_48' : { 'min_time_before_end_for_testing' : 48. }
}

for m in models_to_do:
	for par in model_template:
		if par not in models_to_do[m]:
			models_to_do[m][par] = model_template[par]


vitals_to_use = [ u'J02_RR', u'J03_HR', u'J04_DBP', u'J05_SBP', 
                 u'J06_TEMP', u'J07_AnionGap', u'J08_BUN', u'J09_CO2', 
                 u'J10_CalciumLvl', u'J11_Chloride', u'J12_Creatinine', 
                 u'J13_GlucoseLvl', u'J14_GlucoseLvlPOC', u'J15_Hct', 
                 u'J16_Hgb', u'J17_MCH', u'J18_MCHC', u'J19_MCV', u'J20_INR', 
                 u'J21_PT', u'J22_PTT', u'J23_Platelet', u'J24_PotassiumLvl', 
                 u'J25_RBC', u'J26_RDW', u'J27_SodiumLvl', u'J28_WBC', u'J29_pCO2Art']
def save_obj(obj, filename):
    with open(filename,'wb') as ofile:
        pickle.dump(obj, ofile, pickle.HIGHEST_PROTOCOL)
        
def load_obj(filename):
    with open(filename, 'rb') as infile:
        return pickle.load(infile)
   
class HMM_and_RF:
	def __init__(self, min_nobs = 6, min_npats = 10, nstates_all = 3, 
		covariance_type = "full",  n_estimators=30, vitals=vitals_to_use):
		self.min_nobs = min_nobs
		self.min_npats = min_npats
		self.nstates_all = nstates_all
		self.covariance_type = covariance_type
		
		self.nstates_alive = {}
		self.nstates_dead = {}
		self.observations_dead = {}
		self.observations_alive = {}

		for v in vitals_to_use:
			self.nstates_alive[v] = nstates_all
		for v in vitals_to_use:
			self.nstates_dead[v] = nstates_all
		self.observations = None
		self.vitals_to_use = vitals
		self.n_estimators = n_estimators
		self.classifier= RandomForestClassifier(n_estimators=self.n_estimators)
		

	#def fit_Nstates(self, apans, dpans):
	#	max_nstates = 7
	#	for v in self.vitals_available:
	#		max_score = None
	#		for nd in range(1, max_nstates + 1):
	#			
	#	for nd in range(1,max_nstates+1):
	#		for na in range(1, max_nstates + 1):		
	#
	#
	#
	#
	#
	#
	def get_observations(self, apans, dpans):
		self.pans = []
		self.apans= apans
		self.dpans = dpans
		self.observations_alive = {}
		self.observations_dead = {}

		
		for v in self.vitals_to_use: 
		    self.observations_alive[v] = []
		    self.observations_dead[v] = []

		print "Getting data for %d alive patients"%(len(apans))
		for pan in self.apans:
		    observation = load_data(pan)
		    for v in self.vitals_to_use:
		        if len(observation[v]) >= self.min_nobs:
		        	self.observations_alive[v].append(np.array(observation[v]))

		print "Getting data for %d dead patients"%(len(dpans))
		for pan in self.dpans:
		    observation = load_data(pan)
		    for v in self.vitals_to_use:
		        if len(observation[v]) >= self.min_nobs:
		        	self.observations_dead[v].append(np.array(observation[v]))

		print "Finding the vitals for which there is enough data"
		self.vitals_available = []
		for v in self.vitals_to_use:
			if len(self.observations_alive[v]) >= self.min_npats and\
			   len(self.observations_dead[v]) >= self.min_npats:
				self.vitals_available.append(v)

		print len(self.vitals_available)," vitals available"
	def fit_HMMs(self, apans=None, dpans=None):
		if apans is not None: self.get_observations(apans, dpans)
		# gather data

		self.HMMs_dead = {}
		self.HMMs_alive = {}
		self.risk_vectors_dead  = []
		self.risk_vectors_alive = []

		print "Training HMM's"
		for v in self.vitals_available:
			self.HMMs_dead[v] = GaussianHMM(self.nstates_dead[v], self.covariance_type ).fit(self.observations_dead[v])
			self.HMMs_alive[v] = GaussianHMM(self.nstates_alive[v], self.covariance_type ).fit(self.observations_alive[v])
	def get_risk_vectors(self,pans):

		risk_vectors = []
		for pan in pans:
		    observation = load_data(pan, testing=True)
		    risk_vector = []
		    for v in self.vitals_available:
		    	observation[v] = np.array(observation[v])
		    	if len(observation[v]) < self.min_nobs: risk_vector.extend([0.0, 0.0])
		    	else: 
		    		risk_vector.append(self.HMMs_dead[v].score(observation[v]))
		    		risk_vector.append(self.HMMs_alive[v].score(observation[v]))
		    risk_vectors.append(risk_vector)
		return risk_vectors
	def get_training_risk_vectors(self, apans=None, dpans=None):
		
		if apans is None:
			apans = self.apans
		if dpans is None:
			dpans = self.dpans 

		print "Building risk vectors (alive)"
		for pan in apans:
		    observation = load_data(pan)
		    risk_vector = []
		    for v in self.vitals_available:
		    	observation[v] = np.array(observation[v])
		    	if len(observation[v]) < self.min_nobs: risk_vector.extend([0.0, 0.0])
		    	else: 
		    		risk_vector.append(self.HMMs_dead[v].score(observation[v]))
		    		risk_vector.append(self.HMMs_alive[v].score(observation[v]))
		    self.risk_vectors_alive.append(risk_vector)

		print "Building risk vectors (dead)"
		for pan in dpans:
		    observation = load_data(pan)
		    risk_vector = []
		    for v in self.vitals_available:
		    	observation[v] = np.array(observation[v])
		    	if len(observation[v]) < self.min_nobs: risk_vector.extend([0.0, 0.0])
		        else: 
		        	risk_vector.append(self.HMMs_dead[v].score(observation[v]))
		        	risk_vector.append(self.HMMs_alive[v].score(observation[v]))
		    self.risk_vectors_dead.append(risk_vector)

		self.risk_vectors_all = []
		self.risk_vectors_all.extend(self.risk_vectors_dead)
		self.risk_vectors_all.extend(self.risk_vectors_alive)
		
		self.risk_labels = []
		self.risk_labels.extend([ 1 for rv in self.risk_vectors_dead ])
		self.risk_labels.extend([ 0 for rv in self.risk_vectors_alive ])

	def fit_classifier(self):
		print "Fitting Random Forest classifier"
		inds = np.arange(0,len(self.risk_vectors_all))
		np.random.shuffle(inds)

		risk_vectors_shuffled = np.array([ self.risk_vectors_all[i] for i in inds ])
		labels_shuffled = np.array([ self.risk_labels[i] for i in inds ])

		#print risk_vectors_shuffled.shape, labels_shuffled.shape
		self.classifier.fit(risk_vectors_shuffled, labels_shuffled)

	def fit(self, apans=None, dpans=None):

		self.fit_HMMs(apans,dpans)
		self.get_training_risk_vectors(apans,dpans)
		self.fit_classifier()

	def predict_proba(self, pans):
		risk_vectors = self.get_risk_vectors(pans)
		return self.classifier.predict_proba(risk_vectors)


    
# load patient list -- you need to check which of these are Final and Admitting heart condition diagnoses.
pan_list = np.loadtxt('../hypertensive_and_ischemic_heart_disease_pans.dat', dtype='S50')
pan_list = np.unique([ int(p) for p in pan_list ])

# load dictionaries
pan_to_index = load_obj('pan_to_index.dict')
pan_dead_alive = load_obj('pan_dead_alive.dict')

pan_final_codes = load_obj('pan_final_codes.dict.list')
pan_admitting_codes = load_obj('pan_admitting_codes.dict.list')

# load lists of alive/dead pans
alive_heart_pans = load_obj('alive_heart_pans.list')
dead_heart_pans = load_obj('dead_heart_pans.list')
fate_unknown_heart_pans = load_obj('fate_unknown_heart_pans.list')

def patient_file(pan):
    return "/data/penn_new/lifangc/single/patient_%06d.csv"%(pan_to_index[pan])




for model_name in models_to_do:
	model = models_to_do[model_name]
	keep_ratio_half = model['keep_ratio_half']
	ntest_alive = model['ntest_alive']
	ntest_dead = model['ntest_dead']
	Npatients = model['Npatients']
	max_time_before_end = model['max_time_before_end']
	min_time_before_end_for_testing = model['min_time_before_end_for_testing']
	min_duration = model['min_duration']
	min_nobs = model['min_nobs']

		
	def load_data(pan, testing=False, tmax=None):
	    data_raw = pd.read_csv(patient_file(pan))
	    
	    # Select times that are within "max_time_before_end" hours from the end of the hospital visit
	    times = data_raw['J01_TAKEN_DATE'].values
	    times -= times[0]
	    times = np.array([ float(t) /3600. for t in times ]) 
	    if tmax is None: tmax = max_time_before_end
	    if testing:
	    	okinds = np.array([ i for i in np.arange(0, len(times)) if ((max(times) - times[i]) < tmax and (max(times) - times[i]) > min_time_before_end_for_testing ) ])
	    else:
	    	okinds = np.array([ i for i in np.arange(0, len(times)) if (max(times) - times[i]) < tmax ])
	    try:
	    	data_raw = data_raw.iloc[okinds]
	    except:
	    	print "LOAD DATA ERROR", okinds
	    
	    # Translate to a dictionary
	    data_final = {}
	    for v in vitals_to_use:
	        data_final[v] = []
	        for d in data_raw[v].values:
	            if not np.isnan(d):
	                data_final[v].append([d])
	    return data_final

	def pan_is_ok(pan, testing=False, tmax=None):
		# Can we read the bloody file?
	    try:
	        data = pd.read_csv(patient_file(pan))
	    except:
	        return False
	    
	    # Yay we can read the file. But are there enough observations in the time window?
	    times = data['J01_TAKEN_DATE'].values
	    times -= times[0]
	    times = np.array([ float(t) /3600. for t in times ]) 
	    if tmax is None: tmax = max_time_before_end
	    if testing:
	    	okinds = np.array([ i for i in np.arange(0, len(times)) if ((max(times) - times[i]) < tmax and (max(times) - times[i]) > min_time_before_end_for_testing) ])
	    else:
	    	okinds = np.array([ i for i in np.arange(0, len(times)) if (max(times) - times[i]) < tmax ])
	    
	    if len(okinds) < min_nobs: return False

	    try:
	    	getdat = data.iloc[okinds]
	    except:
	    	print "Can't select ",okinds,"indices from data!"
	    	return False
	    

	    # What about the duration -- were they in the hospital for long enough?
	    if max(times) < min_duration: return False
	    
	    # OK, what about this: is there at least one vital for which there are enough observations?
	    data = data.iloc[okinds]
	    
	    nvals = {}
	    for v in vitals_to_use:
	        nvals[v] = len([ val for val in data[v].values if not np.isnan(val) ] )
	    if max([ nvals[v] for v in nvals]) < min_nobs: return False
	    
	    # Yes!! 
	    return True
	Model = HMM_and_RF(min_nobs = model['min_nobs'], min_npats = model['min_npats'], nstates_all = model['nstates'], 
		covariance_type = model['covariance_type'],  n_estimators=model['n_estimators'], vitals=vitals_to_use)
	# get list of pans to use
	if not keep_ratio_half:
		dead_to_alive_ratio = float(len(dead_heart_pans))/float(len(alive_heart_pans))
	else:
		dead_to_alive_ratio = 1.0

	# Keep the alive/dead ratios roughly the same.
	nalive = int(float(Npatients)/(1. + dead_to_alive_ratio))
	ndead = Npatients - nalive

	# Build list of alive/dead pans
	alive_indices = np.arange(0,len(alive_heart_pans))
	dead_indices = np.arange(0,len(dead_heart_pans))

	np.random.shuffle(alive_indices)
	np.random.shuffle(dead_indices)

	idead = 0
	ialive = 0
	alive_pans = []
	dead_pans = []
	alive_pans_test = []
	dead_pans_test = []


	print "Building list of alive patients"
	while len(alive_pans) < nalive and ialive < len(alive_heart_pans):
	    pan = alive_heart_pans[alive_indices[ialive]]
	    if pan_is_ok(pan):
	        alive_pans.append(pan)
	    ialive+=1

	while len(alive_pans_test) < ntest_alive and ialive < len(alive_heart_pans):
	    pan = alive_heart_pans[alive_indices[ialive]]
	    if pan_is_ok(pan, testing=True):
	        alive_pans_test.append(pan)
	    ialive+=1

	print "Building list of dead patients"
	while len(dead_pans) < ndead and idead < len(dead_heart_pans):
	    pan = dead_heart_pans[dead_indices[idead]]
	    if pan_is_ok(pan):
	        dead_pans.append(pan)
	    idead+=1  

	while len(dead_pans_test) < ntest_dead and idead < len(dead_heart_pans):
	    pan = dead_heart_pans[dead_indices[idead]]
	    if pan_is_ok(pan, testing=True):
	        dead_pans_test.append(pan)
	    idead+=1  

	 

	alive_pans = np.array(alive_pans)
	dead_pans = np.array(dead_pans)


	#Model.predict_proba(alive_pans, dead)

	alive_pans_test = np.array(alive_pans)
	dead_pans_test = np.array(dead_pans)


	all_pans = []
	labels = []
	all_pans.extend(alive_pans)
	all_pans.extend(dead_pans)

	for a in all_pans:
	    if a in alive_pans: labels.append(0)
	    elif a in dead_pans: labels.append(1)
	all_pans = np.array(all_pans)
	labels = np.array(labels)

	all_test_pans = []
	test_labels = []
	all_test_pans.extend(alive_pans_test)
	all_test_pans.extend(dead_pans_test)
	for a in all_test_pans:
		if a in alive_pans_test: test_labels.append(0)
		else: test_labels.append(1)

	all_test_pans = np.array(all_test_pans)

	Model.fit(alive_pans, dead_pans)
	test_scores = Model.predict_proba(all_test_pans)[:,1]
	'''
	printer = zip(test_scores, test_labels)
	#print printer
	for p in printer:
		if p[1] == 1:
			result = p[0] > 0.5
		else:
			result = p[0] < 0.5
		print p[0], p[1], result 
	
	print fpr, tpr
	'''
	fpr, tpr, _ = roc_curve(test_labels, test_scores)

	f = plt.figure(figsize=(12,4),tight_layout=True, dpi=500)
	ax_roc = f.add_subplot(232)
	ax_dist = f.add_subplot(235)
	ax_importance = f.add_subplot(233)
	ax_pca = f.add_subplot(236)

	summary_text = ""
	translate = {
		'min_time_before_end_for_testing' : 'delta_t'
	}
	#skip_params = [ 'min_npats', 'ntest_dead', 'ntest_alive', 'keep_ratio_half'  ]
	skip_params = []
	for param in model:
		if param in skip_params: continue
		if param in translate: param_text = translate[param]
		else: param_text = param
		summary_text="%s%-20s %-10s\n"%(summary_text, param_text, str(model[param]))

	f.text(0.1, 0.9, "Parameters\n(%s)\n----------"%(model_name),fontsize=14,family='monospace',horizontalalignment='left', verticalalignment='top')
	f.text(0.1, 0.75, summary_text, fontsize=10, family='monospace',horizontalalignment='left', verticalalignment='top')
	ax_roc.set_xlabel("False Positive Rate")
	ax_roc.set_ylabel("True Positive Rate")
	ax_roc.plot(fpr, tpr, color='r')
	ax_roc.plot([0,1],[0,1],ls='--',color='k')
	ax_roc.set_title("ROC (area = %.2f)"%(auc(fpr, tpr)))

	ax_dist.set_xlabel("score")
	test_scores_dead = [ test_scores[i] for i in range(len(test_labels)) if test_labels[i] == 1 ]
	test_scores_alive = [ test_scores[i] for i in range(len(test_labels)) if test_labels[i] == 0 ]
	ax_dist.hist(test_scores_dead, 20, range=(0,1), color='k', alpha=0.3, label="Expired")
	ax_dist.hist(test_scores_alive, 20, range=(0,1), color='g', alpha=0.3, label="Alive")
	ax_dist.set_ylabel("Number of patients")
	
	#max_num_importances = 
	importances = Model.classifier.feature_importances_
	inds = np.arange(0, len(Model.vitals_available))
	#inds_sorted = np.argsort(importances)
	width=0.35
	yd = importances[::2]
	ya = importances[1::2]
	#print yd
	#print ya
	
	ax_importance.bar(inds, ya, width, color='g', alpha=0.3)
	ax_importance.bar(inds+width, yd, width, color='k', alpha=0.3)
	ax_importance.set_xticks(inds+width)
	ax_importance.set_xticklabels(())
	for i,v in enumerate(Model.vitals_available):
		vtext = v[4:]
		ax_importance.text(i + width, 0.01, vtext, horizontalalignment='center', verticalalignment='bottom', rotation='vertical',color='r',fontsize=9 )
	#ax_importance.set_xticklabels(tuple([ v for v in Model.vitals_available ]))
	#plt.setp( ax_importance.xaxis.get_majorticklabels(), rotation=70 )
	ax_importance.set_ylabel("Importance")
	#save_obj(alive_pans,'alive_pans.list')
	#save_obj(dead_pans, 'dead_pans.list')


	pca = PCA()
	riskvects = np.array(Model.risk_vectors_all)
	rv_alive = np.array(Model.risk_vectors_alive)
	rv_dead = np.array(Model.risk_vectors_dead)


	mus = np.array([  np.mean([ rv[i] for rv in riskvects  ]) for i in range(len(riskvects[0])) ])
	sigmas =  np.array([  np.std([ rv[i] for rv in riskvects  ]) for i in range(len(riskvects[0])) ])

	def zscore(rv):
		rv -= mus 
		for i in range(len(rv)):
			if sigmas[i] > 0: rv[i]/=sigmas[i]
		#print rv
		return rv
	def normalize(rvs):
		rv_normed = np.array([ zscore(rv) for rv in rvs ])
		return rv_normed

	
	riskvects = normalize(riskvects)
	for i,r in enumerate(riskvects):
		for j,p in enumerate(r):
			if np.isnan(p): print i,j
	#print riskvects[0]
	rv_alive = normalize(rv_alive)
	rv_dead = normalize(rv_dead)


	pca.fit(riskvects)

	pc_alive = pca.transform(rv_alive) 
	pc_dead = pca.transform(rv_dead)

	x_comp = 0
	y_comp = 1
	ax_pca.scatter([ pc[x_comp] for pc in pc_dead ],[ pc[y_comp] for pc in pc_dead ], c='k',alpha=0.3 )
	ax_pca.scatter([ pc[x_comp] for pc in pc_alive ],[ pc[y_comp] for pc in pc_alive ], c='g',alpha=0.3 )

	ax_pca.set_xlabel("PCA %d"%(x_comp))
	ax_pca.set_ylabel("PCA %d"%(y_comp))
	save_obj(Model,'model_%s.obj'%(model_name))

	f.savefig("model_%s_results_Npatients%d.png"%(model_name,model['Npatients']), dpi=300)
	plt.show()

	#pan1 = alive_pans[0]
	#pan2 = dead_pans[0]




#save_obj(risk_vectors_dead, "risk_vectors_dead.list")
#save_obj(risk_vectors_alive, "risk_vectors_alive.list")
#save_obj(vitals_available, "model_vitals_used.list")
#save_obj(HMMs_dead,"HMMs_dead.model")
#save_obj(HMMs_alive, "HMMs_alive.model")
#save_obj(HMMs_alive, "HMMs_alive.model")
#save_obj(HMMs_alive, ".model")



