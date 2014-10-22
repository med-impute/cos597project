import sys
from sklearn import hmm
from med_data import *
from model_list import *
from medutils_new import *
from model_utils import *
import sys
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 
overwrite = False
random_patient = False
verbose = False
Npatients = 2000
nexpired = Npatients/2
nalive = Npatients - nexpired


all_vitals = []
for m in models:
    for v in models[m]['vital_signs_to_use']:
        if not v in all_vitals: all_vitals.append(v)


def check_patient(patient, vitals_to_use):
    for v in vitals_to_use:
        if not v in patient.vitals:
            if verbose: print "%s isn't in the patient's list of vital signs..."%(v)
            return False
        if len(patient.vitals[v]['time']) < 2: 
            if verbose: print "Patient's %s signs are too few"%(v)
            return False
    dfs = [ '0.0', '1.0' ]
    DF = patient.info['DEATH_FLAG']
    
    if isnan(float( DF )): 
        if verbose: print "Death flag = '%s'"%(DF)
        return False
    return True

def isalive(patient):
	return death_class(patient) == 0

if overwrite:
    print "Getting alive patients"
    alive_patients = get_n_patients(nalive, all_vitals, 
                                    condition_function=check_patient, verbose=verbose, pans=np.array(ALIVE_PANS),
                                    random=random_patient,overwrite=overwrite)
    print "Getting dead patients"
    expired_patients = get_n_patients(nexpired, all_vitals, 
                                      condition_function=check_patient, verbose=verbose, pans=np.array(DEAD_PANS),
                                      random=random_patient,overwrite=overwrite)


for modelname in models:
    print "MODEL: ",modelname
    model = models[modelname]
    vital_signs_to_use = model['vital_signs_to_use']
    NSTATES = model['NSTATES']
    NMIXTURES = model['NMIXTURES']
    MODEL = model['MODEL']
    DT = model['DT']
    max_time = model['max_time']
    min_time = model['min_time']
    do_smoothing = model['do_smoothing']
    do_imputation = model['do_imputation']
    
    print "Getting alive patients"
    alive_patients = get_n_patients(nalive, vital_signs_to_use, 
                                    condition_function=check_patient, verbose=verbose, pans=np.array(ALIVE_PANS),
                                    random=random_patient)
    print "Getting dead patients"
    expired_patients = get_n_patients(nexpired, vital_signs_to_use, 
                                      condition_function=check_patient, verbose=verbose, pans=np.array(DEAD_PANS),
                                      random=random_patient)
    
    
    all_obs = []
    all_classes = []
    
    for selected_patients in [ alive_patients, expired_patients ]:
        if not min_time is None:
            pruned_patients = []
            for p in selected_patients:
                if check_patient(p,vital_signs_to_use): pruned_patients.append(p)
            print len(selected_patients) - len(pruned_patients)," patients pruned."
            selected_patients = [ p for p in pruned_patients ]
            cut_early_times(selected_patients, vital_signs_to_use, min_time)
            pruned_patients = []
            for p in selected_patients:
                if check_patient(p,vital_signs_to_use): pruned_patients.append(p)
            print len(selected_patients) - len(pruned_patients)," patients pruned."
            selected_patients = [ p for p in pruned_patients ]
            
        if len(selected_patients) == 0: continue
        pruned_again_idunno = []
        for s in selected_patients:
            bad = False
            for v in vital_signs_to_use:
                if v not in s.vitals:
                    bad = True
            if not bad: pruned_again_idunno.append(s)
        selected_patients = [ p for p in pruned_again_idunno ] 
        if do_smoothing:
            print "Smoothing"
            smooth_vitals(selected_patients, vital_signs_to_use, dt= DT)
    
        if do_imputation:
            print "Imputing"
            impute_missing_values(selected_patients, vital_signs_to_use)
        print "Getting obs/classes"
        observations, classes = observations_and_classes(selected_patients, vital_signs_to_use)
        all_obs.extend(observations)
        all_classes.extend(classes)
    all_obs = [  np.array(a) for a in all_obs ]
    if len(all_obs) == 0: continue
    should_continue = False
    for a in all_obs:
        if len(a) == 0: 
            should_continue = True
            break
    if should_continue: continue
    try:
        MODEL.fit(all_obs)
    except:
        continue
    statenames = np.arange(0,NSTATES)
    
    all_dstf = [ ]
    all_astf = [ ]
    for i in range(NSTATES):
        all_dstf.append([])
        all_astf.append([])
    
    death_state_freqs = np.zeros(NSTATES)
    alive_state_freqs = np.zeros(NSTATES)
    
    death_state_freqs_std = np.zeros(NSTATES)
    alive_state_freqs_std = np.zeros(NSTATES)
    
    for o,c in zip(all_obs,all_classes):
        states = MODEL.predict(o)
        for i in range(0,NSTATES):
            frq = float(len([ s for s in states if s == statenames[i] ]))/float(len(states))
            if c == 1:
                all_dstf[i].append(frq)
            else:
                all_astf[i].append(frq)
    
    for i in range(NSTATES):
        death_state_freqs[i] = np.mean(all_dstf[i])
        death_state_freqs_std[i] = sqrt(np.mean(np.power(np.array(all_dstf[i]) - death_state_freqs[i],2)))
    
        alive_state_freqs[i] = np.mean(all_astf[i])
        alive_state_freqs_std[i] = sqrt(np.mean(np.power(np.array(all_astf[i]) - alive_state_freqs[i],2)))
    figname_eps = "%s_frqs.eps"%(modelname)
    figname_png = "%s_frqs.png"%(modelname)
    width = 0.35
    f = plt.figure()
    ax = f.add_subplot(111)
    inds = np.arange(0,NSTATES)
    alive_bar = ax.bar(inds - width, alive_state_freqs, width, color='g', alpha=0.5)
    dead_bar = ax.bar(inds, death_state_freqs, width, color='r',alpha=0.5)
    ax.set_xlabel("State")
    ax.set_ylabel("Frequency")
    if min_time is None:
        ax.set_title("All time for %d patients"%(Npatients))
    else:
        ax.set_title("%.1f hrs before end for %d patients"%(min_time,Npatients))
    ax.legend((alive_bar[0],dead_bar[0]), ('Alive', 'Expired'))
    ax.set_ylim(0,1)
    f.savefig(figname_eps)
    f.savefig(figname_png)
    plt.show()

