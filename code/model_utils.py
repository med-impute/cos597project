from math import *
from pandas import *
from pylab import *
import numpy as np
DEFAULT_DEATH_FLAG = 0 # What to transform the '' death flags into 


def smooth(ts, vitals, MAXT = None, MINT = None, DT = 0.5):
    if len(ts) < 2: return ts, vitals
    times = np.array(ts)
    if MINT is None: MINT = 0
    times -= MINT
    if not MAXT is None: patient_max_time = MAXT - MINT
    else: patient_max_time = max(times)

    NBINS = int(round(patient_max_time/DT))
    for t in times: 
        if int(t/DT) > NBINS:
            print "NBINS = ",NBINS," MAXT = ",MAXT," max(times) = ",max(times), " MINT = ",MINT, "min(times) = ",min(times)

    smtimes = np.array([ (i + 0.5)*DT for i in range(NBINS) ])
    smoothed_vitals = np.empty(NBINS)
    num_meas = np.zeros(NBINS)
    smoothed_vitals[:] = np.nan
    
    for t,v in zip(times,vitals):
        i = int(t/DT)
        if i >= len(smtimes): i-=1
        if num_meas[i] == 0:
            smoothed_vitals[i] = v
        else:
            smoothed_vitals[i] += v
        num_meas[i] += 1
    for i in range(NBINS):
        if num_meas[i] > 0:
            smoothed_vitals[i]/=num_meas[i]
    return smtimes, smoothed_vitals

def indices_of_missing_values(vitals):
    inds = []
    for i,v in enumerate(vitals):
        if isnan(v): inds.append(i)
    return inds

def remove_nans(times, vitals):
    ts = []
    vs = []
    for i in range(len(vitals)):
        if not isnan(vitals[i]): 
            vs.append(vitals[i])
            ts.append(times[i])
    return np.array(ts), np.array(vs)

def impute_linear(times,vitals):
    # very simple imputation method -- basically linear interpolation
    # or, if the rest of the values are missing, fill them in with the last value
    ivitals = np.zeros(len(vitals))
    ivitals[:] = vitals
    
    for i in range(len(ivitals)):
        if isnan(ivitals[i]):
            # Test if the rest of the values are NaN's
            rest_are_nans = True
            for j in range(i, len(ivitals)):
                if not isnan(ivitals[j]):
                    rest_are_nans = False
                    J = j
                    break
                    
            if i == 0 and rest_are_nans:
                # All values are NaN's; nothing to do here...
                return None
            elif i == 0:
                for j in range(J):
                    ivitals[j] = ivitals[J]
            elif rest_are_nans:
                for j in range(i,len(ivitals)):
                    ivitals[j] = ivitals[i-1]
            else:
                y0 = ivitals[i-1]
                
                dy = ivitals[J] - ivitals[i-1]
                dt = times[J] - times[i-1]
                for j in range(i,J):
                    T = times[j] - times[i-1]
                    ivitals[j] = y0 + T*(dy/dt)
    return ivitals


def cut_early_times(patients, vitals_to_use, delta_t):
    for p in patients:
        MINT = max([ max(p.vitals[v]['time']) for v in vitals_to_use]) - delta_t
        for v in vitals_to_use:
            times = []
            vitals = []
            for T, V in zip(p.vitals[v]['time'], p.vitals[v]['measurement']):
                if T > MINT:
                    times.append(T)
                    vitals.append(V)
            p.vitals[v]['time'], p.vitals[v]['measurement'] = times, vitals
    

def smooth_vitals(patients, vitals_to_use, dt):
    for p in patients:
        MINT = min([ min(p.vitals[v]['time']) for v in vitals_to_use])
        MAXT = max([ max(p.vitals[v]['time']) for v in vitals_to_use])
        for v in vitals_to_use:
            p.vitals[v]['time'], p.vitals[v]['measurement'] = smooth(p.vitals[v]['time'], p.vitals[v]['measurement'], 
                                                                     MAXT = MAXT, MINT = MINT, DT = dt)
    
    
def impute_missing_values(patients, vitals_to_use, impute=impute_linear):
    for p in patients:
        for v in vitals_to_use:
            p.vitals[v]['measurement'] = impute(p.vitals[v]['time'], p.vitals[v]['measurement'])

def check_all_vitals_are_same_size(patient, vitals_to_use):
    i = None
    for v in vitals_to_use:
        if i is None: i = len(patient.vitals[v]['measurement'])
        else:
            if len(patient.vitals[v]['measurement']) != i: return False
    return True
       
def death_class(patient):
    if patient.dead: return 1
    else: return 0
    
def observations_and_classes(patients, vitals_to_use):
    observations = []
    classes = []
    for p in patients:
        observation = []
        mult_vitals_OK = check_all_vitals_are_same_size(p,vitals_to_use)
        if not mult_vitals_OK and len(vitals_to_use) > 1: 
            print "ERROR -- you can't use multiple vitals that have different dimensions. Impute them!"
            return None
        elif not mult_vitals_OK:
            observation =[ [vs] for vs in p.vitals[vitals_to_use[0]]['measurement'] ]
        else:
            for t in range(len(p.vitals[vitals_to_use[0]]['measurement'])):
                feature = []
                for v in vitals_to_use:
                    feature.append(p.vitals[v]['measurement'][t])
                observation.append(feature)
        observations.append(observation)
        classes.append(death_class(p))
    return observations, classes

