import matplotlib
import cPickle as pickle
from med_data import *
from pandas import *
from pylab import *
import numpy as np
import re,os
from matplotlib.ticker import MultipleLocator
from scipy.optimize import minimize, show_options
from datetime import datetime, date, time
# import raw data
dead_encounters = Encounters.iloc[np.where(Encounters['DEATH_FLAG'] == 1.0)]
alive_encounters = Encounters.iloc[np.where(Encounters['DEATH_FLAG'] == 0.0)]
ALIVE_PANS = alive_encounters.PAN.tolist()
DEAD_PANS = dead_encounters.PAN.tolist()

def get_random_PANS(npans):
    inds = np.arange(0,len(PANS))
    np.random.shuffle(inds)
    return PANS[inds[:npans]]



class Patient:
    def __init__(self, encounter_entry, vitals_to_use):
        self.info = {}
        for key in ENCOUNTER_COLUMNS: self.info[key] = encounter_entry[key]
        if float(self.info['DEATH_FLAG']) == 1.0: self.dead = True
        else: self.dead = False
        self.vitals = {}
        for vs in vitals_to_use:
            self.vitals[vs] = {'time' : [], 'measurement' : []}
    def add_vital(self, E):
        if E['VITAL_DESCRIPTION'] in self.vitals:
            vs = E['VITAL_DESCRIPTION']
            #print vs
            meas = float(E['VITAL_SIGN_VALUE'])
            tmeas = (E['VITAL_TAKEN_DATE'] - self.info['ADM_DATE'])/(60.0*60.0)
            #print E['VITAL_TAKEN_DATE'] , self.info['ADM_DATE'], tmeas
            #print E['VITAL_TAKEN_DATE'], self.info['ADM_DATE']
            #if vs == "Heart Rate (beats/min)":
            #    print "  Adding %s of %f at %s (t-t0 = %.2f hrs)"%(vs,meas,E['VITAL_TAKEN_DATE'],tmeas)
                
            self.vitals[vs]['measurement'].append(meas)
            self.vitals[vs]['time'].append(tmeas)
            seq = self.vitals[vs]['time']
            #print seq
            argsort_inds = sorted(range(len(seq)), key=seq.__getitem__)
            if not all(argsort_inds[:] == range(len(self.vitals[vs]['time']))):
                self.vitals[vs]['measurement'] = np.array(self.vitals[vs]['measurement'])[argsort_inds].tolist()
                self.vitals[vs]['time'] = np.array(self.vitals[vs]['time'])[argsort_inds].tolist()
    def add_vitals(self,vital_entries):
        for i in range(len(vital_entries)): self.add_vital(vital_entries.iloc[i])
    def add_labs(self,entries):
        self.lab_orders = {}
    def plot_vitals(self,vital_name,ax=None,fmt='b^-',label=None,title=None,ylabel=None,xlabel="Hours after admission"):
        if ax is None:
            f = plt.figure()
            ax = f.add_subplot(111)
        if not ylabel is None: ax.set_ylabel(ylabel)
        if not xlabel is None: ax.set_xlabel(xlabel)
        if not title is None: ax.set_title(title)
        if not label is None: ax.plot(self.vitals[vital_name]['time'], \
                                      self.vitals[vital_name]['measurement'],fmt,alpha=0.5,label=label)
        else: ax.plot(self.vitals[vital_name]['time'], self.vitals[vital_name]['measurement'],fmt,alpha=0.5)
        mind = 0.5
        xmax = ax.get_xlim()[1]
        xmin = ax.get_xlim()[0]
        ymax = ax.get_ylim()[1]
        ymin = ax.get_ylim()[0]
        
        D = (ymax - ymin)
        majloc = mind
        while D/majloc > 3: majloc+=mind
        minloc = majloc/5.
        
        ax.yaxis.set_major_locator(MultipleLocator( majloc ))
        ax.yaxis.set_minor_locator(MultipleLocator( minloc ))
        
        D = (xmax - xmin)
        majloc = mind
        while D/majloc > 3: majloc+=mind
        minloc = majloc/5.
        
        ax.xaxis.set_major_locator(MultipleLocator( majloc ))
        ax.xaxis.set_minor_locator(MultipleLocator( minloc ))
        
        
        return ax
    def summary_string(self):
        
        return "PAN: %d"%(self.info['PAN'])

def fname(pan):
    return "patients/patient_PAN%d.pkl"%(int(pan))

def save_patient(patient):
    print "Saving patient ", patient.info['PAN']
    with open(fname(patient.info['PAN']), 'wb') as output:
        pickle.dump(patient, output, pickle.HIGHEST_PROTOCOL)

def is_saved(pan):
    return os.path.exists(fname(pan))

def load_patient(pan):
    with open(fname(pan), 'rb') as inputf:
        return pickle.load(inputf)

def get_n_patients(n, vitals_to_use, 
  random=True, condition_function = None, verbose=True, pans=PANS,
  overwrite=False):
    patients = []
    
    if random: 
        inds = np.arange(0,len(pans))
        np.random.shuffle(inds)
        pans = pans[inds]
    
    for p in pans:
        if len(patients) == n: return patients
        if is_saved(p) and not overwrite:
            patients.append(load_patient(p))
            if verbose: print " Added patient ", len(patients)
            continue
        encounters = Encounters.iloc[np.where(Encounters['PAN'] == p)]
        try:
            vinds = np.where(Vitals['PAN'] == p)
            vitals = Vitals.iloc[vinds]
        except:
            continue
        if len(encounters) > 1: 
            if verbose: print "More than one encounter for PAN=",p, len(encounters)
            continue
        encounter = encounters.iloc[0]
        P = Patient(encounter, vitals_to_use)
        P.add_vitals(vitals)
        if condition_function is not None:
            if not callable(condition_function): 
                print "Condition function is not callable!"
                return None
            else: 
                if not condition_function(P, vitals_to_use): continue
        
        save_patient(P)
        patients.append(P)
        if verbose: print " Added patient ", len(patients)
    






        
