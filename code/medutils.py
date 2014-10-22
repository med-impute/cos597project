import csv
import matplotlib
import datetime
import time
from pandas import *
from pylab import *
from numpy import *
import re
from matplotlib.ticker import MultipleLocator
from scipy.optimize import minimize, show_options
from datetime import datetime, date, time
# import raw data
COHORT = read_csv('../../../../../data/penn/large/COHORT.csv')
LABS = read_csv('../../../../../data/penn/large/LABS.csv')
DIAGNOSIS = read_csv('../../../../../data/penn/large/DIAGNOSIS.csv')


vital_signs_to_use = [ 
                      'BP Noninvasive Diastolic (mm Hg)',
                      'BP Noninvasive Systolic (mm Hg)',
                      'Respirations (breaths/min)',
                      'Temperature (degrees F)',
                      'Heart Rate (beats/min)'
                      ]
#time_fmt = '%m/%d/%Y %H.%M.%S %p'  # Pandas does not currently recognize AM/PM. Manually fixing this.
time_fmt = '%m/%d/%Y %H.%M.%S'
date_re_str = "([0-9]{2})/([0-9]{2})/([0-9]{4}).([0-9]{2})\.([0-9]{2})\.([0-9]{2}).([APM]{2})"
date_re = re.compile(date_re_str)

def convert_tstr(tstr_old):
    timing = date_re.search(tstr_old)

    month = int(timing.group(1))
    day = int(timing.group(2))
    year = int(timing.group(3))
    
    hours = int(timing.group(4))
    minutes = int(timing.group(5))
    seconds = int(timing.group(6))
    AMPM = timing.group(7)
    if AMPM == "PM" and hours < 12: hours += 12
    elif AMPM == "AM" and hours == 12: hours = 0
    new_time = "%02d/%02d/%d %02d.%02d.%02d"%(month,day,year,hours,minutes,seconds)
    return new_time

def time_subt(t1str, t2str):
    t1 = convert_tstr(t1str)
    t2 = convert_tstr(t2str)
    
    dts = to_datetime(t1,format=time_fmt) - to_datetime(t2,format=time_fmt)
    dt = dts.total_seconds()
    return float(dt)/3600. # dt in hours
    
npatients = 1
patients = []
class Patient:
    def __init__(self, entries):
        self.info = {
                     'VISIT_ID' : None,
                     'AGE' : None,
                     'GENDER_DESCRIPTION' : None,
                     'RACE_DESCRIPTION' : None,
                     'HEIGHT_INCHES' : None,
                     'WEIGHT_LBS' : None,
                     'ADMIT_DATE' : None,
                     'DISHCATGE_DATE' : None
                     }
        for key in self.info:
            self.info[key] = entries[0][key]
        #print "Patient %s"%(self.info['VISIT_ID'])             
        self.vitals = {}
        for vs in vital_signs_to_use:
            self.vitals[vs] = {'time' : [], 'measurement' : []}
        for E in entries:
            if E['VITAL_DESCRIPTION'] in self.vitals:
                
                vs = E['VITAL_DESCRIPTION']
                
                meas = float(E['VITAL_SIGN_VALUE'])
                tmeas = time_subt(E['VITAL_TAKEN_DATE'], E['ADMIT_DATE'])
                #if vs == "Heart Rate (beats/min)":
                #    print "  Adding %s of %f at %s (t-t0 = %.2f hrs)"%(vs,meas,E['VITAL_TAKEN_DATE'],tmeas)
                    
                self.vitals[vs]['measurement'].append(meas)
                self.vitals[vs]['time'].append(tmeas)
    def add_labs(self,entries):
	# This only enters the lab order name and the time of the order.
	# Results 
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
        ht_ft = int(self.info['HEIGHT_INCHES']/12.0)
        ht_in = self.info['HEIGHT_INCHES'] - 12*ht_ft
        return "ID: %d, AGE: %.1f, SEX: %s, HT/WT: %d'%d''/%dlbs"%(
                                                                       self.info['VISIT_ID'],
                                                                       self.info['AGE'],
                                                                       self.info['GENDER_DESCRIPTION'],
                                                                       ht_ft, ht_in,
                                                                       self.info['WEIGHT_LBS']
                                                                       )
        
