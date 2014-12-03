import numpy as np
from sklearn.gaussian_process import GaussianProcess
from matplotlib import pyplot as pl

patient_data = np.genfromtxt('/memex/gdarnell/patient_000001.csv',delimiter=',') # import data, set missing values to nan
HR = patient_data[:,3] # heart rate vector
HR = filter(lambda obs: not np.isnan(obs), HR) # remove all nan entries
x = np.atleast_2d(np.linspace(1,len(HR),len(HR))).T # create x-vector of same length as HR

our_kernel = lambda x, y: 0.5 * np.exp( -0.99 * 0.15 * np.sum((x - y)**2)) + np.pi # custom kernel
gp = GaussianProcess(corr=our_kernel,theta0=50,nugget = 1,regr='linear')
gp.fit(x,HR) # fit GP using custom kernel and parameters
y_pred, MSE = gp.predict(x, eval_MSE=True) # predict given x points
sigma = np.sqrt(MSE)

fig = pl.figure()
pl.plot(x,HR,'r*')
pl.plot(x, y_pred, 'b-', label=u'Prediction') # plot predicted values
pl.fill(np.concatenate([x, x[::-1]]),
        np.concatenate([y_pred - 1.9600 * sigma,
                       (y_pred + 1.9600 * sigma)[::-1]]),
        alpha=.5, fc='b', ec='None', label='95% confidence interval') # fill in confidence interval
gp.score(x,HR) # output overall score based on correlation
