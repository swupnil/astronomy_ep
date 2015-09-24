
import pystan
import matplotlib
import pyfits
import numpy as np
import pickle
import sys
from pystan import StanModel

import matplotlib.pyplot as plt

#parse parameters from command line
if len(sys.argv)!=5:
    print '\n Error! Incorrect Usage! \n'
    print 'Proper Usage: python fit_by_long.py <gl_split> <gb_slit> <gl_wind> <fitType> \n'
    print 'Example: python fit_by_long.py 100 10 100 MLE \n'
    sys.exit(2)
else:
    gl_split = int(sys.argv[1])   #number of longitude windows to create
    gb_wind = int(sys.argv[2])   #number of latitude bins to create
    gl_wind = int(sys.argv[3])    #number of longitude windows to evaluate
    fitType = sys.argv[4]         #MLE or MCMC

#model code for MCMC
dust_code = """
    data {
    int<lower=0> N;
    vector[N] FUV;
    vector[N] I100;
    }
    parameters {
    vector[2] beta;
    real<lower=0> sigma;
    }
    model {
    FUV ~ normal(beta[1] + beta[2] * I100, sigma);
    }
    """

#read in the data
filename='../allskyfuv_all.fits'
print "Reading in:",filename
hdulist=pyfits.open(filename)
print "Done reading in:",filename
tbdata=hdulist[1].data
cols=hdulist[1].columns
print cols.names

#parse the variables
fuv=tbdata.field('FUV')
i100=tbdata.field('IR')
gl=tbdata.field('LONGITUDE')
gb=tbdata.field('LATITUDE')

#determine increment size of each window
gl_inc = float(360)/float(gl_split)
gb_inc = float(180)/float(gb_wind)

#initialize QOIs
Error = np.empty([gb_wind,gl_wind],float)
m_beta1 = []
m_beta2 = []
longbin = []

print '\nTotal Windows: ', gb_wind,'by', gl_wind, '\n'
for j in range(gl_wind):
    #print status
    if j%(gl_wind/20)==0:
        print 'Window: [',j, '] [',100*j/gl_wind,'%]'
        
    #determine range of the data we want
    p=np.where((gl>=gl_inc*j) & (gl<gl_inc*(j+1)) & (fuv>0) & (i100>0))

    #take appropriate subsets of the data
    glsub=gl[p]
    gbsub=gb[p]
    fuvsub=fuv[p]
    i100sub=i100[p]
    ndat=glsub.size
        
    #move to the next window if there are no data points
    if len(i100sub)==0 | len(fuvsub)==0:
        continue

    #calculate beta1 and beta2 depending on whether fitType is MLE or MCMC
    if fitType=='MLE':
        regression = np.polyfit(i100sub,fuvsub,1)
        beta1 = regression[1]
        beta2 = regression[0]
    else:
        #write the model
        dust_dat = {'N': ndat,
                'FUV': fuvsub[0:ndat],
                'I100': i100sub[0:ndat]}
        
        #run the model
        if j==1:
            fit = pystan.stan(model_code=dust_code, data=dust_dat,iter=1000, chains=4)
        else:
            fit = pystan.stan(fit=fit, data=dust_dat,iter=1000, chains=4)

        #print the fit
        print(fit)

        #extract simulated coefficients
        beta=fit.extract('beta')
        beta1=beta["beta"][:,0]
        beta2=beta["beta"][:,1]
        nsamp=beta1.size

    #computer errors of each point
    m_beta1.append(np.mean(beta1))
    m_beta2.append(np.mean(beta2))
    longbin.append(gl_inc*j)
    temp_error = []
    for k in np.arange(0,ndat):
        temp_error.append(np.abs(fuvsub[k]-m_beta1[j]-m_beta2[j]*i100sub[k]))

    #take average error in a bin
    for i in range(gb_wind):
        p_i = np.where((gbsub<=90-gb_inc*i) & (gbsub>90-gb_inc*(i+1)))
        if len(p_i[0])==0:
            Error[i][j] = -1000
        else:
            sum = 0
            for k in range(len(p_i[0])):
                sum += temp_error[p_i[0][k]]
            Error[i][j] = sum/len(p_i[0])

print 'All Computations Finished!'

print 'Plotting Solid Heatmap Now...'
#make heatmap
plt.xticks(np.arange(0,390,30))
plt.yticks(np.arange(-90,105,15))
plt.imshow(Error,extent=[0,360,-90,90], interpolation='nearest', origin='upper')
plt.show()

'''
print 'Plotting Slope by Longitude Now...'
#plot slope by longitude
ax = plt.axes()
ax.axis([-10,370,np.min(m_beta2)-5,np.max(m_beta2)+5])
ax.scatter(longbin,m_beta2,s=4,c='black',linewidth=2)
ax.set_xlabel('Longitude')
ax.set_ylabel('Slope')
plt.show()

print 'Plotting Intercept by Longitude Now...'
#plot intercept by longitude
ax = plt.axes()
ax.axis([-10,370,np.min(m_beta1)-5,np.max(m_beta1)+5])
ax.scatter(longbin,m_beta1,s=4,c='black',linewidth=2)
ax.set_xlabel('Longitude')
ax.set_ylabel('Intercept')
plt.show()
'''


