
import pystan
import matplotlib
import pyfits
import numpy as np
import pickle
import sys
from pystan import StanModel

import matplotlib.pyplot as plt

#parse parameters from command line
if len(sys.argv)!=6:
    print '\n Error! Incorrect Usage! \n'
    print 'Proper Usage: python MSE_by_long.py <nsplit> <nwind> <fitType> <plotBool> <mseBool> \n'
    print 'Example: python MSE_by_long.py 2000 10 MLE T T \n'
    sys.exit(2)
else:
    nsplit = int(sys.argv[1])   #number of windows to create
    nwind = int(sys.argv[2])    #number of windows to evaluate
    fitType = sys.argv[3]       #MLE or MCMC
    plotBool = sys.argv[4]      #whether to plot the fits
    mseBool = sys.argv[5]       #whether to compute SSE

#initialize QOIs
MSE = []
long = []

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
filename='allskyfuv_all.fits'
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
inc = np.max(gl)/nsplit

print '\nTotal Windows: ', nwind,'\n'
for i in range(0,nwind):
    if i%10==0:
        print 'Window: ', i, ' [',100*float(i)/float(nwind),'%]'
    #determine range of the data we want
    p=np.where((gl >= inc*i) & (gl < inc*(i+1)) & (fuv > 0) & (i100 > 0))

    #take appropriate subsets of the data
    glsub=gl[p]
    gbsub=gb[p]
    fuvsub=fuv[p]
    i100sub=i100[p]
    ndat=glsub.size

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
        fit = pystan.stan(model_code=dust_code, data=dust_dat,
                  iter=1000, chains=4)

        #print the fit
        print(fit)

        #extract simulated coefficients
        beta=fit.extract('beta')
        beta1=beta["beta"][:,0]
        beta2=beta["beta"][:,1]
        nsamp=beta1.size

    #plot fitting resutls, if asked to do so
    if plotBool=='T':
        ax=plt.axes()
        ax.scatter(i100sub[0:ndat],fuvsub[0:ndat],s=4,c='black',linewidth=0)
        ax.set_xlabel('I100')
        ax.set_ylabel('FUV')
        xax=np.arange(0,np.max(i100sub[0:ndat]))
        #plot single line or simluated lines dependeing on MLE v. MCMC fit type
        if fitType=='MLE':
            ax.plot(xax,beta1+beta2*xax,c='red',alpha=1)
        else:
            for j in np.arange(0,nsamp):
                ax.plot(xax,beta1[j]+beta2[j]*xax,c='black',alpha=0.009)
        plt.show()

    #compute MSE, if asked to do so
    if mseBool=='T':
        m_beta1 = np.mean(beta1)
        m_beta2 = np.mean(beta2)
        MSE.append(0)
        for j in np.arange(0,ndat):
            MSE[i] += (fuvsub[j]-m_beta1-m_beta2*i100sub[j])**2
        MSE[i]/=ndat
        long.append(inc*i)

#print results of SSE
ax=plt.axes()
ax.scatter(long,MSE,s=4,c='black',linewidth=2)
ax.set_xlabel('Longitude')
ax.set_ylabel('MSE')
xax=np.arange(0,np.max(i100sub[0:ndat]))
plt.show()



