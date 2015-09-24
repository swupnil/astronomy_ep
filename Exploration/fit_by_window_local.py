
import pystan
import matplotlib
import pyfits
import numpy as np
import pickle
import sys
from pystan import StanModel

import matplotlib.pyplot as plt

#parse parameters from command line
if len(sys.argv)!=7:
    print '\n Error! Incorrect Usage! \n'
    print 'Proper Usage: python fit_by_window_local.py <gl_split> <gb_split> <gl_wind> <gb_wind> <fitType> <plotBool> \n'
    print 'Example: python fit_by_window_local.py 2000 2000 10 20 MLE T \n'
    sys.exit(2)
else:
    gl_split = int(sys.argv[1])   #number of longitude windows to create
    gb_split = int(sys.argv[2])   #number of latitude windows to create
    gl_wind = int(sys.argv[3])    #number of longitude windows to evaluate
    gb_wind = int(sys.argv[4])    #number of latitude windows to evaluate
    fitType = sys.argv[5]         #MLE or MCMC
    plotBool = sys.argv[6]        #whether to plot the fits

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
gb_inc = float(90)/float(gb_split)

#initialize QOIs
MSE = np.empty([gb_split,gl_split],float)
fig, ax = plt.subplots(6,10)

print '\nTotal Windows: ', gl_wind,'by', gb_wind, '\n'
for i in range(gb_wind):
    for j in range(gl_wind):
        #create overall index
        index = i*gl_split+j+1
        if (j%10==0):
            print 'Window: [', i,',',j, '] [',100*float(index)/float(gb_split)/float(gl_split),'%]'
        
        #determine range of the data we want
        p=np.where((gl>=gl_inc*j) & (gl<gl_inc*(j+1)) & (gb>=90-gb_inc*(i+1)) & (gb<90-gb_inc*i) & (fuv>0) & (i100>0))

        #take appropriate subsets of the data
        glsub=gl[p]
        gbsub=gb[p]
        fuvsub=fuv[p]
        i100sub=i100[p]
        ndat=glsub.size
        
        #move to the next window if there are no data points
        if len(i100sub)==0 | len(fuvsub)==0:
            MSE[i][j] = np.log(.01)
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
            if j==0 & i==0:
                fit = pystan.stan(model_code=dust_code, data=dust_dat, iter=1000, chains=4)
            else:
                fit = pystan.stan(fit=fit, data=dust_dat, iter=1000, chains=4)

            #print the fit
            print(fit)

            #extract simulated coefficients
            beta=fit.extract('beta')
            beta1=beta["beta"][:,0]
            beta2=beta["beta"][:,1]
            nsamp=beta1.size

        #plot slopes, if asked to do so
        if plotBool=='T':
            xax=np.arange(0,5)
            #plot single line or simluated lines dependeing on MLE v. MCMC fit type
            if fitType=='MLE':
                ax[i,j].plot(xax,beta2*xax,c='red',alpha=10)
                ax[i,j].set_ylim([-50,1600])
            else:
                for k in np.arange(0,nsamp):
                    ax[i,j].plot(xax,beta2[k]*xax,c='black',alpha=0.009)
            ax[i,j].set_xlabel('')
            ax[i,j].set_ylabel('')
            ax[i,j].set_xticklabels([])
            ax[i,j].set_yticklabels([])
            ax[i,j].text(0.5, 800, round(np.mean(beta2),2), fontsize = 13, color = 'blue')

        m_beta1 = np.mean(beta1)
        m_beta2 = np.mean(beta2)
        for k in np.arange(0,ndat):
            MSE[i][j] += (fuvsub[k]-m_beta1-m_beta2*i100sub[k])**2
        MSE[i][j]/=ndat
        MSE[i][j] = np.log(MSE[i][j])

#show plots of all regressions
if plotBool=='T':
    plt.show()

#make heatmap
plt.xticks(np.arange(0,390,30))
plt.yticks(np.arange(-90,105,15))
plt.imshow(MSE,extent=[0,360,-90,90],origin='upper')
plt.show()



