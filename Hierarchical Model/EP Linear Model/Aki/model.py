
from __future__ import division
import numpy as np
from scipy import linalg as sc_linalg
import pickle
import matplotlib.pyplot as plt


# Store LAPACK positive definite inverse routine
DPOTRI = sc_linalg.get_lapack_funcs('potri')

# Smoothing factor
SMOOTH = 0.1
# Smoothing is now done for Qt by (and similary for r)
# Qt'_n = (SMOOTH-1) * Qt_n + SMOOTH * Qt_{n-1},
# where Qt_n is the precision matrix of the tilted distribution (for
# one site) estimated from the mcmc samples in the iteration n, and
# Qt'_n is the smoothed precision matrix given back to the master.
# TODO try including more iterations into the smoothing.

# Use seed = None for random seed
RAND = np.random.RandomState(seed=0)

def main():

    # ----------------
    #       Data
    # ----------------
    

    
    J = 6                                # number of groups
    Nj = RAND.randint(80,90,size=J)      # number of observations per group
    N = np.sum(Nj)                       # number of observations
    K = 6                                # number of inputs
    M = 2                                # number of workers

    # Evenly distributed group indices for M parallel jobs
    iiM = tuple(np.arange((J//M+1)*m, (J//M+1)*(m+1))
                if m < J%M else
                np.arange(J//M*m, J//M*(m+1)) + J%M
                for m in range(M))

    # observation index limits for J groups
    iiJ = np.concatenate(([0], np.cumsum(Nj)))
    
    # Model definition:
    # y_i ~ B(alpha_i + x_i * beta)
    # Local parameter alpha_i ~ N(0, sigma2_a)
    # Shared parameter phi = [log(sigma2_a), beta]
    # Hyperparameter sigma2_a
    # Fixed parameters beta
    
    # simulate from the model
    sigma2_a = 2**2
    a_sim = RAND.randn(J)*np.sqrt(sigma2_a)
    b_sim = RAND.randn(K)
    X = RAND.randn(N,K)
    y = X.dot(b_sim)
    for j in range(J):
        y[iiJ[j]:iiJ[j+1]] += a_sim[j]
    y = 1/(1+np.exp(-y))
    y = (RAND.rand(N) < y).astype(int)


    # -------------------------------------------------
    # Load  pre-built stan model for the tilted moments
    # -------------------------------------------------
    model_name = 'hier_log'
    with open(model_name+'.pkl', 'rb') as f:
        sm = pickle.load(f)
    
    
    # -----------------
    # Initialise the EP
    # -----------------
    
    nchains = 1  # Requires nchains cores
    nsamp = 800
    warmup = 100
    thin = 2
    # Using thin > 1 seems to cause segmentation fault on fit.extract().
    # Bug reported and fix is implemented into develop branch.
    
    niter = 20
    dphi = K+1

    # Priors for the shared parameters
    # sigma2_a = exp(phi[0]) ~ lognpdf(log(2),log(5))
    m0_a = np.log(2)
    V0_a = np.log(5)**2

    #t = np.linspace(1e-6,2,100)
    #plt.plot(t, sc.stats.lognorm.pdf(t, s=np.sqrt(V0_a), scale=np.exp(m0_a)))
    #plt.show()

    # prior for coefs
    m0_b = 0
    V0_b = 1**2

    # natural parameters of the prior
    Q0 = np.diag(np.append(1./V0_a, np.ones(K)/V0_b)).T
    r0 = np.append(m0_a/V0_a, np.ones(K)*(m0_b/V0_b))
    
    # natural parameters
    Q = Q0.copy(order='F')
    r = r0.copy(order='F')
    
    # mean cov of (Q,r)
    C_phi = np.empty((dphi,dphi), order='F')
    m_phi = np.empty(dphi)
    
    
    # natural site parameters
    Qi = np.zeros((dphi,dphi,J), order='F')
    ri = np.zeros((dphi,J), order='F')
    # natural site proposal parameters
    Qi2 = np.zeros((dphi,dphi,J), order='F')
    ri2 = np.zeros((dphi,J), order='F')
    
    # site parameter updates
    dQi = np.empty((dphi,dphi,J), order='F')
    dri = np.empty((dphi,J), order='F')
    
    # check positive definitness for each cavity distribution
    posdefs = np.empty(J, dtype=np.bool_)

    df0 = 0.8         # Damping factor for posterior approximation updates
    dfmin = 0.01      # Minimum damping factor
    
    workers = [worker(dphi,
                      X[iiJ[ji]:iiJ[ji+1],:],
                      y[iiJ[ji]:iiJ[ji+1]],
                      sm,
                      nchains,
                      nsamp,
                      warmup,
                      thin)
               for ji in range(J)]
    
    # Utility variable for copying the upper triangular to bottom
    upind = np.triu_indices(dphi,1)
    
    # Convergence analysis
    dm = np.zeros((niter,J))
    dV = np.zeros((niter,J))
    m_phi_s = np.empty((niter, dphi))
    err = np.empty((niter+1, dphi))
    # Plotting variables
    phi_t = np.append(np.log(sigma2_a), b_sim)
    inds = np.arange(1,K+2)
    plot_intermediate = False
    
    # -----------------------------------------------------------------
    # This part can be re-evaluated to continue from the previous stage
    # -----------------------------------------------------------------

    for i1 in range(niter):
        
        print 'Iter {} ...'.format(i1)
        
        df = df0
        while True:
            
            # -------------
            # Try EP update
            
            if i1 > 0:
                # skip this in the first round because required quantities
                # are computed later
                np.add(Qi, np.multiply(df, dQi, out=Qi2), out=Qi2)
                np.add(ri, np.multiply(df, dri, out=ri2), out=ri2)
                
                np.add(Qi2.sum(2, out=Q), Q0, out=Q)
                np.add(ri2.sum(1, out=r), r0, out=r)
            
            # check for positive definiteness
            np.copyto(C_phi, Q)
            try:
                sc_linalg.cho_factor(C_phi, overwrite_a=True)
            except sc_linalg.LinAlgError:
                # not positive definite -> reduce damping factor
                df *= 0.8
                print 'Neg def posterior cov,', \
                      'reducing df to {:.3}'.format(df)
                if i1 == 0:
                    print 'Invalid prior'
                    return False
                if df < dfmin:
                    print 'Damping factor reached minimum'
                    return False
                continue
            
            # -------------------------------
            # cavity distributions (parallel)
            # check positive definitness for each cavity distribution
            for m in range(M):
                # run jobs in parallel
                for j in range(len(iiM[m])):
                    ji = iiM[m][j] # group to update
                    
                    posdefs[ji] = workers[ji].cavity(Q, Qi2[:,:,ji],
                                                          r, ri2[:,ji])

            if np.all(posdefs):
                # All cavity distributions are positive definite.
                
                # Accept step (switch Qi-Qi2 and ri-ri2)
                temp = Qi
                Qi = Qi2
                Qi2 = temp
                temp = ri
                ri = ri2
                ri2 = temp
                
                # N.B. The following inversion could be done while parallel jobs
                # are running, thus saving time. Also this is only needed for
                # convergence analysis.
                # Invert Q
                _, info = DPOTRI(C_phi, overwrite_c=True)
                if info:
                    # Inversion failed
                    print "DPOTRI failed with error code {}".format(info)
                    return False
                # Copy the upper triangular into the bottom
                # This could be made more efficient ... cython?
                C_phi.T[upind] = C_phi[upind]
                # Calculate mean
                C_phi.dot(r, out=m_phi)
                
                break
                
            else:
                # Not all cavity distributions are positive definite
                # reduce the damping factor
                df *= 0.8
                ndef_ind = np.nonzero(~posdefs)[0]
                print 'Neg def cavity cov in site(s) {},'.format(ndef_ind), \
                      'reducing df to {:.3}'.format(df)
                if i1 == 0:
                    print 'Invalid prior'
                    return False
                if df < dfmin:
                    print 'Damping factor reached minimum'
                    return False
        
        # Check approximation        
        m_phi_s[i1] = m_phi
        np.subtract(m_phi, phi_t, out=err[i1])
        if plot_intermediate:
            cint = 1.96*np.sqrt(np.diag(C_phi))
            plt.plot(inds, phi_t, 'b',
                     inds, m_phi, 'r',
                     inds, m_phi+cint, 'r--',
                     inds, m_phi-cint, 'r--')
            plt.show()
        
        # -------------------------------
        # tilted distributions (parallel)
        for m in range(M):
            for j in range(len(iiM[m])):
                # compute tilted moments for each group
                ji = iiM[m][j] # group to update
                
                # print 'Iter {}, sample tilted moments for group {} of {}' \
                #       .format(i1, ji+1, J)
                
                if i1 != niter-1:
                    dm[i1,ji], dV[i1,ji] = \
                        workers[ji].tilted(dQi[:,:,ji], dri[:,ji], C_phi, m_phi)
                else:
                    # Final iteration ... memorise samples
                    workers[ji].tilted_final()
                
        print 'Iter {} done, max diff in tilted mean {:.4}, in cov {:.4}' \
              .format(i1, np.max(dm[i1]), np.max(np.sqrt(dV[i1])))
#        print dm[i1,:]
#        print np.sqrt(dV[i1,:])
        
    
    # Form final approximation by mixing samples from all the sites
    # N.B. calculated variable by variable for memory reasons
    prcs_p = [2.5, 25, 50, 75, 97.5]
    prcs = np.empty((dphi,len(prcs_p)))
    means = np.empty(dphi)
    for i in range(dphi):
        comb = np.concatenate([w.samp[:,i] for w in workers])
        means[i] = np.mean(comb)
        np.percentile(comb, prcs_p, overwrite_input=True, out=prcs[i])
    np.subtract(means, phi_t, out=err[-1])
    
    print '\nResults:'
    print (' '*7+'{:>10}'*2+'{:>9}%'*len(prcs_p)) \
          .format('real', 'mean', *prcs_p)
    for i in range(dphi):
        print ('{:7}'+'{:>10.3}'*2+'{:>10.3}'*len(prcs_p)) \
              .format('phi[{}]'.format(i+1), phi_t[i], means[i], *prcs[i])
    
#    print ('\nFinal error\n'+'{:>10.3}'*dphi).format(*err[-1])
    
#    plt.plot(np.arange(niter+1), err)
#    plt.hold(True)
#    plt.plot(np.arange(niter+1), np.zeros(niter+1), 'k')
#    plt.show()
    print '\nmax diff for each group tilted mean'
    for i in xrange(niter-1):
        print ('iter {:>2}: '+'{:>10.3}'*J).format(i, *dm[i])
    print '\nmax diff for each group tilted cov'
    for i in xrange(niter-1):
        print ('iter {:>2}: '+'{:>10.3}'*J).format(i, *np.sqrt(dV[i]))
    
    print '\neach parameters error to the real'
    for i in xrange(niter+1):
        print ('iter {:>2}: '+'{:>10.3}'*dphi).format(i, *err[i])

class worker(object):
    
    
    def __init__(self, dphi, X, y, sm, nchains, nsamp, warmup, thin):
        
        # Allocate space for calculations
        # N.B. these arrays are used for various variables
        self.Qc = np.empty((dphi,dphi), order='F')
        self.rc = np.empty(dphi)
        self.vtemp = np.empty(dphi)
        # The following are commented out because shared memory is used
        # self.Q = np.empty((dphi,dphi), order='F')
        # self.r = np.empty(dphi)
        
        # Data for stan model in method tilted
        self.data = dict(N=X.shape[0],
                         K=X.shape[1],
                         X=X,
                         y=y,
                         mu_cavity=self.rc,
                         Sigma_cavity=self.Qc.T)
                         # Qc transposed in order to get C-order
        
        # Stan model
        self.sm = sm
        self.nchains = nchains
        self.nsamp = nsamp
        self.warmup = warmup
        self.thin = thin
        
        # Utilities
        self.upind = np.triu_indices(dphi,1)
        self.dphi = dphi
        
        # Memorise previous tilted distributions
        self.iteration = 0
        self.prev_Qt = np.empty((dphi,dphi), order='F')
        self.prev_rt = np.empty(dphi)
        
        # Final tilted samples
        self.samp = None
    
    
    def cavity(self, Q, Qi, r, ri):
        
        # >>> Simulating parallel programs with separate memory
        # np.copyto(self.Q, Q)
        # np.copyto(self.Qc, Qi)
        # np.copyto(self.r, r)
        # np.copyto(self.rc, ri)
        # np.subtract(self.Q, self.Qc, out=self.Qc)
        # np.subtract(self.r, self.rc, out=self.vtemp)
        # <<< Simulating parallel programs with separate memory
        # >>> ... or using shared memory
        self.Q = Q
        self.r = r
        np.subtract(self.Q, Qi, out=self.Qc)
        np.subtract(self.r, ri, out=self.vtemp)
        # <<< ... or using shared memory
        
        
        
        # use mean-cov parameters because of Stan
        # it would be nice if Stan would support natural parameters
        # ... and upper triangulars only
        try:
            sc_linalg.cho_factor(self.Qc, overwrite_a=True)
        except sc_linalg.LinAlgError:
            # Not positive definite
            # print "Neg def cavity cov"
            return False
        
        # Positive definite
        # Invert self.Qc
        _, info = DPOTRI(self.Qc, overwrite_c=True)
        if info:
            # Inversion failed
            print "DPOTRI failed with error code {}".format(info)
            return False
        # Copy the upper triangular into the bottom
        # This could be made faster and more memory efficient ... cython?
        self.Qc.T[self.upind] = self.Qc[self.upind]
        # Calculate mean and store it in self.rc
        np.dot(self.Qc, self.vtemp, out=self.rc)
        
        return True
        
        
    def tilted(self, dQi, dri, C_phi, m_phi):
        
        # Reuse the memory
        Qt = self.Qc
        rt = self.rc
        mt = self.vtemp
        
        # Sample from the model
        with suppress_stdout():
            samp = self.sm.sampling(data=self.data, chains=self.nchains,
                                     iter=self.nsamp, warmup=self.warmup,
                                     thin=self.thin, seed=RAND)
        
        # >>> This would be more efficient if samp.extract would not copy data
        # samp = samp.extract(permuted=False)[:,:,:self.dphi]
        # # Stack the chains
        # samp = np.lib.stride_tricks.as_strided(
        #            samp,
        #            (samp.shape[0]*samp.shape[1], samp.shape[2]),
        #            (samp.strides[0], samp.strides[2])
        #        )
        # <<< This would be more efficient if samp.extract would not copy data
        samp = samp.extract(pars='phi')['phi']
        
        
        # Sample mean and covariance
        np.mean(samp, 0, out=mt)
        samp -= mt
        np.dot(samp.T, samp, out=Qt.T)
        # Qt /= samp.shape[0]
        # Qt /= (samp.shape[0] - 1)
        Qt /= (samp.shape[0] - self.dphi)
        
        dm = np.max(np.abs(mt-m_phi))
        dV = np.max(np.abs(Qt-C_phi))
        
        try:
            sc_linalg.cho_factor(Qt, overwrite_a=True)
        except sc_linalg.LinAlgError:
            # Not positive definite -> discard update
            print 'Neg def tilted covariance'
            return (dm, dV)
        
        # Positive definite
        # Invert Qt
        _, info = DPOTRI(Qt, overwrite_c=True)
        if info:
            # Inversion failed
            print "DPOTRI failed with error code {}".format(info)
            dQi.fill(0)
            dri.fill(0)
            return (dm, dV)
        # Copy the upper triangular into the bottom
        # This could be made faster and more memory efficient ... cython?
        Qt.T[self.upind] = Qt[self.upind]
        np.dot(Qt, mt, out=rt)
        
        if self.iteration < 2:
            # First two iterations ... no smoothing
            self.iteration += 1
            
            np.subtract(Qt, self.Q, out=dQi)
            np.subtract(rt, self.r, out=dri)
        
        else:
            # Smooth tilted distribution using the previous ones
            np.multiply(Qt, 1-SMOOTH, out=dQi)
            np.multiply(self.prev_Qt, SMOOTH, out=self.prev_Qt)
            np.add(dQi, self.prev_Qt, out=dQi)
            np.multiply(rt, 1-SMOOTH, out=dri)
            np.multiply(self.prev_rt, SMOOTH, out=self.prev_rt)
            np.add(dri, self.prev_rt, out=dri)
        
            np.subtract(dQi, self.Q, out=dQi)
            np.subtract(dri, self.r, out=dri)
        
        # Switch self.Qc - self.prev_Qt (... and r) so that the current
        # precision matrix is now stored for next iteration
        temp = self.Qc
        self.Qc = self.prev_Qt
        self.prev_Qt = temp
        temp = self.rc
        self.rc = self.prev_rt
        self.prev_rt = temp
        # And update the pointer in self.data
        self.data['mu_cavity'] = self.rc
        self.data['Sigma_cavity'] = self.Qc.T
        
        return (dm, dV)
    
    
    def tilted_final(self):
        # Sample from the model
        with suppress_stdout():
            samp = self.sm.sampling(data=self.data, chains=self.nchains,
                                     iter=self.nsamp, warmup=self.warmup,
                                     thin=self.thin, seed=RAND)
        self.samp = samp.extract(pars='phi')['phi']
        
        

# >>> Temp solution to suppres output from STAN model (remove when fixed)
import os
class suppress_stdout(object):
    '''
    A context manager for doing a "deep suppression" of stdout and stderr in 
    Python, i.e. will suppress all print, even if the print originates in a 
    compiled C/Fortran sub-function.
       This will not suppress raised exceptions, since exceptions are printed
    to stderr just before a script exits, and after the context manager has
    exited (at least, I think that is why it lets exceptions through).      

    '''
    def __init__(self):
        # Open a pair of null files
        self.null_fds =  [os.open(os.devnull,os.O_RDWR) for x in range(2)]
        # Save the actual stdout (1) and stderr (2) file descriptors.
        self.save_fds = (os.dup(1), os.dup(2))

    def __enter__(self):
        # Assign the null pointers to stdout and stderr.
        os.dup2(self.null_fds[0],1)
        os.dup2(self.null_fds[1],2)

    def __exit__(self, *_):
        # Re-assign the real stdout/stderr back to (1) and (2)
        os.dup2(self.save_fds[0],1)
        os.dup2(self.save_fds[1],2)
        # Close the null files
        os.close(self.null_fds[0])
        os.close(self.null_fds[1])
# <<< Temp solution to suppres output from STAN model (remove when fixed)


if __name__ == '__main__':
    main()




