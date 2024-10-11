import numpy as np
import matplotlib.pyplot as plt
from sympy.ntheory import primefactors, isprime
from scipy.special import zeta
from joblib import Parallel, delayed

def generator(n):
    # For prime n, find the primitive root modulo n
    if not isprime(n):
        raise ValueError('n must be prime')
    factorlist = primefactors(n-1)
    g = 2; i = 1
    while i <= len(factorlist):
        if pow(g,int((n-1)/factorlist[i-1]),n) == 1:
            g = g+1; i = 0
        i = i+1
    return g

def fastcbc(s,n,Gammaratio,gamma):
    # Fast CBC construction with POD weights
    m = int((n-1)/2)
    
    # Rader permutation
    g = generator(n)
    perm = np.ones(m,dtype=int)
    for j in range(1,m):
        perm[j] = np.mod(perm[j-1]*g,n)
    perm = np.minimum(perm,n-perm)-1
    permflip = np.flip(perm)

    # Precompute the FFT of the first column (permuted indices)
    bernoulli = lambda x: x**2-x+1/6
    fftomega = np.fft.fft(bernoulli(np.mod((perm+1)*(perm[-1]+1)/n,1)))
    z = np.zeros(s)
    p = np.zeros((s,m))

    for d in range(s):
        pold = np.vstack((np.ones((1,m)),p))
        x = gamma[d]*Gammaratio @ pold[:-1,:]
        if d == 0:
            minind = 1
        else:
            tmp = np.real(np.fft.ifft(fftomega * np.fft.fft(x[permflip])))
            minind = perm[np.argmin(tmp)]+1
        z[d] = minind
        omega = bernoulli(np.mod(minind*np.arange(1,m+1)/n,1))
        for l in range(d+1):
            p[l,:] = pold[l+1,:] + omega * pold[l,:] * Gammaratio[l] * gamma[d]
    return z

if __name__ == '__main__':
    # Discretize the ODE
    h = .01 # mesh size
    x = np.arange(0,1+h,h) # mesh
    ncoord = len(x)
    G = np.tril(np.ones(ncoord))
    G = G - .5*np.eye(ncoord)
    G[:,0] = .5
    X,_ = np.meshgrid(1-x,1-x)
    G = h*X*G
    G = G[50,:] # pick the row corresponding to coordinate x = 0.5

    # Specify the parametric diffusion coefficient
    s = 100 # stochastic dimension
    decay = 2 # decay rate of stochastic fluctuations
    indices = np.arange(1,s+1)
    deterministic = np.sin(np.pi*np.outer(x,indices))/indices**decay # precompute deterministic part
    a = lambda y: 1 + deterministic @ y # diffusion coefficient

    # ODE solution
    u = lambda y: G @ (1/a(y))
    
    # Weight structure
    amin = 1-zeta(decay)/2
    b = np.arange(1,s+1,dtype=float)**(-decay)/amin
    delta = .05
    lam = 1/(2-2*delta)
    Gammaratio = np.arange(1,s+1)**(2/(1+lam))
    gamma = (b/np.sqrt(2*zeta(2*lam)/(2*np.pi**2)**lam))**(2/(1+lam))

    # Solve the expected value of the quantity of interest
    nlist = [17,31,67,127,257,503,1009,2003,4001,8009,16007,32003,64007]
    R = 16 # number of random shifts
    rms = [] # store the computed RMS errors
    with Parallel(n_jobs=-1) as parallel:
        for n in nlist:
            print('n = ' + str(n))
            z = fastcbc(s,n,Gammaratio,gamma) # find generating vector
            shift = np.random.uniform(0,1,s)
            results = []
            for r in range(R):
                shift = np.random.uniform(0,1,s) # random shift
                tmp = parallel(delayed(u)(np.mod(i*z/n+shift,1)-1/2) for i in range(n))
                results.append(np.mean(tmp)) # compute the QMC mean for each random shift
            qmcavg = np.mean(results) # average over the random shifts
            rmserror = np.linalg.norm(qmcavg-results)/np.sqrt(R*(R-1)) # RMS errors
            rms.append(rmserror) # save the RMS errors for each n
            
    # Least squares fit for the error
    A = np.ones((len(nlist),2))
    nscaled = [R*i for i in nlist]
    A[:,1] = np.log(nscaled)
    lsq = np.linalg.solve(A.T @ A, A.T @ np.log(rms))
    lsq[0] = np.exp(lsq[0])

    # Visualize the results
    fig, ax = plt.subplots()
    ax.loglog(nscaled,lsq[0]*nscaled**lsq[1],'--b',linewidth=2,label='slope: ' + str(lsq[1]))
    ax.loglog(nscaled,rms,'.r','QMC error')
    ax.set_title('QMC error ($s = ' + str(s) + '$)',fontsize=15)
    ax.set_xlabel('number of nodes $nR$',fontsize=13)
    ax.set_ylabel('R.M.S. error',fontsize=13)
    ax.legend()
    plt.show()