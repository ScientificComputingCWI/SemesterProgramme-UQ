import numpy as np
from numpy.fft import fft,ifft
import matplotlib.pyplot as plt

def kernel(N):
    s = 100
    gaussian = lambda x, s: np.exp(-s*x**2.)
    K = gaussian(np.linspace(-1./2,1./2,N),s)
    K = N*np.fft.ifftshift(K)/np.sum(K)
    return K

def convolution(K,f):
    return ifft(fft(K)*fft(f))

def deconvolution(K,g,a):
    ftK = fft(K)
    ind = [i for i,v in enumerate(ftK) if abs(v) > a]
    L = np.zeros(len(ftK))
    L[ind] = 1./ftK[ind]
    L = ifft(L)
    return convolution(L,g)

def signal(N):
    x = np.linspace(0,1,N)
    b0 = 0
    b1 = int(np.ceil(N/9.))
    b2 = int(np.ceil(N*2/9.))
    b3 = int(np.ceil(N*3/9.))
    b4 = int(np.ceil(N*2/3.))
    b5 = N

    f = np.zeros(N)
    f[b0:b1] = x[b0:b1]*(1./3-x[b0:b1])*40
    f[b1:b2] = 0.3
    f[b2:b3] = x[b2:b3]*(1./3-x[b2:b3])*40
    f[b3:b4] = np.sin(3*10.5*np.pi*(x[b3:b4]-1./3))*x[b3:b4]
    f[b4:b5] = 6*(1-x[b4:b5])**2

    return f*3
    

if __name__=="__main__":
    N = 1000
    f = signal(N)
    K = kernel(N)
    delta = 1e-12

    g = convolution(K,f)
    gn = g+np.random.randn(len(g))*delta
    fs = deconvolution(K,g,1e-11)
    fsn = deconvolution(K,gn,1e-11)

    myf = open("signal.txt", "w")
    myf.write(' '.join([str(x) for x in f]))
    myf.close()
    
    plt.plot(f)
    plt.show()
        
    plt.plot(g)
    plt.show()

    plt.plot(fs)
    plt.plot(fsn)    
    plt.show()
