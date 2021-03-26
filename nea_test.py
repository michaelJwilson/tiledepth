#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

x=np.linspace(-3,3,7)

# PSF value in each of the 7 pixels
# here a simple Gaussian
sigma=1.
P=np.exp(-x**2/2/sigma**2)
P/=np.sum(P)

flux=np.linspace(0,100,10)
nea=np.zeros(flux.size)

# typical readnoise variance
var_readnoise = 3**2

for i,f in enumerate(flux):
    # variance in pixel is readnoise + Poisson noise
    var_p = var_readnoise + f * P

    # variance of flux (inverse of Fisher)
    var_f = 1./(np.sum( 1/var_p * P**2))

    # because we define nea such that var_f = f + nea*var_readnoise
    nea[i] = (var_f - f)/var_readnoise

plt.plot(flux,nea,"-")
plt.xlabel("flux (electrons)")
plt.ylabel(r"$\alpha \cdot$ npix (NEA)")
plt.grid()
plt.show()
