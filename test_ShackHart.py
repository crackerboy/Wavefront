# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 23:03:30 2018

@author: rfetick
"""

from numpy import sqrt, log, array
import matplotlib.pyplot as plt

#z=Zernike(200,J=[6,7],unit=95)
#z.coeffs=array([1000,500])
z=Zernike(200,J=4,unit=95)
z.coeffs=array([5e3])
plt.figure()
plt.subplot(131)
z.look()

plt.subplot(132)
sh = ShackHart(200,20,RON=10,photonNoise=True)
sh.wavefront = 4.*z.toWF()
sh.look(sqrt)

plt.subplot(133)
plt.pcolormesh((sh.estimateWF()*sh.getValidMask()).T)
plt.title("Reconstructed")
plt.axis('equal')