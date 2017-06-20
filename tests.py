import sys, os, math, re
import numpy as np
# import ExtendedTSC
# import MDAnalysis as md

import matplotlib.pyplot as plt

x = np.array([1.,2,3,4,5])
y = np.array([1.,2,3,4,5])
s = np.array([0.1,0.2,0.4,0.8,1.6])

f,ax = plt.subplots()
ax.plot(x,y,c='black',lw=3)
ax.fill_between(x,y+s,y-s,facecolor='gold')
plt.show()
