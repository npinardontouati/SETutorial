# Author: Jorge L. García
# Date: 7/1/2014

import matplotlib.pyplot as mplot
import numpy as np 

#Symmetric Intervals 

#Parameters
X = 100

#TeX Options
mplot.rc('text', usetex=True)
mplot.rc('font', family='serif')

#Plots
for i in range(2,12,2):
    mplot.plot([x for x in range(-X,X)],[((1/(1+np.exp(-x))) - .5)*2*i for x in range(-X,X)])

#Labels
mplot.xlabel(r'$x$')

#Legends
mplot.legend(('a = 2', 'a = 4','a = 6', 'a = 8', 'a = 10'), loc='upper left')

#Title
mplot.title(r"$f(x)$ ="
            r"$\left[ \frac{1}{1+\exp(-x)} - .5\right] \times 2 \times a$",
          fontsize=16, color='black')
#space
mplot.subplots_adjust(top=0.9)

#Save and close
mplot.savefig("symmetric_intervals.png")
mplot.close("all")