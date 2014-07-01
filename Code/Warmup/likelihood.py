# Author: Jorge L. García
# Date: 7/1/2014

#Import packages
import numpy as np
from scipy.optimize import minimize 
#Set seed
np.random.seed(0)

#2.1 Basic optimization
#Define the Rosenbrock Function
def rbfun(x):
    return sum(100.0*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0)

#Initial Condition
x0 = np.array(10*[0.5])

#Optimze
opt = minimize(rbfun, x0, method='BFGS', options={'disp': True, 'maxiter': 50000})

#Print optimization information
print (opt.x)
print (opt.success)
print (opt.message)

#2.3 (Estimation I)
#Generate data
N = 10000
x = np.array([np.random.normal(0,1) for i in range(0,N)])
n = float(x.size)
mle_1 = (1/n)*sum(x**2)
print mle_1

#2.4 (Estimation II)
#Define the negative of the likelihood function (the routine minimizes)
def llk(theta):
    return .5*(log(2.0*pi)  + log(theta) + (1.0/(n*theta))*sum(x**2))

#Initial Condition
theta0 = 1.7
 
#Optimize
opt = minimize(llk, theta0, method='BFGS', options={'disp': True, 'maxiter': 50000})

#Print optimization information
mle_2 = opt.x
print mle_2
print (opt.success)
print (opt.message)