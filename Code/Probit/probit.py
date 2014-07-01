# Author: Jorge L. García
# Date: 7/1/2014

#Import packages
import numpy as np
from scipy.stats import norm
from scipy.optimize import minimize 
#Set seed
np.random.seed(0)

#Exercise 1.3


#Generate data
N  = 1000
X1 = np.random.normal(0,1.5,N)
X2 = np.random.randint(0,1,N)
X  = zip(X1,X2)
E  = np.random.normal(0,1,N)

#Set arbitrary beta and construct latent
beta   = np.array([-.3,2])
Xbeta  = np.dot(X,beta)
y_star = Xbeta + E

#Construct latent and observed decision
d = y_star >=0


#Define the negative of the likelihood function (the routine minimizes)
def llk(B):
    XB = np.dot(X,B)
    l_0 = np.log(norm.cdf(-XB))
    l_1 = np.log(norm.cdf(XB))
    for i in range(0,N):
        if l_0[i] == -np.inf:
            l_0[i] = np.log(1e-3)
        elif l_1[i] == -np.inf:
            l_1[i] = np.log(1e-13)
    return -sum(d*l_1 + (1-d)*l_0)


# Solving the minimization#
#Initial Condition
beta0 = np.array([-.3,2.0])
 
#Optimize
opt = minimize(llk, beta0, method='BFGS', options={'disp': True, 'maxiter': 50000})

#Print optimization information
print (opt.x)
print (opt.success)
print (opt.message)
