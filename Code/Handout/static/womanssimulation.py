# Author: Jorge L. García
# Date: 7/1/2014

#Import packages
import numpy as np
#Set seed
np.random.seed(0)

#Model Parameters
N = 1000; T = 6; betak = .5; betan = .2; sigmaeps = 0.4; pi = .2; 
gamma = .8; sigmaeta = 1; covepseta = .3; 

#Simimulation Parameters
y_lb = 0.0; y_up = 10.0;
k_lb = 0.0; k_up = 5.0;
z_lb = 0.0; z_up = 5.0;
n_lb = 0; n_up = 3

# Simulate observed variables
#y
y = np.random.uniform(low=y_lb,high=y_up,size=(N,T))

#z,k,n
z = np.zeros((N,T))
kappa = np.zeros((N,T))
n = np.zeros((N,T))

z_N = np.random.uniform(low=z_lb,high=z_up,size=N)
kappa_N = np.random.uniform(low=k_lb,high=k_up,size=N)
n_N = np.random.random_integers(low=n_lb,high=n_up,size=N)

for t in range(T):
    z[:,t] = z_N
    kappa[:,t] = kappa_N
    n[:,t] = n_N

#Simulate unobserved variables
epseta = np.random.multivariate_normal([0,0],[[sigmaeps**2, covepseta],
[covepseta, sigmaeta**2]],(N,T))
eps = epseta[:,:,0]
eta = epseta[:,:,1]

#Construct latent and observed decision
w  = gamma*z + eta
U1 = y + gamma*z + eta - pi*n
U0 = y + betak*kappa + betan*n + eps
v = U1 - U0
d = v >= 0

np.savetxt('womansdata', np.column_stack((y, z, kappa, n, d, w)), fmt= '%8.3f')