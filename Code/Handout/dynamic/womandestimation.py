# Author: Jorge L. García
# Date: 7/1/2014

#Import packages
import numpy as np
from scipy.stats import norm
from scipy.optimize import minimize 

#Set seed
np.random.seed(0)

#Parameters
T = 6; N = 1000

#Import data
data = np.genfromtxt('womanddata', dtype = 'float')
y = data[:,0:T]
z = data[:,T:2*T]
kappa = data[:,2*T:3*T]
n = data[:,3*T:4*T]
d = data[:,4*T:5*T]
w = data[:,5*T:6*T]
h = data[:,6*T:7*T]

#Maximization
#Initial Condition (the initial conditions are over unconstrained parameters)
theta0 = np.array([0.5, np.log(0.4), 0.2, 0.8, 0.9, 0.85, np.log(1.0),
-np.log((2*.4*1.0)/(.3+.4*1.0) - 1),0.2])
# order: betak,sigmaeps,betan,gamma1,gamma2,delta,sigmaeta,covepseeta,pi

#Define the negative of the likelihood function (the routine minimizes)
#Non-stochastic component of the utility function
def W1(z,h,n,theta):
    gamma1 = theta['gamma1']; gamma2 = theta['gamma2']; pi = theta['pi']
    return gamma1*z + gamma2*h - pi*n
    
def W0(kappa,n,theta):
    betak = theta['betak']; betan = theta['betan']
    return betak*kappa + betan*n

#xi_star @T 
def xi_starTm0(z,h,n,kappa,theta):
    return W1(z,h,n,theta) - W0(kappa,n,theta)

#Emax @T   
def EmaxTm0(z,h,n,kappa,theta):
    sigmaxi = np.sqrt(theta['sigmaeta']**2 + theta['sigmaeps']**2 - 2*theta['covepseta'])
    cdf1 = norm(0,1).cdf(xi_starTm0(z,h,n,kappa,theta)/sigmaxi)
    pdf  = norm(0,1).pdf(xi_starTm0(z,h,n,kappa,theta)/sigmaxi)
    cdf0 = norm(0,1).cdf(-xi_starTm0(z,h,n,kappa,theta)/sigmaxi)
    return W1(z,h,n,theta)*cdf1 + W0(kappa,n,theta)*cdf0 + sigmaxi*pdf

#xi_star @t = 1,...,T-1
for t in range(1,T):
    exec("def xi_starTm" + str(t) + "(z,h,n,kappa,theta):\
    delta = theta['delta'];\
    return W1(z,h,n,theta) - W0(kappa,n,theta)   \
    + delta*EmaxTm" + str(t-1) + "(z,h+1.0,n,kappa,theta)\
    - delta*EmaxTm" + str(t-1) + "(z,h,n,kappa,theta)")

#Emax @ t = 2,...,T-1
for t in range(1,T-1):
    exec("def EmaxTm" + str(t) + "(z,h,n,kappa,theta):\
    sigmaxi = np.sqrt(theta['sigmaeta']**2 + theta['sigmaeps']**2 - 2*theta['covepseta']);\
    delta = theta['delta'];\
    cdf1 = norm(0,1).cdf(xi_starTm" + str(t) + "(z,h,n,kappa,theta)/sigmaxi);\
    pdf = norm(0,1).pdf(-xi_starTm" + str(t) + "(z,h,n,kappa,theta)/sigmaxi);\
    cdf0 = norm(0,1).cdf(-xi_starTm" + str(t) + "(z,h,n,kappa,theta)/sigmaxi);\
    return (W1(z,h,n,theta) + delta*EmaxTm" + str(t-1) + "(z,h+1.0,n,kappa,theta))*cdf1 + \
    (W0(kappa,n,theta) + delta*EmaxTm" + str(t-1) + "(z,h,n,kappa,theta))*cdf0 + \
    sigmaxi*pdf")


#Define the negative of the likelihood function (the routine minimizes)
def logllk(thetae):

    #Estimands (declare and transfer to constrained parameters)
    theta = {}
    theta['betak'] = thetae[0]
    theta['sigmaeps'] = np.exp(thetae[1])
    theta['betan'] = thetae[2]
    theta['gamma1'] = thetae[3]
    theta['gamma2'] = thetae[4]
    theta['delta'] = thetae[5]
    theta['sigmaeta'] = np.exp(thetae[6])
    #theta['covepseta'] = thetae[7]
    theta['covepseta'] = (1/(1+np.exp(-thetae[7])) - .5)*2*np.exp(thetae[1])*np.exp(thetae[6])   
    theta['pi'] = thetae[8]
    

    #Locals
    theta['sigmaxi'] = np.sqrt(theta['sigmaeta']**2 + theta['sigmaeps']**2 - 2.0*theta['covepseta'])   
    theta['covxieta'] = theta['sigmaeta']**2 - theta['covepseta']
    xi_star = np.zeros((N,T))
    for t in range(0,T):
        exec("xi_star[:," + str(t) + "] = xi_starTm" + str(T-1-t) 
        + "(z[:," + str(t) + "], h[:," + str(t) + "]\
        , n[:," + str(t) + "], kappa[:," + str(t) + "],theta)")

    #l_0
    l_0 = norm(0,1).cdf(-xi_star/theta['sigmaxi'])
    
    #l_1   
    pdf_l1 = norm(0, 1).pdf((w -z*theta['gamma1'] -h*theta['gamma2'])/theta['sigmaeta'])   
    arg1_cdf_l1 = xi_star + (theta['covxieta']/(theta['sigmaeta']**2))*(w -z*theta['gamma1'] -h*theta['gamma2'])
    arg2_cdf_l1 = np.sqrt(theta['sigmaxi']**2 -(theta['covxieta']**2)/(theta['sigmaeta']**2))
    cdf_l1 = norm(0,1).cdf(arg1_cdf_l1/arg2_cdf_l1) 
    
    l_1 = (1/theta['sigmaeta'])*pdf_l1*cdf_l1 
    
    return -np.sum(d*np.log(l_1) +(1.0-d)*np.log(l_0))
    
    
#Optimize
opt = minimize(logllk, theta0, method='Nelder-Mead', options={'disp': True, 'maxiter': 100000})
thetaopt = opt.x
#(map back unconstrained estimates to constrained estimates)
thetaopt[1] = np.exp(thetaopt[1])
thetaopt[6] = np.exp(thetaopt[6])
thetaopt[7] = (1/(1+np.exp(-thetaopt[7])) - .5)*2*thetaopt[1]*thetaopt[4]

#Print optimization information
print (thetaopt)
print (opt.success)
print (opt.message)