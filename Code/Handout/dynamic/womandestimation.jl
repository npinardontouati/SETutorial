# Author: Jorge L. García
# Date: 7/1/2014

# Import packages
using Distributions
using Optim
using DataFrames

# Set a seed
srand(1)

# Read data
womandsample = Array(readtable("womandsample.csv"))

# Define variables
T = 6; N = 1000
y = womandsample[:,1:T]
z = womandsample[:,T+1:2*T]
kappa = womandsample[:,2*T+1:3*T]
n = womandsample[:,3*T+1:4*T]
d = womandsample[:,4*T+1:5*T]
w = womandsample[:,5*T+1:6*T]
h = womandsample[:,6*T+1:7*T]


# Define functions
# non-stochastic terms of utility
function W1(z,h,n,theta)
        gamma1 = theta[4]; gamma2 = theta[5]; pi = theta[9]
        return gamma1*z + gamma2*h - pi*n
end

function W0(kappa,n,theta)
        betak = theta[1]; betan = theta[3]
        return betak*kappa + betan*n
end

# xistar, t = T - 0
function xistarTm0(z,h,n,kappa,theta)
        return W1(z,h,n,theta) - W0(kappa,n,theta)
end

# Emax, t = T - 0
function EmaxTm0(z,h,n,kappa,theta)
        sigmaeps  = exp(theta[2]); sigmaeta  = exp(theta[7]);
        covepseta = (1/(1 + exp(-theta[8])) - .5)*2*exp(theta[2])*exp(theta[7]);
        sigmaxi  = sqrt(sigmaeta^2 + sigmaeps^2 - 2*covepseta);
        cdf1 = cdf(Normal(0,1), xistarTm0(z,h,n,kappa,theta)/sigmaxi)
        pdf0 = pdf(Normal(0,1), -xistarTm0(z,h,n,kappa,theta)/sigmaxi)
        cdf0 = cdf(Normal(0,1), -xistarTm0(z,h,n,kappa,theta)/sigmaxi)
        return W1(z,h,n,theta).*cdf1 + W0(kappa,n,theta).*cdf0 + sigmaxi*pdf0
end

# xistar, t = T - 1,..., 1
for t = 1:T-1
        eval(parse("function xistarTm" * string(t) * "(z,h,n,kappa,theta)"
        * "delta = theta[6];"
        * "W1(z,h,n,theta) - W0(kappa,n,theta)"
        * "+" * "delta*EmaxTm" * string(t-1) * "(z,h+1,n,kappa,theta)" * "-" * "delta*EmaxTm" * string(t-1) * "(z,h,n,kappa,theta)"
        * "end"))
end

# EmaxTmt,   t = T -1, ..., 2
for t = 1:T-1
        eval(parse("function EmaxTm" * string(t) * "(z,h,n,kappa,theta)"
        * "sigmaeps  = exp(theta[2]); sigmaeta  = exp(theta[7]);"
        * "covepseta = (1/(1 + exp(-theta[8])) - .5)*2*exp(theta[2])*exp(theta[7]);"
        * "sigmaxi  = sqrt(sigmaeta^2 + sigmaeps^2 - 2*covepseta);"
        * "delta = theta[6];"
        * "W1(z,h,n,theta) - W0(kappa,n,theta)"
        * "+" * "delta*EmaxTm" * string(t-1) * "(z,h+1,n,kappa,theta)" * ".*"
        * "cdf(Normal(0,1), xistarTm" * string(t) * "(z,h,n,kappa,theta)/sigmaxi)"
        * "+" * "delta*EmaxTm" * string(t-1) * "(z,h,n,kappa,theta)" * ".*"
        * "cdf(Normal(0,1), -xistarTm" * string(t) * "(z,h,n,kappa,theta)/sigmaxi)"
        * "+" * "sigmaxi*pdf(Normal(0,1), -xistarTm" * string(t) * "(z,h,n,kappa,theta)/sigmaxi)"
        * "end"))
end

#xistar vector function
xistarvf = [xistarTm0]
for t = 1:T-1
        xistarvf = [eval(parse("xistarTm" * string(t)));xistarvf]
end

# Likelihood
# function as function of data and parameters
function llkd(z,h,n,kappa,theta)
        # transform parameter to avoid evaluating outside function domain
        betak     = theta[1]
        sigmaeps  = exp(theta[2])
        betan     = theta[3]
        gamma1    = theta[4]
        gamma2    = theta[5]
        delta     = theta[6]
        sigmaeta  = exp(theta[7])
        pi        = theta[9]
        covepseta = (1/(1 + exp(-theta[8])) - .5)*2*exp(theta[2])*exp(theta[7])

        # define locals
        covxieta = sigmaeta^2 - covepseta
        sigmaxi  = sqrt(sigmaeta^2 + sigmaeps^2 - 2*covepseta)
        xistar = zeros(N,T)

        for t = 1:T
                xistar[:,t] = xistarvf[t](z[:,t], h[:,t], n[:,t], kappa[:,t],theta)
        end

        #l_0
        l_0 = logcdf(Normal(0,1), -xistar./sigmaxi)

        #l_1
        pdfl1     = pdf(Normal(0,1), (w - z*gamma1 - h*gamma2)/sigmaeta)
        arg1cdfl1 = xistar + (covxieta/(sigmaeta^2))*(w - z*gamma1 - h*gamma2)
        arg2cdfl1 = sqrt(sigmaxi^2 - (covxieta^2)/(sigmaeta^2))
        cdfl1     = cdf(Normal(0,1), arg1cdfl1/arg2cdfl1)

        l_1 = log((1/sigmaeta)*pdfl1.*cdfl1)

        return l = -sum(d.*l_1 + (1-d).*l_0)
end

# function wrapper
function wllkd(theta)
        return llkd(z,h,n,kappa,theta)
end

# Optimization
# initial condition
theta0 = [0.5, log(0.4), 0.2, 0.8, 0.9, 0.85, log(1.0),-log((2*.4*1.0)/(.3+.4*1.0) - 1),0.2]

# display optimization
optimal = optimize(wllkd, theta0, NelderMead())
mled = optimal.minimizer
mled[2] = exp(mled[2])
mled[7] = exp(mled[7])
mled[8] = (1/(1+exp(-mled[8])) - .5)*2*mled[2]*mled[7]
println(mled)
