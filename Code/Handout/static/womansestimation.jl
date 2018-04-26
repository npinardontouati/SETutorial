# Author: Jorge L. García
# Date: 7/1/2014

# Import packages
using Distributions
using Optim
using DataFrames

# Set a seed
srand(1)

# Read data
# womanssample = matrix(readtable("womanssample.csv"))

# Define variables
T = 6
y = womanssample[:,1:T]
z = womanssample[:,T+1:2*T]
kappa = womanssample[:,2*T+1:3*T]
n = womanssample[:,3*T+1:4*T]
d = womanssample[:,4*T+1:5*T]
w = womanssample[:,5*T+1:6*T]

# Likelihood
function llks(theta)
        # transform parameter to avoid evaluating
        # outside function domain
        betak     = theta[1]
        sigmaeps  = exp(theta[2])
        pibetan   = theta[3]
        gamma     = theta[4]
        sigmaeta  = exp(theta[5])
        covepseta = (1/(1 + exp(-theta[6])) - .5)*2*exp(theta[2])*exp(theta[5])

        # define locals
        covxieta = sigmaeta^2 - covepseta
        sigmaxi  = sqrt(sigmaeta^2 + sigmaeps^2 - 2*covepseta)
        xistar   = z*gamma - n*pibetan - kappa*betak

        #l_0
        l_0 = logcdf(Normal(0,1), -xistar./sigmaxi)

        #l_1
        pdfl1     = pdf(Normal(0,1), (w - z*gamma)/sigmaeta)
        arg1cdfl1 = xistar + (covxieta/(sigmaeta^2))*(w - z*gamma)
        arg2cdfl1 = sqrt(sigmaxi^2 - (covxieta^2)/(sigmaeta^2))
        cdfl1     = cdf(Normal(0,1), arg1cdfl1/arg2cdfl1)

        l_1 = log((1/sigmaeta)*pdfl1.*cdfl1)

        return l = -sum(d.*l_1 + (1-d).*l_0)
end

# Optimization
# initial condition
theta0 = [.5, log(.4), .4, .8, log(1), -log((2*.4*1)/(.3+.4*1) -1)]

# display optimization
optimal = optimize(llks, theta0, NelderMead())
mles = optimal.minimizer
mles[6] = (1/(1 + exp(-mles[6])) - .5)*2*exp(mles[2])*exp(mles[5])
mles[2] = exp(mles[2])
mles[5] = exp(mles[5])
println(mles)
