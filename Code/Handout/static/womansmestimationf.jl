# Author: Jorge L. García
# Date: 7/1/2014

# Import packages
using Distributions
using Optim
using DataFrames

# Set a seed
srand(1)

# Read data
womansmsamples = matrix(readtable("womansmsamples.csv"))

# Define variables: individuals, parameters, periods, variables
N = 1000; P = 6; T = 6; V = 6

# Optimize for each sample
# define likelihood
function llks(theta,z,kappa,n,d,w)
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

# initial condition
theta0 = [.5, log(.4), .4, .8, log(1), -log((2*.4*1)/(.3+.4*1) -1)]

function mestimationmles(H)

        # preallocate space
        mlesm = zeros(P,H)

        for h = 1:H

                #define sample h
                hh = h - 1
                i  = (hh*N)+1
                j  = h*N

                y = womansmsamples[i:j,1:T]
                z = womansmsamples[i:j,T+1:2*T]
                kappa = womansmsamples[i:j,2*T+1:3*T]
                n = womansmsamples[i:j,3*T+1:4*T]
                d = womansmsamples[i:j,4*T+1:5*T]
                w = womansmsamples[i:j,5*T+1:6*T]

                # wrapper function
                function wllks(theta)
                        return llks(theta,z,kappa,n,d,w)
                end

                # optimization
                optimal = optimize(wllks, theta0, NelderMead())
                mles = optimal.minimizer
                mles[6] = (1/(1 + exp(-mles[6])) - .5)*2*exp(mles[2])*exp(mles[5])
                mles[2] = exp(mles[2])
                mles[5] = exp(mles[5])

                # store h sample optimum
                mlesm[:,h] = mles
        end
        return mlesm'
end
