# Author: Jorge L. García
# Date: 7/1/2014

# Import packages
using Distributions
using Optim
using DataFrames

# Set a seed
srand(1)

# Parameters
# model
N = 1000; T = 6; betak = .5; betan = .2; sigmaeps = .4
pi = .2; gamma = .8; sigmaeta = 1; covepseta = .3

# simulation
y_lb = 0; y_ub = 10
z_lb = 0; z_ub = 5
k_lb = 0; k_ub = 5
n_lb = 0; n_ub = 3

# samples, parameters, variables
H = 1000 ; P = 6; V = 6

# preallocate space
womansmsamples = zeros(1,T*V)

for h = 1:H
        # simulate
        # y
        y = zeros(N,T)
        for t = 1:T
                y[:,t] = [rand(Uniform(y_lb,y_ub)) for i = 1:N]
        end

        # z,kappa,n (static over time)
        z = zeros(N,T)
        kappa = zeros(N,T)
        n = zeros(N,T)

        z_N = [rand(Uniform(z_lb,z_ub)) for i = 1:N]
        kappa_N = [rand(Uniform(k_lb,k_ub)) for i = 1:N]
        n_N = [rand(Uniform(n_lb,n_ub)) for i = 1:N]

        # z,kappa,n (wide)
        for t = 1:T
                z[:,t] = z_N
                kappa[:,t] = kappa_N
                n[:,t] = n_N
        end

        # unobserved
        eps = zeros(N,T)
        eta = zeros(N,T)

        for t = 1:T
                epseta =  [rand(MvNormal([sigmaeps^2 covepseta; covepseta sigmaeta^2])) for i = 1:N]
                for i = 1:N
                        eps[i,t] = epseta[i,1][1]
                        eta[i,t] = epseta[i,1][2]
                end
        end

        w = gamma*z + eta
        U1 = y + gamma*z + eta - pi*n
        U0 = y + betak*kappa + betan*n + eps
        v  = U1 - U0

        d = zeros(N,T)
        for t = 1:T
                for i = 1:N
                        d[i,t] = v[i,t] > 0
                end
        end

        #store sample h
        sampleh = [y z kappa n d w]
        womansmsamples = [womansmsamples, sampleh]

end

# Save to a .csv file
womansmsamples = womansmsamples[2:end,:]
womansmsamples = DataFrame(womansmsamples)
writetable("womansmsamples.csv",womansmsamples,header=true)
