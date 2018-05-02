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
N = 1000; T = 6; betak = .5; betan = .2; sigmaeps = .4; pi = .2;
gamma1 = .8; sigmaeta = 1; covepseta = .3; gamma2 = .9; delta = .85;
covepsxi = covepseta - sigmaeps^2
covetaxi = sigmaeta^2 - covepseta
sigmaxi  = sqrt(sigmaeta^2 + sigmaeps^2 - 2*covepseta)

# simulation
y_lb = 0; y_ub = 10
z_lb = 0; z_ub = 5
k_lb = 0; k_ub = 5
n_lb = 0; n_ub = 3

# Simulate
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

# non-stochastic terms of utility
function W1(z,h,n)
        return gamma1*z + gamma2*h - pi*n
end

function W0(kappa,n)
        return betak*kappa + betan*n
end

# xistar, t = T - 0
function xistarTm0(z,h,n,kappa)
        return W1(z,h,n) - W0(kappa,n)
end

# Emax, t = T - 0
function EmaxTm0(z,h,n,kappa)
        cdf1 = cdf(Normal(0,1), xistarTm0(z,h,n,kappa)/sigmaxi)
        pdf0 = pdf(Normal(0,1), -xistarTm0(z,h,n,kappa)/sigmaxi)
        cdf0 = cdf(Normal(0,1), -xistarTm0(z,h,n,kappa)/sigmaxi)
        return W1(z,h,n).*cdf1 + W0(kappa,n).*cdf0 + sigmaxi*pdf0
end

# xistar, t = T - 1,..., 1
for t = 1:T-1
        eval(parse("function xistarTm" * string(t) * "(z,h,n,kappa)"
        * "W1(z,h,n) - W0(kappa,n)"
        * "+" * "delta*EmaxTm" * string(t-1) * "(z,h+1,n,kappa)" * "-" * "delta*EmaxTm" * string(t-1) * "(z,h,n,kappa)"
        * "end"))
end

# EmaxTmt,   t = T -1, ..., 2
for t = 1:T-1
        eval(parse("function EmaxTm" * string(t) * "(z,h,n,kappa)"
        * "W1(z,h,n) - W0(kappa,n)"
        * "+" * "delta*EmaxTm" * string(t-1) * "(z,h+1,n,kappa)" * ".*"
        * "cdf(Normal(0,1), xistarTm" * string(t) * "(z,h,n,kappa)/sigmaxi)"
        * "+" * "delta*EmaxTm" * string(t-1) * "(z,h,n,kappa)" * ".*"
        * "cdf(Normal(0,1), -xistarTm" * string(t) * "(z,h,n,kappa)/sigmaxi)"
        * "+" * "sigmaxi*pdf(Normal(0,1), -xistarTm" * string(t) * "(z,h,n,kappa)/sigmaxi)"
        * "end"))
end

#xistar vector function
xistarvf = [xistarTm0]
for t = 1:T-1
        xistarvf = [eval(parse("xistarTm" * string(t)));xistarvf]
end

# xistar, d, h
xistar = zeros(N,T); d = zeros(N,T); h = zeros(N,T)

for t = 1:T
        xistar[:,t] = xistarvf[t](z[:,t], h[:,t], n[:,t], kappa[:,t])
        for i = 1:N
                d[i,t] = eta[i,t] - eps[i,t] + xistar[i,t] > 0
        end
        if t < T
                h[:,t+1] = h[:,t] + d[:,t]
        end
end

# w
w = d.*(gamma1*z + gamma2*h + eta)

# Save to a .csv file
womandsample = [y z kappa n d w h]
womandsample = DataFrame(womandsample)
writetable("womandsample.csv",womandsample,header=true)
