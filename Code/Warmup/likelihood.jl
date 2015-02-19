# Author: Jorge L. García
# Date: 7/1/2014

# Import packages
using Distributions
using Optim

# Set a seed
srand(1)

# Ronsensbrock function
function rosenbrock(x)
        xf = x[2:end]
        xl = x[1:end-1]
        return sum((1-x).^2) + sum(100*(xf - xl.^2).^2)
end

# optimize
optimal = optimize(rosenbrock, zeros(10), method = :l_bfgs)
rosenbrockmin = optimal.minimum
println(rosenbrockmin)

# Likelihood (normal mean zero unknown variance = 1)
# generate data
N = 10000
x = rand(Normal(0,1),N)

# closed form mle
mle = sum(x.^2)/N
println(mle)

# negative of likelihood
function llknorm(theta)        
        # transform parameter to avoid evaluating
        # outside function domain
        theta  = exp(theta)
        l = log(theta) + (1./(N*theta))*sum(x.^2)
        return l[1]
end

optimal = optimize(llknorm, [1.7], method = :l_bfgs)
mleopt = exp(optimal.minimum)
println(mleopt)
