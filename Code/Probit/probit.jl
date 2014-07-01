# Author: Jorge L. García
# Date: 7/1/2014

# Import packages
using Distributions
using Optim

# Set a seed
srand(1)

# Probit
# data generation
N = 10000
beta = [3,1.5]

x1 = [rand(Normal(0,1.5)) for i = 1:N]
x2 = [rand(Uniform(0,1)) for i = 1:N]
X = [x1 x2]

E = [rand(Normal(0,1)) for i = 1:N]

Ystar = X*beta + E
d = [Ystar[i] > 0 for i = 1:N]

# negatve of the likelihood
function llkprobit(B)
        l_0 = logcdf(Normal(0,1), -X*B)
        l_1 = logcdf(Normal(0,1),  X*B)
        for i = 1:N
                if l_0[i] == -Inf
                        l_0[i] = log(1e-300)
                elseif l_1[i] == -Inf
                        l_1[i] = log(1e-300)
                end
        end
        return l = -sum(d.*l_1 + (1-d).*l_0)
end

optimal = optimize(llkprobit, beta, method = :l_bfgs)
mleprobit = optimal.minimum
println(mleprobit)
