# Author: Jorge L. García
# Date: 7/1/2014

# Set seed
srand(1)

# Call number of processors
procs = 4
addprocs(procs)

# Define "to parallelize process"
require("womansmestimationf.jl")
H = 1000
h = 250

# Store
MPsamples = pmap(mestimationmles,[h,h,h,h])
mmles = vcat(MPsamples[1],MPsamples[2],MPsamples[3],MPsamples[4])

# Describe estimation in H samples
# still need to figure out kernel estimation

# describe optima via a matrix

# true parameters
tp = [.5, .4, .4, .8, 1, .3]

# mean, median, min, max, and std in H samples
tpH = zeros(5,1)

for p = 1:P
        tpH[p,1] = mean(mmlesm[:,p])
        tpH[p,2] = median(mmlesm[:,p])
        tpH[p,3] = minimum(mmlesm[:,p])
        tpH[p,4] = maximum(mmlesm[:,p])
        tpH[p,5] = std(mmlesm[:,p])
end

tpH = [tp tpH]
