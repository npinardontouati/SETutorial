# Author: Jorge L. García
# Date: 7/1/2014

# Set seed
srand(1)

# Call number of processors
procs = 4
addprocs(procs)

# Define "to parallelize process"
require("womansbootstrapf.jl")
B = 1000
b = 250

# Store
MPsamples = pmap(bootstrapmles,[b,b,b,b])
bmles = vcat(MPsamples[1],MPsamples[2],MPsamples[3],MPsamples[4])

# Obtain standard errors through sd of B samples
bsemles = zeros(P,1)
for p = 1:P
        bsemles[p,1] = std(bmles[:,p])
end
