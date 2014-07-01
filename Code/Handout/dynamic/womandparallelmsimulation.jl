# Author: Jorge L. García
# Date: 7/1/2014

# Import packages
using DataFrames


# Set seed
srand(1)

# Call number of processors
procs = 4
addprocs(procs)

# Define "to parallelize process"
require("womandmsimulationf.jl")
H = 1000
h = 250

# Store
MPsamples = pmap(womandmsimulation,[h,h,h,h])
womandmsamples = vcat(MPsamples[1],MPsamples[2],MPsamples[3],MPsamples[4])

# Save to a .csv file
womandmsamples = DataFrame(womandmsamples)
writetable("womandmsamples.csv",womandmsamples,header=true)
