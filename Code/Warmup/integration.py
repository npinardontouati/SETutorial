# Author: Jorge L. García
# Date: 7/1/2014

#Import packages
import numpy as np
import matplotlib.pyplot as mplot
import os
#Set seed
np.random.seed(0)

#Change to my principal path
os.chdir(os.environ["jorge"])
os.chdir("econ350/StructuralEstimation/Warmup/Solution/Output")

#1.1 Basic continous function through a grid
n = 1000
x = np.linspace(0, 2*np.pi, n+1)
s = np.sin(x)

mplot.plot(x,s)
mplot.title("$\sin(x)$")
mplot.xlabel('x'); 
mplot.savefig("sin.png")
mplot.close()

#1.4 Trapezoidal integration
def trapezoidalint(f, a, b, n):
    h = (b-a)/float(n)
    I = f(a) + f(b)
    for k in xrange(1, n, 1):
        x = a + k*h
        I += 2*f(x)
    I *= h/2
    return I

#1.5 Trapezoidal integration of the sine function
# Calculate the closed form integral
cfi_sin = -np.cos(np.pi) + np.cos(0)
print cfi_sin
#Calculate the approximations
for n in 10, 50, 100, 500:
    approx_trapint_sin = trapezoidalint(np.sin, 0, np.pi, n)
    print approx_trapint_sin

#1.6 Monte Carlo Integration
def mcint(f, a, b, n):
    s = 0
    for i in range(n):
        x = np.random.uniform(a, b)
        s += f(x)
    I = (float(b-a)/n)*s
    return I

#1.7 Monte Carlo Integration of the sine function
for n in 10, 50, 100, 500, 1000:
    approx_mcint_sin = mcint(np.sin, 0, np.pi, n)
    print approx_mcint_sin