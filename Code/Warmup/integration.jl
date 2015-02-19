# Author: Jorge L. García
# Date: 7/1/2014

# Import packages
using Distributions

# Set a seed
srand(1)

# Plot 1000 points of sine(x) for one cycle of the unit circle
n = 1000
x = linspace(0, 2*pi, n)
s = sin(x)

# Trapezoidal integration
# function
function trapezoidalint(f,a,b,n)
        h = (b - a)/n
        l = f(a) + f(b)
        for k = 1:n
                x = a + k*h
                l += 2*f(x)
        end
        l *= h/2
        return l
end

# compare trapezoidal with closed form integration
# closed form
l = -cos(pi) + cos(0)
println(l)

# approximations
for n in [10,50,100,500,10000]
        l = trapezoidalint(sin,0,pi,n)
        println(l)
end

# Monte Carlo Integrarion
# function
function mcint(f,a,b,n)
    l = 0
    for k = 1:n
            x = rand(Uniform(a,b))
            l += f(x)
    end
    return l *= (b-a)/n
end

# approximations
for n in [10,50,100,500,10000]
        l = mcint(sin,0,pi,n)
        println(l)
end
