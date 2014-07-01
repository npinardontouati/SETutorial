# Author: Jorge L. García
# Date: 7/1/2014

# Import packages
using Plotly                                 #this package needs to be registred
Plotly.signin("jorgelgarcia", "h2ajfsbinv")  #place yout username and your password here once you register

# Still need to figure out LaTeX #
# Is currently crashing with $ marks #

# Transformation for symmetric intervals

# Parameters
x = linspace(-100, 100, 1000)
f10 = (1/(1 + exp(-x)) - .5)*2*10
f20 = (1/(1 + exp(-x)) - .5)*2*20

# Plot

# simple plot
layout = [
  "title" => "Symmetric Intervals Transformation: f(x) = 2a ( 1/(1+exp(-x) -.5 )",
  "xaxis" => [
    "title" => "x",
    "showgrid" => true,
    "zeroline" => false
  ],
  "yaxis" => [
    "title" => "f(x)",
    "showline" => false
  ]
]

trace10 = [
  "x" => [x],
  "y" => [f10],
  "name" => "a = 10",
  "type" => "scatter"
]

trace20 = [
  "x" => [x],
  "y" => [f20],
  "name" => "a = 20",
  "type" => "scatter"
]

trace = [trace10, trace20]

plot = Plotly.plot([trace], ["layout" => layout])
ploturl = plot["url"]
