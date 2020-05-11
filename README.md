# pbpkme
This package is designed for classic compartmental or physiologically-based pharmacokinetic models. 
These types of models are most frequently described by systems of ordinary differential equations ODEs. This package
takes a matrix-based approach by finding the eigenvalue/eigenvector solution to a linear system of ODEs.
Matrix manipulation is more computationally fast, making it ideal for model fitting/parameter estimation.

Structure:
It is designed to take an initial condition matrix (single column vector of initial concentrations) 
and a rate coefficient matrix (flow/volume), dose profile (number of doses, start and stop time of each dose)
along with any parameters needed to solve the system, and outputs a matrix of the concentrations in each 
compartment (rows) versus the time course (columns) for a given step size "h". Each column represents a 
predicted data point at every "2*h". 
