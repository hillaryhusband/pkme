# pkme
This package is designed for classic compartmental or physiologically-based pharmacokinetic models. These types of models are most frequently described by systems 
of linear ordinary differential equations ODEs. This package takes a matrix-based approach by finding the eigenvalue/eigenvector solution to a linear system of 
ODEs. Matrix manipulation is computationally fast in R, making it ideal for model fitting/parameter estimation, particularly for large models. This function most 
benefits models/systems with high dimensionality and/or problems with stiffness. The manner in which the exact solutions are solved for and integrated allows for 
relatively accurate estimates while using larger step sizes than traditional methods - increasing computation speed. We do not claim that this method is faster than
traditional methods for all models, but it is a nice option for difficult models.  

Structure:
It is designed to take an initial condition matrix (single column vector of initial concentrations) 
and a rate coefficient matrix (flow/volume), dose profile (number of doses, start and stop time of each dose)
along with any parameters needed to solve the system, and outputs a matrix of the concentrations in each 
compartment (rows) versus the time course (columns) for a given step size `h`. Each column represents a 
predicted data point at every `2*h`. 

Additional details:

**Rate Coefficient Matrix** `(rate_coeff)` : your flow per volume and any other parameters and microconstants included in the equations – as you’re looking at the matrix, rows represent compartments and columns represent concentrations

 **Initial Conditions Matrix** `(init_condition)`: just your initial conditions or in this case, initial concentrations in each compartment. Normally a vector of zeroes, but if you did have existing drug concentrations in some organs, it could be represented here
 
**Dose Vector** `(dose)` : where the dose is designated – the nice thing about this function is you can dose into any compartment just by specifying the row in the vector corresponding to that compartment. Let’s say you’re dosing IM and your muscle compartment is #4, you would place your dose in this vector in row 4. Homogenous and non-homogenous linear systems are handled using the same function. For homogenous systems, the "g vector" corresponding to the dose profile can be replaced by setting the total dose as the initial condition for the first compartment and the dose amount in the function to 0.


Also, some first order linear systems produce complex conjugate pairs in the eigenvalues or eigenvectors of the matrix solution. This function handles those complex conjugate pairs automatically.

# Installation Shortcut

Install the latest release on CRAN with:
```
install.packages("pkme")
```

Install the current development version from GitHub with:
```
# install.packages("devtools")
devtools::install_git("git://github.com/hillaryhusband/pkme.git")
```

# Examples:

## One-compartment IV bolus:

```
library(pkme)
```

Model structure:
dC/dt = -k10*C

Designate PK parameters
```
CL = 2.08 #(mL/min)
V = 7.88 #(mL)
k10 = CL/V #(1/min)
# create parameter vector to feed into solver
params <- c(CL, V, k10)
```

Create matrix form of model
```
pmat = matrix(0, nrow=1, ncol=1)
pmat[1,1] = -k10
```

Initialize initial conditions matrix
```
init = matrix(0, nrow=1, ncol=1)
```

Designate time:
```
tfinal = 10  # end of analysis (min)
dt = 0.5     # time step (min)
```

Set-up dose:
```
begin <- 0
end <- dt # for bolus, set infusion rate to length of time step
dosen <- 1
doseamt <- 30
```

Solve system:
```
model1 <- pkme(minute = tfinal,
               h = dt,
               rate_coeff = pmat,
               init_condition = init,
               pars = params,
               number_of_doses = dosen,
               dose = doseamt,
               start_time = begin,
               stop_time = end)
```

Output prints concentration at each time step for each compartment in the system. This code is available in the R directory as `R/pkme-bolus-example`. 

***NOTE: Three-compartment PK model vignette is available in the R directory as `R/pkme-infusion-example.R`.***
