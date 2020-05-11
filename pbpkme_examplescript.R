rm(list=ls())

## run this chunk once to install package needed to pull 'pbpkme' from Github
install.packages("devtools")
library(devtools)
install_github("hillaryhusband/pbpkme")
############
library(pbpkme)

# example describes 3 compartment model
# dC1 = dose - (k10 + k12 + k13) * C1 + k21 * C2 + k31 * C3
# dC2 = k12 * C1 - k21 * C2
# dC3 = k13 * C1 - k31 * C3

#parameters
k10 <- 0.2 # clearance of drug from main compartment
k12 <- 0.5 # volume of distribution of main compartment
k13 <- 0.1
k21 <- 0.14
k31 <- 0.4

#package parameters into one term
parameters <- data.frame(k10 = 0.2, k12 = 0.5, k13 = 0.1, k21 = 0.14, k31 = 0.4)

# as matrix
p_matrix <- matrix(0, nrow = 3, ncol = 3)

p_matrix[1,1] = -(k10 + k12 + k13)
p_matrix[1,2] = k21
p_matrix[1,3] = k31
p_matrix[2,1] = k12
p_matrix[2,2] = -k21
p_matrix[3,1] = k13
p_matrix[3,3] = -k31

p_matrix

# define dose
begin <- 0
end <- 5
dose_number <- 1
dose_amt <- 3


# initial conditions - amount of drug in compartment at time 0
initial_condition_x <- matrix(0, nrow = 3, ncol = 1)
initial_condition_x[1,1] = 0 #C1
initial_condition_x[2,1] = 0 #C2
initial_condition_x[3,1] = 0 #C3

initial_condition_x






model <- pbpkme(50, 0.1, p_matrix, initial_condition_x, parameters, 1, 60, 0, 1)
model
