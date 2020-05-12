rm(list=ls())

## run this chunk once to install package needed to pull 'pbpkme' from Github
install.packages("devtools")
library(devtools)
install_github("hillaryhusband/pbpkme")
############
library(pbpkme)

# Pharmacokinetics of Remifentanil: three-compartmental modeling approach, Cascone et al., Translational Medicine, 2013, 7(4), 18-22
# example describes 3 compartment model
# dCp/dt = - CL1*C1 + k21*V2*C2 + k31*V3*C3 - (k12+k13+k10)*C1)*V1 + dose) / V1
# dC2/dt = (k12 *C1*V1 - k21*C2*V2 - CL2*C2) / V2
# dC3/dt = (k13*C1*V1 - k31*C3*V3 - CL3*C3) / V3

#parameters
k10 <- 0.172    # min^-1
k12 <- 0.373    # min^-1
k13 <- 0.0367   # min^-1
k21 <- 0.103    # min^-1
k31 <- 0.0124   # min^-1

V1 <- 7.88 #mL
V2 <- 23.9 #mL
V3 <- 13.8 #mL

CL1 <- 2.08   # mL/min
CL2 <- 0.828  # mL/min
CL3 <- 0.0784 # mL/min

body_weight <- 70 # kg

#package parameters into one term
parameters <- data.frame(k10 = 0.172, k12 = 0.373, k13 = 0.0367, k21 = 0.103, k31 = 0.0124, V1 = 7.88, V2 = 23.9, V3 = 13.8, CL1 = 2.08, CL2 = 0.828, CL3 = 0.0784)

# example describes 3 compartment model
# dCp/dt = - CL1*C1 + k21*V2*C2 + k31*V3*C3 - (k12+k13+k10)*C1)*V1 + dose) / V1
# dC2/dt = (k12 *C1*V1 - k21*C2*V2 - CL2*C2) / V2
# dC3/dt = (k13*C1*V1 - k31*C3*V3 - CL3*C3) / V3

# as matrix
p_matrix <- matrix(0, nrow = 3, ncol = 3)

p_matrix[1,1] = -((k10 + k12 + k13)*V1 + CL1) / V1
p_matrix[1,2] = k21*V2 / V1
p_matrix[1,3] = k31*V3 / V1
p_matrix[2,1] = k12*V1 / V2
p_matrix[2,2] = - (k21*V2 + CL2) / V2
p_matrix[3,1] = k13*V1 / V3
p_matrix[3,3] = - (k31*V3 + CL3) / V2

p_matrix

# define dose
begin <- 0   # min
end <- 1    # min
dose_number <- 1
dose_amt <- 5 / 10 * body_weight #micrograms per kg

# initial conditions - amount of drug in compartment at time 0
initial_condition_x <- matrix(0, nrow = 3, ncol = 1)
initial_condition_x[1,1] = 0 #Cp
initial_condition_x[2,1] = 0 #C2
initial_condition_x[3,1] = 0 #C3

initial_condition_x


pred_time <- c(0.6836160344,
               1.2201418806,
               2.2681860729,
               3.5826361957,
               5.0317748779,
               8.0800919402,
               11.797147418,
               16.9744835334,
               28.9233821642,
               36.0935631794,
               45.7879153612,
               52.8278645992,
               60.0015120005)
pred_conc <- c(25.6697717548,
               33.3290050888,
               14.0424352764,
               6.0819984034,
               2.9001146879,
               1.7004870746,
               1.2424387791,
               0.8028673946,
               0.3132759634,
               0.186624412,
               0.0997417216,
               0.0690991259,
               0.0498866695)

pred_data <- data.frame(pred_time, pred_conc)



model <- pbpkme(60, 0.05, p_matrix, initial_condition_x, parameters, 1, dose_amt, begin, end)
model

time <- seq(from = 0, to = 60, by = 0.05)
Cp <- model[1,]
Cp_data <- data.frame(time,Cp)
Cp_data <- Cp_data[-1,]

pl1 <- ggplot(data = Cp_data, aes(x = time, y=Cp)) + geom_line() + scale_y_log10() + annotation_logticks(sides = "l")
pl2 <- pl1 + geom_line(data = pred_data, aes(x = pred_time, y = pred_conc, colour = "paper predictions"))
pl2
