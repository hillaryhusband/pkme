rm(list=ls())

## run this chunk once to install package needed to download 'pbpkme' from Github
## skip if you've already done this
install.packages("devtools")
# library(devtools)
# install_github("hillaryhusband/pkme")
############
# library(pkme)
source("R/pkme.R")
library(pracma)
library(readxl)
library(lattice)
library(ggplot2)
library(tidyverse)
library(ggthemes)

dt = 0.01
tfinal = 60

### Parameters and model based on Remifentanil model by Cascone et al.:
##  Cascone S, Lamberti G, Titomanlio G, Piazza O. Pharmacokinetics of Remifentanil:
## a three-compartmental modeling approach. Translational Medicine. 2013;7(4):18–22

#parameters
# k10 <- 0.172   min^-1
# k12 <- 0.373   min^-1
# k13 <- 0.0367  min^-1
# k21 <- 0.103   min^-1
# k31 <- 0.0124  min^-1
#
# V1 <- 7.88 mL
# V2 <- 23.9 mL
# V3 <- 13.8 mL
#
# CL1 <- 2.08    mL/min
# CL2 <- 0.828   mL/min
# CL3 <- 0.0784  mL/min

body_weight <- 79 # kg


params <- c(k10 = 0.172, k12 = 0.373, k13 = 0.0367, k21 = 0.103, k31 = 0.0124, V1 = 7.88,
            V2 = 23.9, V3 = 13.8, CL1 = 2.08, CL2 = 0.828, CL3 = 0.0784)
# drug injection
# R_in = 5 * body_weight

# C1 = 0
# C2 = 0
# C3 = 0

init_cond_matrix = matrix(0, nrow = 3, ncol = 1)

init_cond_matrix[1,1] = 0 #C1
init_cond_matrix[2,1] = 0 #C2
init_cond_matrix[3,1] = 0 #C3

#dC1_dt = (((-CL1*C1 + k21*V2*C2 + k31*V3*C3 - (k12+k13+k10)*C1)*V1) + R_in) / V1
#dC2_dt = (k12 *C1*V1 - k21*C2*V2 - CL2*C2) / V2
#dC3_dt = (k13*C1*V1 - k31*C3*V3 - CL3*C3) / V3


p_matrix = matrix(0, nrow=3, ncol=3)
p_matrix[1,1] = -((CL1/V1) + k12 + k13 + k10)
p_matrix[1,2] = (k21 * V2) / V1
p_matrix[1,3] = (k31 * V3) / V1
p_matrix[2,1] = (k12 * V1) / V2
p_matrix[2,2] = -((CL2/V2) + k21)
p_matrix[3,1] = (k13 * V1) / V3
p_matrix[3,3] = -((CL3/V3) + k31)


begin <- 0
end <- 5
dose_number <- 1
dose_amt <- 30


model1 <- pkme(minute = tfinal, h = dt, rate_coeff = p_matrix, init_condition = init_cond_matrix,
                 pars = params, number_of_doses = dose_number, dose = dose_amt, start_time = begin,
                 stop_time = end)

## create continuous time vector and pull out compartments from final matrix
time <- seq(from=0, to = tfinal, by = dt)
conc1 <- model1[1,]
conc2 <- model1[2,]
conc3 <- model1[3,]


## merge time vector and concentrations for plotting in ggplot
C1 <- data.frame(time, conc1)
C2 <- data.frame(time, conc2)
C3 <- data.frame(time, conc3)


plot1 <- ggplot(data = C1, aes(x = time, y = conc1))+ geom_line() + theme_base() +
   xlab("Time (minutes)") + ylab("Concentration (microMolar)") + ggtitle("Concentration vs Time", subtitle = "Compartment 1") + scale_y_log10()
plot2 <- ggplot(data = C1, aes(x = time, y = conc2))+ geom_line() + theme_base() +
   xlab("Time (minutes)") + ylab("Concentration (microMolar)") + ggtitle("Concentration vs Time", subtitle = "Compartment 2") + scale_y_log10()
plot3 <- ggplot(data = C1, aes(x = time, y = conc3))+ geom_line() + theme_base() +
   xlab("Time (minutes)") + ylab("Concentration (microMolar)") + ggtitle("Concentration vs Time", subtitle = "Compartment 3") + scale_y_log10()
