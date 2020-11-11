rm(list=ls())

## run this chunk once to install package needed to pull 'pbpkme' from Github
install.packages("devtools")
library(devtools)
install_github("hillaryhusband/pbpkme")
############
library(pbpkme)

working.directory <- ("C:/Users/HH/Dropbox/Research - Hillary")
data1 <- read_excel(file.path(working.directory, "DataIndex.xlsx"), "3cptex2")
data1 <- as.data.frame(data1[,1:2])

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

body_weight <- 79 # kg (mean value)

#package parameters into one term
parameters <- data.frame(k10 = 0.172, k12 = 0.373, k13 = 0.0367, k21 = 0.103, k31 = 0.0124, V1 = 7.88, V2 = 23.9, V3 = 13.8, CL1 = 2.08, CL2 = 0.828, CL3 = 0.0784)

# example describes 3 compartment model
# dCp/dt = - CL1*C1 + k21*V2*C2 + k31*V3*C3 - (k12+k13+k10)*C1)*V1 + dose) / V1
# dC2/dt = (k12 *C1*V1 - k21*C2*V2 - CL2*C2) / V2
# dC3/dt = (k13*C1*V1 - k31*C3*V3 - CL3*C3) / V3

# as matrix
p_matrix <- matrix(0, nrow = 3, ncol = 3)

 # p_matrix[1,1] = -((k10 + k12 + k13)*V1 + CL1) / V1
 # p_matrix[1,2] = (k21*V2) / V1
 # p_matrix[1,3] = (k31*V3) / V1
 # p_matrix[2,1] = (k12*V1) / V2
 # p_matrix[2,2] = - (k21*V2 + CL2) / V2
 # p_matrix[3,1] = (k13*V1) / V3
 # p_matrix[3,3] = - (k31*V3 + CL3) / V3

p_matrix = matrix(0, nrow=3, ncol=3)
p_matrix[1,1] = -((CL1/V1) + k12 + k13 + k10)
p_matrix[1,2] = (k21 * V2) / V1
p_matrix[1,3] = (k31 * V3) / V1
p_matrix[2,1] = (k12 * V1) / V2
p_matrix[2,2] = -((CL2/V2) + k21)
p_matrix[3,1] = (k13 * V1) / V3
p_matrix[3,3] = -((CL3/V3) + k31)

p_matrix

# define dose
begin <- 0   # min
end <- 0.5     # min
dose_number <- 1
dose_amt <- 30 *body_weight / V1 # 1 * body_weight / V1 #micrograms per kg

# initial conditions - amount of drug in compartment at time 0
initial_condition_x <- matrix(0, nrow = 3, ncol = 1)
initial_condition_x[1,1] = 0 #Cp
initial_condition_x[2,1] = 0 #C2
initial_condition_x[3,1] = 0 #C3

initial_condition_x
#
#
# pred_time <- c(1.10080446,
#                2.0084653157,
#                2.9121645872,
#                3.972341131,
#                5.0325176749,
#                6.0921990207,
#                7.0186774013,
#                8.0793491432,
#                9.9340390974,
#                11.9224272149,
#                14.9742106634,
#                19.8858836453,
#                25.0634673598,
#                29.9803399209,
#                39.8138374443,
#                59.8863837631)
# pred_conc <- c(70.9108628188,
#                19.7896272624,
#                4.4337123319,
#                3.660451177,
#                3.0220505563,
#                2.4274174178,
#                1.9230709498,
#                1.6318745996,
#                1.1275187146,
#                0.8118522503,
#                0.5769069087,
#                0.3779003773,
#                0.2475758807,
#                0.2163643909,
#                0.1629965866,
#                0.1340392725)
#
# pred_data <- data.frame(pred_time, pred_conc)


tic()
model <- pbpkme(100, 0.05, p_matrix, initial_condition_x, parameters, 1, dose_amt, begin, end)
toc()


Cp <- model[1,]
pred_time <- seq(0,100, by = 0.05)
pred <- cbind(pred_time,Cp)
pred<-as.data.frame(pred)
obs_time <- as.matrix(data1$Minute)
obs_conc <- as.matrix(data1$microMolar)
obs <- cbind(obs_time,obs_conc)
obs <- as.data.frame(obs)

pl1 <- ggplot(data = obs, aes(x=V1, y=V2, colour = "Observed Concentration")) + geom_point(size = 2) +
   scale_y_log10() + ggtitle("Conc vs Time Curve of Remifentanil", subtitle = "30 microgram/kg bolus") +
   theme_classic(base_size = 14) + xlab("Time (minutes)") + ylab("Concentration (microMolar)") +
   theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) +
   theme(legend.title = element_blank(), legend.position = c(0.7,0.8)) + annotation_logticks(sides = "l")
pl2 <- pl1 + geom_line(data = pred, aes(x=pred_time, y=Cp, colour = "Model Predictions"), size = 0.8)
pl2
