#library(pkme)

# Model structure:

# dC/dt = -k10*C

# Designate PK parameters
CL = 2.08 #(mL/min)
V = 7.88 #(mL)
k10 = CL/V #(1/min)

params <- c(CL, V, k10)

# Create matrix form of model

pmat = matrix(0, nrow=1, ncol=1)
pmat[1,1] = -k10

# Initialize initial conditions matrix

init = matrix(0, nrow=1, ncol=1)

# Designate time:
tfinal = 10  # end of analysis (min)
dt = 0.5     # time step (min)


# Set-up dose:
begin <- 0
end <- dt # for bolus, set infusion rate to length of time step
dosen <- 1
doseamt <- 30


# Solve system:

model1 <- pkme(minute = tfinal,
               h = dt,
               rate_coeff = pmat,
               init_condition = init,
               pars = params,
               number_of_doses = dosen,
               dose = doseamt,
               start_time = begin,
               stop_time = end)
