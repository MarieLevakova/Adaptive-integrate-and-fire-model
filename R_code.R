#-----------------------------------------------------------------------------------------------
# Generates the spike train in response to a step stimulus of an odorant of given concentration
# and plots the kernel smoothed firing
#-----------------------------------------------------------------------------------------------

rm(list=ls())

source("Cpp_functions.R") # Contains a Cpp function to simulate the model

# Initialization of parameters
receptor.pars <- c(ki = 1e6, 
                   k1 = 0.209, 
                   km1 = 7.9, 
                   k2 = 16.8, 
                   km2 = 98, 
                   k3 = 100, 
                   km3 = 98.9, 
                   k4 = 4e4, 
                   Rtot = 1.64, 
                   Ntot = 1)

LIF.pars = c(Cm = 0.00144, 
             gL = 1.44, 
             EL = -62, 
             ER = 0, 
             Vreset = -62, 
             theta0 = -55)

fit.pars = c(tau = 0.59512792, 
             Delta = 0.77496971, 
             gamma = 98.38484201,
             n = 0.05517291)

concentration <- c(1e-7, 1e-6, 1e-5, 1e-4) # Concentration of the odorant in micro molar

# Simulation of the model for a step stimulus
time.step <- 2e-6
time.grid <- seq(0, 0.5, by = time.step)
L.air <- stepfun(0, c(0,1)) # Assuming a step stimulus with onset at 0
valve.states <- L.air(time.grid) 

plot(c(-0.05, 0.5), c(0, 85), type = "n", xlab = "Time [s]", ylab = "Firing rate")

for(i in 1:length(concentration)){
  # Generate a vector of 0s and 1s indicating occurrence of a spike at each time
  trace <- simAdp(conc = concentration[i],
                  receptorPars = receptor.pars, LIFPars = LIF.pars, fitPars = fit.pars, 
                  times = time.grid, 
                  initial = c(0, 0), 
                  valveStates = valve.states)
  
  # Extract the times of spikes
  spike.times <- time.grid[which(trace == 1)]
  
  # Kernel estimate of the firing rate
  rate.kernel <- density(spike.times, bw = 0.03, from = 0, to = 0.5)
  lines(rate.kernel$x, rate.kernel$y*sum(trace), lty = 2, col = i)
}

