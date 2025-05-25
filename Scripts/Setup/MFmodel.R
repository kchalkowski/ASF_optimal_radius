# Load package
#library(deSolve)

# Define the SEIR model with control
seir_model <- function(thyme, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Force of infection
    lambda <- beta * I / N
    
    # ODEs with removal applied to E and I
    dS <- -lambda * S
    dE <- lambda * S - sigma * E - removal_rate * E
    dI <- sigma * E - gamma * I - removal_rate * I
    dR <- gamma * I
    
    return(list(c(dS, dE, dI, dR)))
  })
}

# Parameters
parameters <- c(
  beta = 0.6,         # transmission rate
  sigma = 1/5.2,      # incubation rate (1/mean latent period)
  gamma = 1/7,        # recovery rate (1/mean infectious period)
  removal_rate = 0.1  # control removal rate applied to E and I
)

# Initial state
init_state <- c(
  S = 990,
  E = 5,
  I = 5,
  R = 0
)

# Total population (for force of infection calculation)
parameters["N"] <- sum(init_state)

# Time vector
times <- seq(0, 100, by = 1)

# Solve the model
output <- ode(y = init_state, times = times, func = seir_model, parms = parameters)

# Convert to data frame for plotting
output_df <- as.data.frame(output)

# Plot
matplot(output_df$thyme, output_df[ , 2:5], type = "l", lty = 1, lwd = 2,
        xlab = "Time (days)", ylab = "Population", col = c("blue", "orange", "red", "green"))
legend("right", legend = c("S", "E", "I", "R"),
       col = c("blue", "orange", "red", "green"), lty = 1, lwd = 2)
