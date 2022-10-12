############Generating Fake F2_FL data for glm 
#(until get real dataset/code to fit actual model in R)

n_sim <- 10000           # for the initial dataset
#set.seed(2)            # for replicability
xx <- runif(n_sim)     # predictor values
coefficients <- c(-0.80094,-1.9128) # my assumption
prob <- 1/(1+exp(-(coefficients[1]+coefficients[2]*xx)))

yy <- runif(n_sim)<prob

F2_FL <- glm(yy~xx,family="binomial")

############Generating Fake F2i_FL data for glm 
#(until get real dataset/code to fit actual model in R)

n_sim <- 10000           # for the initial dataset
#set.seed(2)            # for replicability
xx <- runif(n_sim)     # predictor values
coefficients <- c(2.3215,-6.3723) # my assumption
prob <- 1/(1+exp(-(coefficients[1]+coefficients[2]*xx)))

yy <- runif(n_sim)<prob

F2i_FL <- glm(yy~xx,family="binomial")

############Generating Fake F2_SC data for glm 
#(until get real dataset/code to fit actual model in R)

n_sim <- 10000           # for the initial dataset
#set.seed(2)            # for replicability
xx <- runif(n_sim)     # predictor values
coefficients <- c(0.63955,-1.1612) # my assumption
prob <- 1/(1+exp(-(coefficients[1]+coefficients[2]*xx)))

yy <- runif(n_sim)<prob

F2_SC <- glm(yy~xx,family="binomial")

############Generating Fake F2i_SC data for glm 
#(until get real dataset/code to fit actual model in R)

n_sim <- 10000           # for the initial dataset
#set.seed(2)            # for replicability
xx <- runif(n_sim)     # predictor values
coefficients <- c(9.82,-4.107) # my assumption
prob <- 1/(1+exp(-(coefficients[1]+coefficients[2]*xx)))

yy <- runif(n_sim)<prob

F2i_SC <- glm(yy~xx,family="binomial")
