# install the package if you haven't done so yet
install.packages('clustermq')

# load the library and create a simple function
library(clustermq)
fx = function(x) x * 2

# queue the function call on your scheduler
Q(fx, x=1:3, n_jobs=1)
# list(2,4,6)