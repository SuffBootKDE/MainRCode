# OBJECTIVES
# Make it so bandwidth scales

h.opt.fn(n = 10)

ddnorm <- function(x) (x^2-1)*dnorm(x)

r.int(ddnorm)


(r.int(dnorm)/(1^2*r.int(ddnorm)*10))^(1/5)

