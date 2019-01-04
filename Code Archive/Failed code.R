# FAILED VERSION #
#### Bootstrap Bandwidth Selection (Farraway-Jhun 1990) ####
# Calculate bootstrap based bandwidth.
# x is sample data
# h0 is pilot bandwidth
# B is bootstrap sample size.
# upper.limit is the upperbound for the optimize call.

# THIS PROGRAM IS REALLY SLOW AND PROBABLY MESSED UP
boot.FJ.bw <-function(x, h0 = NA, B = 1000, m = NA, suff = F, upper.limit = 5){
	
	# Initial Smoothing Parameter
	if(is.na(h0)) h0 <- std.bw(x, "SJ")
	
	# Smoothed bootstrap sample)
	xx <- boot.smooth(x, B = B, m = m, h = h0, suff = suff)
	fx <- fhat(x, h0)
	#fxx <- lapply(xx, fhat, h = h0) # creating this list is being verrrrry problematic, kept as
	#legacy code
	
	
	bigfun<-function(x,h){
		rowSums(sapply(fxx, function(f) (f(x,h) - fx(x))^2) )} #ACTUALLY SLOW PART HEREbw
	
	
	bmise.fj <- Vectorize(function(h){ 
		
		integrate(bigfun, lower = -Inf, upper = Inf, h)$value #SLOW PART ALSO HERE
		
	}, "h")
	
	optimize(bmise.fj, c(0,1))$minimum # SLOW PART IS HERE
	
	
}


suff.bw <-function(x, h0 = NA, B = 1000, suff = F, upper.limit = 5){
	
}