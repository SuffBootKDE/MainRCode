#### THE GOLD CODE #### HOW TO SOLVE THE PROBLEM #RUNS VERY SLOW...
funs <- fxx

bigfun<-function(x,funs){
	rowSums(sapply(funs, function(f) f(x)))}

bigfun(x=c(-1,0,1),funs=funs)

integrate(bigfun, lower= -Inf, upper = Inf, funs)
#### END GOLD CODE ####


#### SUBSIDIARY CODE ####
bigfun2<-function(vec,funs){vf<-function(vec,funs){funs(vec)} 
sum(sapply(1:length(vec),function (i) vf(vec[i],funs[[i]])))}

optim(par=initvec,fn=bigfun2,funs=funs)

#### Single MISE calculation formula? ####
# MAYBE NOT RIGHT NOW