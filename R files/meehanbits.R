
# Meehan approximation notes - used in NmixtureS and NmixtureV from nimble 
# This was for me to understand how it worked 

# original paper - check appendix B
# https://www.jstatsoft.org/article/view/v095i02

lambda=8 # set a lambda
p=0.8 # set a p
# say you saw 10 individuals, then for the upper 
# bound of the summation you can 
# use a high quartile of the distribution of individuals 
# you didnt see given a set p 
10 + qpois(0.999, lambda * (1-p))
hist(qpois(seq(0.1,0.999,0.001), lambda * (1-p))) # e.g conditional distribution of missed individuals


# trial of approximation
meehan = function(lamda){
 p=0.6 # fixed p
 fac = 1
 lambda=lamda; 
 ff = lambda * (1-p) 
 y = 10 #  count is 10
 N.max = y + qpois(0.999, lambda * (1-p))

 for(i in (N.max - y):1) fac <- 1 + fac * ff / i # loop through the recursive sum

 log.L <- dpois(y, lambda, log = TRUE) + 
  + dbinom(y, y, p, log = TRUE) + log(fac) # calc log likelihood 

 return(log.L)
}

lat =c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)
plot(lat,sapply(lat,meehan))
abline(v=lat[which.max(sapply(lat,meehan))]) # likelihood maximized at 17 as true value

16*0.6
17*0.6 # closest to 10
18*0.6

# without the meehan component 
meehan2 = function(lamda){
  p<-0.6
  fac <- 1
  lambda<- lamda; 
  ff <- lambda * (1-p) 
  y <- 10 #  count
  N.max <- y + qpois(0.999, lambda * (1-p))
  
  for(i in (N.max - y):1) fac <- 1 + fac * ff / i
  log.L <- dpois(y, lambda, log = TRUE) +
    + dbinom(y, y, p, log = TRUE) 
  
  return(log.L)
}

plot(lat,sapply(lat,meehan2))
abline(v=lat[which.max(sapply(lat,meehan2))]) 
