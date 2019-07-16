
#This script takes simulation output and prepares data for stan models

library(rethinking)

# Load data

setwd("")  

Result <- read_excel("")


#Prepare data and create new variables from simulation output

dat <- as.list(Result)
dat$farm <- as.integer( ceiling( dat$trial / 25 ) )
dat$new_farm <- sapply( 1:nrow(Result) , function(i) {
	if ( i>1 ) {
		return( ifelse(dat$farm[i]!=dat$farm[i-1],1,0) )
	} else {
		return(0)
	}
} )


# prep id's
dat$group_id <- dat$group.id
dat$group.id <- NULL
dat$sid <- dat$id
dat$id <- (dat$Session_ID - 1)*8 + dat$id
dat$N_id <- max(dat$id)
dat$N <- nrow(Result)

#prep matrices with choices and ages
nmat <- cbind( dat$n1 , dat$n2 , dat$n3 , dat$n4 )
choice_models <- cbind(dat$Choice1, dat$Choice2, dat$Choice3)
age_models    <- cbind(dat$Age1, dat$Age2, dat$Age3)

choice_models[is.na(choice_models)] <- 0
age_models[is.na(age_models)] <- 0

dat$nmat <- nmat
dat$choice_models <- choice_models
dat$age_models <- age_models


#Experience matric for Gaussian processes
dat$expmat <- matrix(nrow = 20, ncol = 20)
for (i in 1:20) {
  for (j in 1:20) {
    dat$expmat[i,j] <- abs(i-j)
  }
}

# Run different stan models

m1 <- stan( file="ewa_model1.stan" , data=dat , chains=1 ) #Individual learning only

m2 <- stan( file="ewa_model2.stan" , data=dat , chains=1 )  #Simplest model with individual and (unbiased) social learning

m3 <- stan( file="ewa_model3.stan" , data=dat , chains=1 )  #Social Learning with frequency and experience bias

m4 <- stan( file="ewa_model10.stan" , data=dat , chains=1 )  #Full Gaussian process regression

m5 <- stan( file="ewa_model13.stan" , data=dat , chains=1 ) # All functional

m6 <- stan( file="ewa_model14.1.stan" , data=dat , chains=1 )  # All Monotonic effects



