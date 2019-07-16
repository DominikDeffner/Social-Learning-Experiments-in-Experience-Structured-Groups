
# Code to compute and plot implied probabilities of choosing option 1 of 4 according 
# to 3 different ways to combine conformity and experience (here called "Age")
# see preregistration document for details


# First we construct available social information

Info <- data.frame(Ind1=rep(0,2), Ind2=rep(0,2), Ind3=rep(0,2), row.names = c("Choices", "Ages"))

# Which of options 1-4 did 3 models choose in previous round?

Info["Choices",1] <- 1
Info["Choices",2] <- 1
Info["Choices",3] <- 4


#What is their experience ("Age") in the present region?

Info["Ages",1] <- 5
Info["Ages",2] <- 10
Info["Ages",3] <- 4

#Frequency of each option among models

n1 = length(which(Info["Choices",]==1))
n2 = length(which(Info["Choices",]==2))
n3 = length(which(Info["Choices",]==3))
n4 = length(which(Info["Choices",]==4))


#3 options to combine conformity and experience/age bias

# 1) All combined into one expression with experience per option as cue similar to frequency of that option
  
functComb <- function(x){
  
  #Compute average experience PER OPTION if it was chosen
  
  if (n1 != 0){
    age1 = sum(Info["Ages",][Info["Choices",]==1]) / n1
  } else {
    age1 = 0
  }
  
  if (n2 != 0){
    age2 = sum(Info["Ages",][Info["Choices",]==2]) / n2
  } else {
    age2 = 0
  }
  
  if (n3 != 0){
    age3 = sum(Info["Ages",][Info["Choices",]==3]) / n3
  } else {
    age3 = 0
  }
  
  if (n4 != 0){
    age4 = sum(Info["Ages",][Info["Choices",]==4]) / n4
  } else {
    age4 = 0
  }
  
  (n1^x)*exp(b*age1) / ((n1^x)*exp(b*age1)+(n2^x)*exp(b*age2)+(n3^x)*exp(b*age3)+(n4^x)*exp(b*age4))
} 

# 2) Pseudocounts for individuals

functPseudo <- function(x) {

#Fist we calculate pseudo counts for all 4 options. 
  
Pseudo_n1 <- 0
for (i in Info["Ages",][Info["Choices",]==1]) {
  Pseudo_n1 <-  Pseudo_n1 + exp(b*i)
}

Pseudo_n2 <- 0
for (i in Info["Ages",][Info["Choices",]==2]) {
  Pseudo_n2 <-  Pseudo_n2 + exp(b*i)
}

Pseudo_n3 <- 0
for (i in Info["Ages",][Info["Choices",]==3]) {
  Pseudo_n3 <-  Pseudo_n3 + exp(b*i)
}

Pseudo_n4 <- 0
for (i in Info["Ages",][Info["Choices",]==4]) {
  Pseudo_n4 <-  Pseudo_n4 + exp(b*i)
}

#Pseudocounts are then used in standard conformity expression
Pseudo_n1^x / (Pseudo_n1^x + Pseudo_n2^x+ Pseudo_n3^x+ Pseudo_n4^x)

}


#3 Convex combination of conformity and experience

#Assign relative strength of conformity and experience bias; 1 is all conformity
  
k = 0.4
  
functConvex <- function(x){
  
  Ages1 <- 0
  for (i in Info["Ages",][Info["Choices",]==1]) {
    Ages1 <-  Ages1 + exp(b*i)
  }

  Ages2 <- 0
  for (i in Info["Ages",][Info["Choices",]==2]) {
    Ages2 <-  Ages2 + exp(b*i)
  }

  Ages3 <- 0
  for (i in Info["Ages",][Info["Choices",]==3]) {
    Ages3 <-  Ages3 + exp(b*i)
  }

  Ages4 <- 0
  for (i in Info["Ages",][Info["Choices",]==4]) {
    Ages4 <-  Ages4 + exp(b*i)
  }

  
     k * (n1^x / ((n1^x)+(n2^x)+(n3^x)+(n4^x))) +
  (1-k)* (Ages1 / (Ages1 + Ages2 + Ages3 + Ages4))
}


par(mfrow=c(1,3), 
    oma = c(2,2,0,0), 
    mar= c(2,2,3,2))

#plot implied probabilities

par(new=FALSE)
for (b in seq(-1,1, by=0.2)) {
  curve( functComb, from = 1, to=5, ylim=c(0,1),xlim=c(1,5.2),xlab= "", ylab = "")
  text(2+0.2, functComb(x=2), b)
  par(new=TRUE)
}
abline(h=n1/3, lty=2)
mtext(side=3, "Combined", line = 1.2)



par(new=FALSE)
for (b in seq(-1,1, by=0.2)) {
  curve( functPseudo, from = 1, to=5, ylim=c(0,1),xlim=c(1,5.2),xlab= "", ylab = "")
  text(2+0.2, functPseudo(x=2), b)
  par(new=TRUE)
}
abline(h=n1/3, lty=2)
mtext(side=3, "Pseudo-counts", line = 1.2)



par(new=FALSE)
for (b in seq(-1,1, by=0.2)) {
  curve( functConvex(x), from = 1, to=5, ylim=c(0,1),xlim=c(1,5.2),xlab= "", ylab = "")
  text(2+0.2, functConvex(x=2), b)
  par(new=TRUE)
}
abline(h=n1/3, lty=2)
mtext(side=3, paste("Convex Combination ( k =", k,")"), line = 1.2)

mtext(side=1, "Conformity exponent f", line = 1, outer = TRUE)
mtext(side=2, "Probability of choosing option 1", line = 1, outer = TRUE)


