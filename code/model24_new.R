# load the cleaned dataset and required packages --------------------------
library("nimble")
library(readxl)

data_cleaned = read.csv("data_input.csv")
for (i in 1:ncol(data_cleaned)){
  data_cleaned[,i] = as.numeric(data_cleaned[,i])
}
Tage = max(data_cleaned[,3])
Ttime = max(data_cleaned[,5] + data_cleaned[,3])
tbdata_adj <- read_excel("../data/tbdata_adj.xlsx")
# number of spatial locations
N = length(tbdata_adj$Adj_ID)
adj_info = read.csv("../data/Adj.txt",header = F)
adj_num = read.csv("../data/Num.txt",header = F)[,1]
adj = rep(0, sum(adj_num))
count = 1
for (i in 1:nrow(adj_info)){
  for (j in 1:ncol(adj_info)){
    if (!is.na(adj_info[i,j])){
      adj[count] = adj_info[i,j]
      count = count+1
    }
  }
}
weights = rep(1, sum(adj_num))
num = adj_num


# optional, select a subset of the data to do the analysis
# data_cleaned = data_cleaned[1:100,]



# Set initial values for the model ----------------------------------------
initialize_variables <- function(){
  mu = runif(1, -3, -1)
  mbeta0 = 0
  fbeta0 = runif(1,-1,1)
  mcar_age = rnorm(Tage, sd=0.1)
  fcar_age = rnorm(Tage, sd=0.1)
  mhprec <- rgamma(1, 2, 1)
  fhprec <- rgamma(1, 2, 1)
  mhetero = rnorm(N, sd=0.1)
  fhetero = rnorm(N, sd=0.1)
  mprec = rgamma(1, 2, 1)
  fprec = rgamma(1, 2, 1)
  u = matrix(0, nrow = 2, ncol = N)
  omega = matrix(c(1,0,0,1), 2, 2)
  tmprec = rgamma(1,0.01, 0.01)  
  tfprec = rgamma(1,0.01, 0.01)
  mcar_time = rnorm(Ttime, sd=0.1)
  fcar_time = rnorm(Ttime, sd=0.1)
  
  # set to zero constraint
  mcar_age[1] = 0
  fcar_age[1] = 0
  mhetero[1] = 0
  fhetero[1] = 0
  mcar_time[1] = 0
  fcar_time[1] = 0
  
  return(list(mu = mu, mbeta0 = mbeta0, fbeta0 = fbeta0, mcar_age = mcar_age, fcar_age = fcar_age, 
              mhprec = mhprec, fhprec = fhprec, mhetero = mhetero, fhetero = fhetero,
              omega=omega, mprec = mprec, fprec = fprec, u = u, tmprec = tmprec, tfprec = tfprec,
              mcar_time = mcar_time, fcar_time = fcar_time))
}

# To avoid the term in the exponential function explode, we need to check if the initial value is appropriate
check_initial_condition <- function(){
  m = nrow(data_cleaned)
  mu = runif(1,-3, -1)
  mbeta0 = 0
  fbeta0 = runif(1,-1,1)
  mcar_age = rnorm(Tage, sd=0.1)
  fcar_age = rnorm(Tage, sd=0.1)
  mhetero = rnorm(N, sd=0.1)
  fhetero = rnorm(N,sd=0.1)
  mvspace = matrix(0,nrow = 2, ncol = N)
  mcar_time = rnorm(Ttime, sd = 0.1)
  fcar_time = rnorm(Ttime, sd = 0.1)
  
  mcar_age[1] = 0
  fcar_age[1] = 0
  mhetero[1] = 0
  fhetero[1] = 0
  mcar_time[1] = 0
  fcar_time[1] = 0
  data = data_cleaned
  p = rep(NA,m)
  dayhaz = matrix(0,nrow = m, ncol = Tage)
  for (j in 1:m) {
    for (k in 1:data[j,3]) {
      dayhaz[j,k]<-1/2* exp(mu + (1-data[j,2])*(fbeta0+fcar_age[k]+fcar_time[data[j,5]+k])+data[j,2]*(mbeta0+mcar_age[k]+mcar_time[data[j,5]+k])+mvspace[1,data[j,1]]*data[j,2]+mvspace[2,data[j,1]]*(1-data[j,2])+
                              mhetero[data[j,1]]*data[j,2]+fhetero[data[j,1]]*(1-data[j,2]))
    } 
    #print(sum(dayhaz[j,1:data[j,3]]))
    #print(data[j,2])
    p[j] <-1-exp(-sum(dayhaz[j,1:data[j,3]]))
  } 
  return(p)
}

check_initial_condition()

# Main model --------------------------------------------------------------
model <- nimbleCode( { 
  # grand mean
  #mu ~ dflat()
  mu ~ dnorm(-2, 1)
  
  # sex effect set to zero
  mbeta0 ~ dnorm(0,100000)
  #fbeta0 ~ dflat()
  fbeta0 ~ dnorm(0, 1)
  
  # MYBYM model on space effect
  # spatially independent part
  # under this prior, mhprec and fhprec will have a large probability around 0
  #mhprec ~ dgamma(0.5, 0.0025)   # these numbers are directly copied from the paper
  #fhprec ~ dgamma(0.5, 0.0025)   # these numbers are directly copied from the paper
  mhprec ~ dgamma(2, 1)
  fhprec ~ dgamma(2, 1)
  
  # apply set to zero constraint, set the first treatment effect to 0
  # make the variance extremelly small to achieve this goal
  mhetero[1] ~ dnorm(0,10000)
  fhetero[1] ~ dnorm(0,10000)
  for (i in 2:N){
    mhetero[i] ~ dnorm(0, mhprec)
    fhetero[i] ~ dnorm(0, fhprec)
  }
  
  # spatially dependent part
  omega[1:2,1:2] ~ dwish(R[1:2,1:2],2)    
  Cov[1:2,1:2] <- inverse(omega[1:2,1:2])
  achol[1:2,1:2] <- t(chol(Cov[1:2,1:2]))
  #cor12 <- Cov[1,2]/(sqrt(Cov[1,1])*sqrt(Cov[2,2]))
  for (k in 1:2){
    sprec[k] <- 1
    # set the parameter zero_mean = 1 to apply the sum to zero constraint
    u[k,1:N] ~ dcar_normal(adj[1:length_adj], weights[1:length_adj], num[1:N], sprec[k], zero_mean = 1)
  }
  for (i in 1:N){  
    mvspace_temp[1:2,i] <- achol[1:2,1:2]%*%u[1:2,i] 
  }
  
  mvspace[1, 1:N] <- mvspace_temp[1, 1:N] / sd(mvspace_temp[1, 1:N])
  mvspace[2, 1:N] <- mvspace_temp[2, 1:N] / sd(mvspace_temp[2, 1:N])
  
  
  
  
  
  # RW1 model on age effect
  #mprec~dgamma(0.1,0.1)
  #fprec~dgamma(0.1,0.1)
  mprec ~ dgamma(2, 1)
  fprec ~ dgamma(2, 1)
  
  # apply set to zero constraint, set the first treatment effect to 0
  mcar_age[1]~dnorm(0,10000)  
  fcar_age[1]~dnorm(0,10000)
  for (i in 2:Tage){
    mcar_age[i] ~ dnorm(mcar_age[i-1],mprec)
    fcar_age[i] ~ dnorm(fcar_age[i-1],fprec)
  }
  
  # time effect
  #tmprec  ~ dgamma(0.5,0.0025)  
  #tfprec  ~ dgamma(0.5,0.0025)  
  tmprec ~ dgamma(2, 1)
  tfprec ~ dgamma(2, 1)
  mcar_time[1]~dnorm(0,10000) 
  fcar_time[1]~dnorm(0,10000) 
  for (i in 2:Ttime) {
    mcar_time[i]~dnorm(mcar_time[i-1],tmprec)
    fcar_time[i]~dnorm(fcar_time[i-1],tfprec)
  }
  
  # formulation 
  for (j in 1:m) {
    for (k in 1:data[j,3]) {
      gamma[j,k]<- mu + (1-data[j,2])*(fbeta0+fcar_age[k]+fcar_time[data[j,5]+k])+data[j,2]*(mbeta0+mcar_age[k]+mcar_time[data[j,5]+k])+mvspace[1,data[j,1]]*data[j,2]+mvspace[2,data[j,1]]*(1-data[j,2])+
        mhetero[data[j,1]]*data[j,2]+fhetero[data[j,1]]*(1-data[j,2])
      dayhaz[j,k]<-1/2*exp(gamma[j,k]) # since the interval is 0.5 year
    } 
    pos[j]~dbern(p[j])
    icumhaz[j]<-sum(dayhaz[j,1:data[j,3]])
    p[j] <-1-exp(-icumhaz[j])
  } 
}
)

nimble_data <- list(
  pos = data_cleaned[,4]
)

nimble_constant <- list(
  Tage = Tage,
  data = data_cleaned,
  N = N,
  m = nrow(data_cleaned),
  adj = adj,
  weights = weights,
  num = num,
  length_adj = length(adj),
  R = matrix(c(0.015, 0, 0, 0.2),2, 2)
)

inits <- function(){
  return(initialize_variables())
}

#params <- c("mbeta0","fbeta0","mcar_age","fcar_age","omega","mvspace","cor12")

samples <- nimbleMCMC(
  code = model,
  data = nimble_data,
  constants = nimble_constant,
  inits = inits,
  nchains = 1,
  #monitors = params,
  niter = 25000,
  nburnin = 5000,
  thin = 2,
  summary = TRUE,
  WAIC = TRUE)
print(samples$WAIC)

# save the result
save(samples,file = "model24.RData")


