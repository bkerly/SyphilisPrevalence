# From supplement from
# Spicknall, Ian H. PhD; Kreisel, Kristen M. PhD; Weinstock, Hillard S. MD. Estimates of the Prevalence and Incidence of Syphilis in the United States, 2018. Sexually Transmitted Diseases 48(4):p 247-252, April 2021. | DOI: 10.1097/OLQ.0000000000001364

library(deSolve)
args = commandArgs(trailingOnly=TRUE)
SEED<- as.numeric(args[1])

SP<- as.numeric(args[2]) #sub popID num
SAMPS<- as.numeric(args[3])
FILENAME<- paste("d_", SP,"_", SEED, ".RDS", sep="")

setwd("/scicomp/home/xfu0/syphEst/")
inpD<- read.csv("syph_inp_2008.csv")

setwd(paste("/scicomp/home/xfu0/syphEst/data/", sep=""))

###x = the draw, i= the subpopNum
objFcn<- function(x, i=1){
  times <- seq(2017,2018, by=1)#seq(0, to=1000, by=2)
  p<- getParams(x=x, inp=inpD, i=i)
  states<- c( U= p$initN -p$initI1*p$initN -p$initI2*p$initN, 
              I1=p$initI1*p$initN,  
              I2=p$initI2*p$initN, 
              CumInf=0, CumK1=0, CumK2=0,
              initI1<- p$initI1*p$initN,
              initI2<- p$initI2*p$initN
  )
  sim<- as.data.frame(ode(states, times, mod.fun, p ))[2,]
  val    <- c(sim$CumK1[length(sim$CumK1)],sim$CumK2[length(sim$CumK2)])
  target <- c(inpD$K_E[i],                 inpD$K_L[i])
  ret = list()
  ret$sim= sim
  ret$dist = sum(abs(val-target))#1-sum(abs(val-target))/1.851172
  ret$val = val
  ret
}

mod.fun = function(t, state, parameters){
  with(as.list(c(state, parameters)),{
    dU<-  mu-(omega+lambda)*U +I1*(sigma+tau1) + I2*(sigma+tau2) 
    dI1<-   lambda*U -(omega+sigma+tau1+alpha)*I1 
    dI2<-   alpha*I1 -(omega+sigma+tau2)*I2 
    dCumInf<- lambda*U
    dCumK1<- rho*I1*(sigma+tau1)
    dCumK2<- rho*I2*(sigma+tau2)
    dInitI1<- 0
    dInitI2<- 0
    dxdt<- c( dU, dI1, dI2, dCumInf, dCumK1, dCumK2, dInitI1, dInitI2)
    list(dxdt)  
  })
}
#x=the draw, inp=the inputCSV data, i=the subPopNum
getParams<- function(x, inp, i){
  p<- data.frame(mu=.1)  
  p$omega <-0
  p$mu    <- 0
  p$alpha <- 1
  p$rho   <- 1
  p$initN <- inp$N0[i]
  p$lambda<- (x[1]*(inp$lambda_max[i]^.1- inp$lambda_min[i]^.1)+inp$lambda_min[i]^.1)^10
  p$initI1<- (x[2]*(inp$E_0_max[i]^.1- inp$E_0_min[i]^.1)+inp$E_0_min[i]^.1)^10
  p$initI2<- (x[3]*(inp$L_0_max[i]^.1- inp$L_0_min[i]^.1)+inp$L_0_min[i]^.1)^10
  p$sigma <- (x[4]*(inp$sigma_max[i]^.1- inp$sigma_min[i]^.1)+inp$sigma_min[i]^.1)^10
  p$tau1  <- (x[5]*(inp$tau1_max[i]^.1- inp$tau1_min[i]^.1)+inp$tau1_min[i]^.1)^10
  p$tau2  <- (x[6]*(inp$tau2_max[i]^.1- inp$tau2_min[i]^.1)+inp$tau2_min[i]^.1)^10
  return(p)
}

####################do the while loop to do the fitting
logistic = function(x)1/(1+exp(-x))
logit = function(p)log(p/(1-p))
sd<- 0.035
numToRestart<- 50000
nsamp=SAMPS# stop after nsamp steps
e=10 # error 
####The output object is declared below
ret = matrix(nrow=nsamp, ncol=1+6+1+9)
#1 = subpopID  
#6=number of params varied
#1 = error
#9 = modelSim vector length
p=runif(6, min=0, max=1)
a = objFcn(p, SP)
cur.l = a$dist
cur.p = logit(p)
cur.sim<- a$sim
colnames(ret)<-  c("Population", "lambdaP", "initE0_P", "initL0_P", "sigma_P", "tau1_P", "tau2_P", 
                   "error",names(cur.sim))
orig_err<- cur.l
running = TRUE
ind = count = 0
tries<- 0
while(running){
  while(tries> numToRestart){
    print(paste("RESTART HILL CLIMBING", orig_err, cur.l, logistic(cur.p[1]), logistic(cur.p[2]), logistic(cur.p[3]),logistic(cur.p[4]), logistic(cur.p[5]), logistic(cur.p[6])    ))
    tries<- 0
    p=runif(6, min=0, max=1)
    a = objFcn(p, SP)
    cur.l = a$dist
    cur.p = logit(p)
    cur.sim<- a$sim
    orig_err<- cur.l
  }
  tries<- tries+1
  tmp.p = cur.p + rnorm(6, 0, sd)# 6 = number of parameters being varied
  tmp.p[tmp.p > 9] = 9; tmp.p[tmp.p < -9] = -9 #Used to avoid sink at extreem values, equivalent of Unif(logistic(-9), logisitc(9)) prior
  a = objFcn(logistic(tmp.p), SP)
  if((a$dist < cur.l) | (a$dist < e)){ #do hill climbing
    cur.p = tmp.p
    cur.l = a$dist
    cur.sim= a$sim
    tries<- 0
    if(a$dist < e & (cur.sim$I1+cur.sim$I2) >= (cur.sim$V8+cur.sim$V9)){ 
      ind = ind + 1
      ret[ind,] = as.numeric(c(SP, as.numeric(logistic(cur.p)), cur.l, cur.sim))
      if(ind >= nsamp)running = F
    }
    if(a$dist<e & (cur.sim$I1+cur.sim$I2)< (cur.sim$V8+cur.sim$V9)){
      tries= 1+numToRestart
    }
  }
}
saveRDS(ret, file=FILENAME)
