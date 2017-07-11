
########################################################
################### Simulation PK 2 ####################
########################################################

library(aplpack)
library(deSolve)

plot.PK  <-  function(){
refresh.code  <-  function(...){
V1 <- slider(no=1); V2  <-  slider(no=2)
Q  <- slider(no=3); CL  <-  slider(no=4)
add.err <-  slider(no=5); cv  <-  slider(no=6)

parameters <- c(V1 = V1, V2 = V2, Q = Q, CL = CL)
state <- c(A1 = 100, A2 = 0)
Model	<-	function(t, state, parameters) {
 with(as.list(c(state, parameters)),{

 dA1 <- - state[1]*(CL/V1) - state[1]*(Q/V1) + state[2]*(Q/V2)
 dA2 <-   state[1]*(Q/V1) - state[2]*(Q/V2)

# return the rate of change
 list(c(dA1/V1, dA2/V2))
 }) # end with(as.list ...
 }

times <- seq(0, 24, length = 24)

# Write the simulation start and end number here !!
 simst <- 1
 simend <- 100 

 output <- NULL

 for (i in simst:simend) {
 set.seed(i)

 out <- ode(y = state, times = times, func = Model, parms = parameters)

# Combined error
 out[,2] <- rnorm(1:nrow(out), mean = out[,2], sd = add.err + out[,2]*cv)

 out.met  <- cbind(out, REP=i)
 output <- rbind(output, out.met)
 }

#plot(out, ylim=c(-10, 1100))
df.output <- as.data.frame(output)
avr.conc <- tapply(df.output$A1, df.output$time, mean)


 matplot(df.output$time, df.output[,2:3], type="n", xlab="Time, hour", ylab="Concentration"
         , xlim=range(c(-0.5, 24.5)), ylim=range(c(-1, max(df.output$A1, df.output$A2, na.rm=T)*1.05))
         , xaxs="i", yaxs="i", lab=c(40,8,18), cex=1, xaxt="n")
          
  for(i in unique(df.output$REP))
    {
      j <- df.output$REP==i
      axis(1, at=seq(0,24,4), labels=paste(seq(0,24,4)))
      lines(df.output$time[df.output$REP==i], df.output$A1[df.output$REP==i], col="red", lwd=2)
      points(df.output$time[df.output$REP==i], df.output$A2[df.output$REP==i], pch=1, col="blue", cex=0.9)
    }
    
grid()
}

slider(refresh.code, sl.names=c("V1", "V2", "Q", "CL", "Additive Error", "Proportional Error"),
sl.mins=c(0.1,0.1,0,0,0,0), sl.maxs=c(1000,1000,1000,1000,20,1), sl.deltas=c(0.1,0.1,0.1,.1, .1,.01), sl.defaults=c(10,20,10,2,1,0.1))
}
plot.PK()







########################################################
################ Simulation PK 3 (IIV) #################
########################################################

library(aplpack)
library(deSolve)

plot.PK  <-  function(){
refresh.code  <-  function(...){
TVV1  <- slider(no=1); TVV2  <-  slider(no=2)
TVQ   <- slider(no=3); TVCL  <-  slider(no=4)
OMV1  <- slider(no=5); OMV2  <-  slider(no=6)
OMQ   <- slider(no=7); OMCL  <-  slider(no=8)
add.err <-  slider(no=9); cv  <-  slider(no=10)

# Write the simulation start and end number here !!
 simst <- 1
 simend <- 100 

 output <- NULL

 for (i in simst:simend) {

 etaV1  <- rnorm(1,0, sqrt(OMV1))
 etaV2  <- rnorm(1,0, sqrt(OMV2))
 etaQ   <- rnorm(1,0, sqrt(OMQ))
 etaCL  <- rnorm(1,0, sqrt(OMCL))
 
 Model	<-	function(t, state, parameters) {
  with(as.list(c(state, parameters)),{

  dA1 <- - state[1]*(CL/V1) - state[1]*(Q/V1) + state[2]*(Q/V2)
  dA2 <-   state[1]*(Q/V1) - state[2]*(Q/V2)
 
 list(c(dA1/V1, dA2/V2))
 })
 }

 V1     <- TVV1*exp(etaV1)
 V2     <- TVV2*exp(etaV2)
 Q      <- TVQ*exp(etaQ)
 CL     <- TVCL*exp(etaCL)

parameters <- c(TVV1 = TVV1, TVV2 = TVV2, TVQ = TVQ, TVCL = TVCL, OMV1 = OMV1, OMV2 = OMV2, OMQ = OMQ, OMCL = OMCL)
state <- c(A1 = 500, A2 = 0)

times <- seq(0, 24, length = 24)

# set.seed(i)
 out <- ode(y = state, times = times, func = Model, parms = parameters)

 set.seed(i)
# Combined error
 out[,2] <- rnorm(1:nrow(out), mean = out[,2], sd = add.err + out[,2]*cv)

 out.met  <- cbind(out, REP=i)
 output <- rbind(output, out.met)
 }

df.output <- as.data.frame(output)

 matplot(df.output$time, df.output[,2:3], type="n", xlab="Time, hour", ylab="Concentration"
         , xlim=range(c(-0.5, 24.5)), ylim=range(c(-1, max(df.output$A1, df.output$A2, na.rm=T)*1.05))
         , xaxs="i", yaxs="i", lab=c(40,8,18), cex=1, xaxt="n")
         axis(1, at=seq(0,24,4), labels=paste(seq(0,24,4)))          
  for(i in unique(df.output$REP))
    {
      j <- df.output$REP==i
      lines(df.output$time[df.output$REP==i], df.output$A1[df.output$REP==i], col="red", lwd=2)
      points(df.output$time[df.output$REP==i], df.output$A2[df.output$REP==i], pch=1, col="blue", cex=0.9)
    }
      
      legend("topright", legend=c("Predicted Plasma Concentration", "Predicted Concentration in 2nd Compartment")
         , bty = "n", lty=c(1,NA), lwd=c(2,NA), pch=c(NA,16), col=c("red","blue"), cex=1)

grid()
}

slider(refresh.code, sl.names=c("TVV1", "TVV2", "TVQ", "TVCL", "OMV1", "OMV2", "OMQ", "OMCL", "Additive Error", "Proportional Error"),
sl.mins=c(0.1,0.1,0,0, 0.01,0.01,0.01,0.01, 0,0), sl.maxs=c(1000,1000,1000,1000,10,10,10,10,20,1), sl.deltas=c(0.1,0.1,0.1,0.1, 0.01,0.01,0.01,0.01, 0.1,0.01), sl.defaults=c(10,20,10,2, 0.1,0.1,0.1,0.1, 1,0.1))
}
plot.PK()








########################################################
################ Simulation PK 4 (IIV) #################
###################### POLYGON #########################
########################################################

library(aplpack)
library(deSolve)

plot.PK  <-  function(){
refresh.code  <-  function(...){
TVV1  <- slider(no=1); TVV2  <-  slider(no=2)
TVQ   <- slider(no=3); TVCL  <-  slider(no=4)
OMV1  <- slider(no=5); OMV2  <-  slider(no=6)
OMQ   <- slider(no=7); OMCL  <-  slider(no=8)
add.err <-  slider(no=9); cv  <-  slider(no=10)

# Write the simulation start and end number here !!
 simst <- 1
 simend <- 100 

 output <- NULL

 for (i in simst:simend) {

 etaV1  <- rnorm(1,0, sqrt(OMV1))
 etaV2  <- rnorm(1,0, sqrt(OMV2))
 etaQ   <- rnorm(1,0, sqrt(OMQ))
 etaCL  <- rnorm(1,0, sqrt(OMCL))
 
 Model	<-	function(t, state, parameters) {
  with(as.list(c(state, parameters)),{

  dA1 <- - state[1]*(CL/V1) - state[1]*(Q/V1) + state[2]*(Q/V2)
  dA2 <-   state[1]*(Q/V1) - state[2]*(Q/V2)
 
 list(c(dA1/V1, dA2/V2))
 })
 }


 V1     <- TVV1*exp(etaV1)
 V2     <- TVV2*exp(etaV2)
 Q      <- TVQ*exp(etaQ)
 CL     <- TVCL*exp(etaCL)

parameters <- c(TVV1 = TVV1, TVV2 = TVV2, TVQ = TVQ, TVCL = TVCL, OMV1 = OMV1, OMV2 = OMV2, OMQ = OMQ, OMCL = OMCL)
state <- c(A1 = 500, A2 = 0)

times <- seq(0, 24, length = 24)

# set.seed(i)
 out <- ode(y = state, times = times, func = Model, parms = parameters)

# Combined error
 out[,2] <- rnorm(1:nrow(out), mean = out[,2], sd = add.err + out[,2]*cv)

 out.met  <- cbind(out, REP=i)
 output <- rbind(output, out.met)
 }

df.output <- as.data.frame(output)

addShadedRegion <- function(x, y1, y2, col="#CC330030"){
polY<-c(y1,rev(y2))
polX<-c(x, rev(x))
polygon(polX, polY, col=col, border="lightgray", lty = 2)
}

addShadedRegion2 <- function(x, y1, y2, col="#0066FF30"){
polY<-c(y1,rev(y2))
polX<-c(x, rev(x))
polygon(polX, polY, col=col, border=col)
}


a95   <-  sapply(split(df.output$A1, df.output$time), quantile, probs=c(0.025, 0.5, 0.975))
b     <-  as.numeric(dimnames(a95)[[2]])


a95.2   <-  sapply(split(df.output$A2, df.output$time), quantile, probs=c(0.025, 0.5, 0.975))
b.2     <-  as.numeric(dimnames(a95.2)[[2]])


matplot(df.output$time, df.output[,2:3], type="n", xlab="Time, hour", ylab="Concentration"
         , xlim=range(c(-0.5, 24.5)), ylim=range(c(-1, max(df.output$A1, df.output$A2, na.rm=T)*1.05))
         , xaxs="i", yaxs="i", lab=c(40,8,18), cex=1, xaxt="n")
  axis(1, at=seq(0,24,4), labels=paste(seq(0,24,4)))
  grid()
  addShadedRegion(b, a95[1,], a95[3,])
  addShadedRegion2(b.2, a95.2[1,], a95.2[3,])
  lines(b, a95[2,], col="red", lwd=2)
  lines(b.2, a95.2[2,], col="darkblue", lwd=2)
  legend("topright", legend=c("Median Predicted Plasma Concentration", "95% Prediction Interval for Predicted Plasma Concentration"
        , "Median Predicted Concentration in 2nd Compartment", "95% Prediction Interval for Predicted Concentration in 2nd Compartment")
        , pch=c(16,NA,16,NA), col=c("red","#CC3300", "darkblue","#0066FF"), lwd=c(2,16,2,16), cex=0.9, box.col = "white")
}

slider(refresh.code, sl.names=c("TVV1", "TVV2", "TVQ", "TVCL", "OMV1", "OMV2", "OMQ", "OMCL", "Additive Error", "Proportional Error"),
sl.mins=c(0.1,0.1,0,0, 0.01,0.01,0.01,0.01, 0,0), sl.maxs=c(1000,1000,1000,1000,10,10,10,10,20,1), sl.deltas=c(0.1,0.1,0.1,0.1, 0.01,0.01,0.01,0.01, 0.1,0.01), sl.defaults=c(10,20,10,2, 0.1,0.1,0.1,0.1, 1,0.1))
}
plot.PK()













########################################################
################ Simulation PK 5 (IIV) #################
###################### POLYGON #########################
################### IIV COVARIANCE #####################
########################################################

library(aplpack)
library(deSolve)
library(mvtnorm)

plot.PK  <-  function(){
refresh.code  <-  function(...){
TVV1  <- slider(no=1); TVV2  <-  slider(no=2)
TVQ   <- slider(no=3); TVCL  <-  slider(no=4)
OMV1  <- slider(no=5); OMV2  <-  slider(no=6)
OMQ   <- slider(no=7); OMCL  <-  slider(no=8)
add.err <-  slider(no=9); cv  <-  slider(no=10)


############################################################################
############################################################################

###################### Use correlation value as input ######################

#corV1V2 <- 0.2
#corV1Q  <- 0.2
#corV1CL <- 0.2
#corV2Q  <- 0.2
#corV2CL <- 0.2
#corQCL  <- 0.2

#covV1V2 <- corV1V2*sqrt(OMV1)*sqrt(OMV2)   # covariance between V1 and V2
#covV1Q  <- corV1V2*sqrt(OMV1)*sqrt(OMV2)   # covariance between V1 and V2
#covV1CL <- corV1V2*sqrt(OMV1)*sqrt(OMV2)   # covariance between V1 and V2
#covV2Q  <- corV1V2*sqrt(OMV1)*sqrt(OMV2)   # covariance between V1 and V2
#covV2CL <- corV1V2*sqrt(OMV1)*sqrt(OMV2)   # covariance between V1 and V2
#covQCL  <- corV1V2*sqrt(OMV1)*sqrt(OMV2)   # covariance between V1 and V2

###################### Use correlation value as input ######################

###################### Use covariance value as input #######################

covV1V2 <- 0.02                             # covariance between V1 and V2
covV1Q  <- 0.03                             # covariance between V1 and Q
covV1CL <- 0.01                             # covariance between V1 and CL
covV2Q  <- 0.03                             # covariance between V2 and Q
covV2CL <- 0.02                             # covariance between V2 and CL
covQCL  <- 0.04                             # covariance between Q and CL

###################### Use covariance value as input #######################

LOGTHETA <- c(log(TVV1), log(TVV2), log(TVQ), log(TVCL))

OMEGA <- matrix(0,4,4)
OMEGA[1,1] <- OMV1
OMEGA[2,2] <- OMV2
OMEGA[3,3] <- OMQ
OMEGA[4,4] <- OMCL
OMEGA[1,2] <- covV1V2
OMEGA[2,1] <- covV1V2
OMEGA[1,3] <- covV1Q
OMEGA[3,1] <- covV1Q
OMEGA[1,4] <- covV1CL
OMEGA[4,1] <- covV1CL
OMEGA[1,4] <- covV1CL
OMEGA[4,1] <- covV1CL
OMEGA[2,3] <- covV2Q
OMEGA[3,2] <- covV2Q
OMEGA[2,4] <- covV2CL
OMEGA[4,2] <- covV2CL
OMEGA[3,4] <- covQCL
OMEGA[4,3] <- covQCL

N <- 300                                  # of Replicates in the Simulation
set.seed(20170523)
LOGPARA  <- rmvnorm(N, LOGTHETA, OMEGA)
PARA     <- exp(LOGPARA)

PARA     <- as.data.frame(PARA)
names(PARA) <- c("V1", "V2", "Q", "CL")

############################################################################
############################################################################


# Write the simulation start and end number here !!
 simst <- 1
 simend <- N

 output <- NULL

 for (i in simst:simend) {

 V1     <- PARA$V1[i]
 V2     <- PARA$V2[i]
 Q      <- PARA$Q[i]
 CL     <- PARA$CL[i]
 
 Model	<-	function(t, state, parameters) {
  with(as.list(c(state, parameters)),{

  dA1 <- - state[1]*(CL/V1) - state[1]*(Q/V1) + state[2]*(Q/V2)
  dA2 <-   state[1]*(Q/V1) - state[2]*(Q/V2)
 
 list(c(dA1/V1, dA2/V2))
 })
 }

parameters <- c(TVV1 = TVV1, TVV2 = TVV2, TVQ = TVQ, TVCL = TVCL, OMV1 = OMV1, OMV2 = OMV2, OMQ = OMQ, OMCL = OMCL)
state <- c(A1 = 500, A2 = 0)

times <- seq(0, 24, length = 24)

# set.seed(i)
 out <- ode(y = state, times = times, func = Model, parms = parameters)

 set.seed(i)
# Combined error
 out[,2] <- rnorm(1:nrow(out), mean = out[,2], sd = add.err + out[,2]*cv)

 out.met  <- cbind(out, REP=i)
 output <- rbind(output, out.met)
 }


df.output <- as.data.frame(output)

addShadedRegion <- function(x, y1, y2, col="#CC330030"){
polY<-c(y1,rev(y2))
polX<-c(x, rev(x))
polygon(polX, polY, col=col, border="lightgray", lty = 2)
}

addShadedRegion2 <- function(x, y1, y2, col="#0066FF30"){
polY<-c(y1,rev(y2))
polX<-c(x, rev(x))
polygon(polX, polY, col=col, border=col)
}


a95   <-  sapply(split(df.output$A1, df.output$time), quantile, probs=c(0.025, 0.5, 0.975))
b     <-  as.numeric(dimnames(a95)[[2]])
#c     <-  data.frame(TIME=b,"2.5P"=a95[1,],"50P"=a95[2,],"97.5P"=a95[3,])
#colnames(c) <-  c("TIME", "2.5P", "50P", "97.5P")


a95.2   <-  sapply(split(df.output$A2, df.output$time), quantile, probs=c(0.025, 0.5, 0.975))
b.2     <-  as.numeric(dimnames(a95.2)[[2]])
#c.2     <-  data.frame(TIME=b.2,"2.5P"=a95.2[1,],"50P"=a95.2[2,],"97.5P"=a95.2[3,])
#colnames(c.2) <-  c("TIME", "2.5P", "50P", "97.5P")


matplot(df.output$time, df.output[,2:3], type="n", xlab="Time, hour", ylab="Concentration"
         , xlim=range(c(-0.5, 24.5)), ylim=range(c(-1, max(df.output$A1, df.output$A2, na.rm=T)*1.05))
         , xaxs="i", yaxs="i", lab=c(40,8,18), cex=1, xaxt="n")
  axis(1, at=seq(0,24,4), labels=paste(seq(0,24,4)))
  grid()
  addShadedRegion(b, a95[1,], a95[3,])
  addShadedRegion2(b.2, a95.2[1,], a95.2[3,])
  lines(b, a95[2,], col="red", lwd=2)
  lines(b.2, a95.2[2,], col="darkblue", lwd=2)
  legend("topright", legend=c("Median Predicted Plasma Concentration", "95% Prediction Interval for Predicted Plasma Concentration"
        , "Median Predicted Concentration in 2nd Compartment", "95% Prediction Interval for Predicted Concentration in 2nd Compartment")
        , pch=c(16,NA,16,NA), col=c("red","#CC3300", "darkblue","#0066FF"), lwd=c(2,16,2,16), cex=0.9, box.col = "white")
}

slider(refresh.code, sl.names=c("TVV1", "TVV2", "TVQ", "TVCL", "OMV1", "OMV2", "OMQ", "OMCL", "Additive Error", "Proportional Error"),
sl.mins=c(0.1,0.1,0,0, 0.01,0.01,0.01,0.01, 0,0), sl.maxs=c(1000,1000,1000,1000,10,10,10,10,20,1), sl.deltas=c(0.1,0.1,0.1,0.1, 0.01,0.01,0.01,0.01, 0.1,0.01), sl.defaults=c(10,20,10,2, 0.1,0.1,0.1,0.1, 1,0.1))
}
plot.PK()


