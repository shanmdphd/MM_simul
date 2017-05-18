
setwd("./")
# setwd("./shanmdphd/MM_simul")
######################################################
################## Simulation Plot ###################
######################################################

 library(aplpack)
 library(deSolve)

VM <- 100
KM <- 500

plot.MM  <-  function(){
refresh.code  <-  function(...){
VM <- slider(no=1); KM  <-  slider(no=2)
add.err <-  slider(no=3); cv  <-  slider(no=4)

parameters <- c(VM = VM, KM = KM)
state <- c(S = 1000)
Model	<-	function(t, state, parameters) {
 with(as.list(c(state, parameters)),{
# rate of change
 dS <- - VM*S/(KM+S)

# return the rate of change
 list(c(dS))
 }) # end with(as.list ...
 }

times <- seq(0, 24, by = 1)

# Write the simulation start and end number here !!
 simst <- 1
 simend <- 1000 

 output <- NULL

 for (i in simst:simend) {
 set.seed(i)

 out <- ode(y = state, times = times, func = Model, parms = parameters)

# Additive error
#add.err <- 5
#out[,2] <- rnorm(1:nrow(out), mean = out[,2], sd = add.err)

# Proportional error
#cv    <- 0.10
#out[,2] <- rnorm(1:nrow(out), mean = out[,2], sd = out[,2]*cv)

# Combined error
 out[,2] <- rnorm(1:nrow(out), mean = out[,2], sd = add.err + out[,2]*cv)

 out.met  <- cbind(out, REP=i)
 output <- rbind(output, out.met)
 }

#plot(out, ylim=c(-10, 1100))
df.output <- as.data.frame(output)
avr.conc <- tapply(df.output$S, df.output$time, mean)

plot(x=output[,1], y=output[,2], xlab="Time", ylab="Concentration", ylim=c(-10, 1600), lty=1)
lines(as.numeric(names(avr.conc)), avr.conc, col="red", lwd=2)
grid()
}
slider(refresh.code, sl.names=c("Vmax", "Km", "Additive Error", "Proportional Error"),
sl.mins=c(0,0,0,0), sl.maxs=c(1000,10000,300,1), sl.deltas=c(.1,.1, .1,.01), sl.defaults=c(100,500,5,0.1))
}
plot.MM()





######################################################
################## Individual Plot ###################
######################################################

library(aplpack)
library(deSolve)

VM <- 100
KM <- 500

plot.MM  <-  function(){
refresh.code  <-  function(...){
VM <- slider(no=1); KM  <-  slider(no=2)
add.err <-  slider(no=3); cv  <-  slider(no=4)

parameters <- c(VM = VM, KM = KM)
state <- c(S = 1000)
Model	<-	function(t, state, parameters) {
 with(as.list(c(state, parameters)),{
# rate of change
 dS <- - VM*S/(KM+S)

# return the rate of change
 list(c(dS))
 }) # end with(as.list ...
 }

times <- seq(0, 24, by = 1)

 out <- ode(y = state, times = times, func = Model, parms = parameters)

# Additive error
#add.err <- 5
#out[,2] <- rnorm(1:nrow(out), mean = out[,2], sd = add.err)

# Proportional error
#cv    <- 0.10
#out[,2] <- rnorm(1:nrow(out), mean = out[,2], sd = out[,2]*cv)

# Combined error
 out[,2] <- rnorm(1:nrow(out), mean = out[,2], sd = add.err + out[,2]*cv)
 plot(out, xlab="Time", ylab="Concentration", ylim=c(-10, 1100), lty=1)
 grid()
 }
slider(refresh.code, sl.names=c("Vmax", "Km", "Additive Error", "Proportional Error"),
sl.mins=c(0,0,0,0), sl.maxs=c(1000,10000,300,1), sl.deltas=c(.1, .1, .1, .01), sl.defaults=c(100,500,5,0.1))
 }
plot.MM()





######################################################
############ Generate Simulation Dataset #############
######################################################

library(aplpack)
library(deSolve)

parameters <- c(VM = 100, KM = 500)
state <- c(S = 1000)

Model	<-	function(t, state, parameters) {
 with(as.list(c(state, parameters)),{
# rate of change
 dS <- - VM*S/(KM+S)

# return the rate of change
 list(c(dS))
 }) # end with(as.list ...
}

times <- seq(0, 24, by = 1)

# Write the simulation start and end number here !!
 simst <- 1
 simend <- 1000 

 output <- NULL

 for (i in simst:simend) {
 set.seed(i)

 out <- ode(y = state, times = times, func = Model, parms = parameters)

# Additive error
 add.err <- 5
# out[,2] <- rnorm(1:nrow(out), mean = out[,2], sd = add.err)

# Proportional error
 cv    <- 0.10
# out[,2] <- rnorm(1:nrow(out), mean = out[,2], sd = out[,2]*cv)

# Combined error
 out[,2]  <- rnorm(1:nrow(out), mean = out[,2], sd = add.err + out[,2]*cv)

 out.met  <- cbind(out, REP=i)
 output <- rbind(output, out.met)
}

output <- as.data.frame(output)
output$ID <- 101

F.output <- output[,c(4,1,2,3)]

write.table(F.output, file="MM_sim_data.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=F)


