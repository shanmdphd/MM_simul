################################################################################
###############   Sequential Running, In Silico M-M Experiment   ###############
################################################################################

WD <- "./est"   # Directory containing the orinigal nonmem code

DD <- "./"       # Directory containing the orinigal nonmem dataset


setwd(paste(WD))

codename <- "MM"   					  # Put the name of the code without extension (ctl)

dataname <- "MM_sim_data.csv"   		  

nmif <- paste(codename,".ctl", sep="") 			  # NONMEM input file

nmof <- paste(codename,".out", sep="") 			  # NONMEM output file

# make new directory for the multiple estimation

nwd <- paste(WD, paste("mest_", codename,sep=""),sep="")     # dx directory

dir.create(nwd)

O.data <- read.csv(paste(DD,dataname, sep=""), na.strings=".", stringsAsFactors = F)     # Original dataset


################################ Multiple Runs #################################

  setwd(paste(nwd))                                            

  O.code <- readLines(paste(WD, nmif,sep=""))

  O.code <- O.code[-(min(grep("TABLE ", O.code)):(max(grep("TABLE ", O.code))+1))]    # Remove Tables from the original code

  O.code <- O.code[-grep("COV", O.code)]    # Remove Tables from the original code    # Remove Covariance step from the original code

  mruncode <- NULL

  mrunraw  <- NULL       # Boot raw data

 for (i in unique(O.data$REP)) {

   cat(i,"\n")

   out <- NULL

   outraw <- NULL

   mruncode <- sub(O.code[grep("DATA", O.code)], paste(O.code[grep("DATA", O.code)]," ACCEPT(REP.EQ.",i,")", sep=""), O.code, fixed = TRUE)

   write(mruncode, file="mrun.ctl" , sep="")
   
   mrun <- paste("nmfe73", "mrun.ctl", paste(paste("mrun",i,sep=""),".out", sep=""), sep=" ")	 # NONMEM 7.3

   system(mrun, invisible=T, show.output.on.console=F, wait = TRUE)

   out <- read.table(paste("mrun",".ext",sep=""), header=T, skip=1)

   outraw <- out[which(out$ITERATION==-1000000000),-1]

   OUT <- readLines(paste("mrun",i,".out",sep=""))
  
   outraw$SUCCESS <- as.numeric("0MINIMIZATION SUCCESSFUL" %in% OUT)
 
   outraw$ROUNDING <- as.numeric(" DUE TO ROUNDING ERRORS (ERROR=134)" %in% OUT)
  
   outraw$REPN <- i
 
   mrunraw <- rbind(mrunraw, outraw)
 

  }

  write.csv(mrunraw, "mrunraw.csv", row.names=F)


