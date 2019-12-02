#rm(list = ls())

####################################################################
## Author: 	Dan Lythgoe								##
## Date:	2019-11-27								##
## Description: Example of how to fit a latent class model 		##
##			with a time-to-event distal outcome variable	##
##			using the LCSM() R function. Example does 	##
##			not match those in the thesis and instead 	##
##			serves as as a simple of example.			##
####################################################################

##############
## Packages ##
##############

require(poLCA)
require(survival)
require(mixtools)
require(ggplot2)
require(grid)
require(gridExtra)
require(boot)
require(mclust)
require(matrixcalc)
library(xtable)


###################################
## Need to specify file location ##
###################################

my.loc <- "C:\\Users\\knmf5025\\Dropbox\\PhD\\LCA project\\programs\\"		##<--specify your file location here
source(paste(my.loc,"github_lcsm_functions_v40.r",sep=""))


###################################################################################
## Download data set from http://lib.stat.cmu.edu/datasets/Andrews/ - Table 46.1 ##
###################################################################################

byar <- read.csv(paste(my.loc,"prostate.csv",sep=""),header=TRUE)


##############################
## Some data pre processing ##
##############################

##Sort data set by survival time and add 1 to account for subjects with death time = 0
byar2 <- byar[order(byar$dtime),]
byar2$dtime <- byar2$dtime + 1


##Create new survival status ('died') variable
byar2$died <- rep(NA,length(byar$patno))
  byar2$died[grepl("dead",byar$status)] <- 1
  byar2$died[grepl("alive",byar$status)] <- 0


##Categorise treatments into Untreated (placebo, 0.2) and Treated (1.0 and 5.0)
byar2$trt <- rep(NA,length(byar2$patno))
  byar2$trt[byar2$rx=="placebo" | byar2$rx=="0.2 mg estrogen"] <- 0
  byar2$trt[byar2$rx=="1.0 mg estrogen" | byar2$rx=="5.0 mg estrogen"] <- 1


##Create the y matrix (categorical manifest variables)
byar2.y <- cbind( 
            model.matrix(~-1+factor(byar2$hx)), 
            model.matrix(~-1+factor(byar2$bm))
            )


##Create the w matrix (continuous manifest variables)
byar2.w <- cbind(log(byar2$ap))
colnames(byar2.w) <- "logap"


####################
## Fit LCSM model ##
####################

##First create time and indicator matrices which are used to fit a piecewise exponential model with two partitions either side of 30 months
ind.mat.list1 <- ind.mat.fn12(U=byar2$dtime,delta=byar2$died,part=c(0,30,max(byar2$dtime+1)))


#rm(lcsm.mod)
lcsm.mod <- LCSM(y=byar2.y,
  x=NULL,
  z=byar2$trt,
  t=byar2$dtime,
  delta=byar2$died,
  nclass=2,
  nreps=3,
  w=byar2.w,
  alpha.method="PEM",
  ind.mat=ind.mat.list1$ind.mat,
  t.mat=ind.mat.list1$t.mat,
  y.var.vec=c(1,1,2,2)
  )
lcsm.mod$params





