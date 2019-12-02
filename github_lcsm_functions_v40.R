#rm(list = ls())

####################################################################
## Author: 	Dan Lythgoe					  ##
## Date:	2019-11-27					  ##
## Description: LCSM function and subfunctions  		##
####################################################################


####################
## LCSM functions ##
####################

usetest <- FALSE
if(usetest != FALSE){
  X <- X4
  X2
  beta <- beta.true
  v <- v.true
  #v2 <- 0.5
  kappa.mat <- kappa.mat.true
  alpha.start <- rep(0.01,n)
  p
  delta
 #des <- cbind(1,x)
  nclass <- 2
  sigma.sq.true.matrix <- matrix(0,length(sigma.sq.true),length(sigma.sq.true))
  diag(sigma.sq.true.matrix) <- sigma.sq.true
  mu <- mu.true
  alpha.method <- "PEM"
  #U2 <- log(U)
  #v2 <- NULL
  #GAMMA.true <- 0.05
  alpha.PEM.start <- rep(0.01,ncol(t.mat))
}

##############
## Packages ##
##############

require(matrixcalc) #required for test of is.singular.matrix
#require(mixtools)
#require(mclust)
require(mvtnorm)

###########################################################################
post.fn12 <- function(kappa.mat,alpha,p,beta,v,des,X,z,delta,nclass,lambda1=1,lambda2=1,X2=NULL,mu=NULL,sigma.sq=NULL,alpha.method="Cox",ind.mat=NULL,t.mat=NULL,GAMMA=NULL,z2=NULL){
  v.vec <- c(0,v)

  if(!is.null(GAMMA)){
    print("post.fn: z2-by-class interaction included")
    GAMMA.mat <- as.matrix(cbind(0,GAMMA))
  }

  J <- nclass
  P <- ncol(des)
  fullkappa.mat <- rbind(rep(0,P),kappa.mat)

  if(is.null(X)){
    print("post.fn: No categorical manifest variables")
  }
  if(is.null(X2) | is.null(mu) | is.null(sigma.sq)){
    print("post.fn: No continuous manifest variables")
  }
  if(!is.null(X)){
    n <- nrow(X)   
  } else {
    n <- nrow(X2)
  }

  log.w.num <- w <- matrix(NA,n,J)
  PEM.post.contribution.j <- PEM.post.contribution.i.j <- vector("list",J)

    for(j in 1:J){
      if(!is.null(X)){
        log.w.num.j.X <- log(apply(X,1,function(row) prod(p[j,]^row)))
      } else { 
        log.w.num.j.X <- log(1)
      }
      if(!is.null(X2)){ 
        log.w.num.j.X2 <- mvtnorm::dmvnorm(as.matrix(X2),mu[j,],sigma=as.matrix(sigma.sq),log=TRUE)
      } else {
        log.w.num.j.X2 <- log(1)
      }
      if(alpha.method=="Cox"){
        A <- cumsum(alpha)
        #log.w.num[,j] <- des %*% fullkappa.mat[j,] + log.w.num.j.X + log.w.num.j.X2 + v.vec[j]*delta -A*exp(as.matrix(z) %*% beta + v.vec[j] + if(!is.null(GAMMA)){as.matrix(z2) %*% GAMMA.mat[,j]} else {0})
        log.w.num[,j] <- des %*% fullkappa.mat[j,] + log.w.num.j.X + log.w.num.j.X2 + (v.vec[j] + if(!is.null(GAMMA)){as.matrix(z2) %*% GAMMA.mat[,j]} else {0})*delta -A*exp(as.matrix(z) %*% beta + v.vec[j] + if(!is.null(GAMMA)){as.matrix(z2) %*% GAMMA.mat[,j]} else {0})
      }
      if(alpha.method=="PEM"){
        S <- ncol(t.mat)
        PEM.post.contribution.j[[j]] <- matrix(NA,n,S)

        for(s in 1:S){ ##NB: first part of likelihood cancels out with denominator
          PEM.post.contribution.j[[j]][,s] <- exp(z %*% as.matrix(beta) + v.vec[j] + if(!is.null(GAMMA)){as.matrix(z2) %*% GAMMA.mat[,j]} else {0})^(ind.mat[,s]*delta) *
            exp( -ind.mat[,s]*( alpha[s]*t.mat[,s] + cbind(0,t(apply(t.mat,1,function(x) cumsum(x*alpha))))[,s] )*exp(z %*% as.matrix(beta) + v.vec[j] + if(!is.null(GAMMA)){as.matrix(z2) %*% GAMMA.mat[,j]} else {0} ) ) 
        } ##NB: the double-exponential is problematic here
        PEM.post.contribution.i.j[[j]] <- apply(PEM.post.contribution.j[[j]],1,function(x) prod(x))
     
        log.w.num[,j] <- des %*% fullkappa.mat[j,] + log.w.num.j.X + log.w.num.j.X2 + log( PEM.post.contribution.i.j[[j]] )
      }
    }
  
  #w <- exp(log.w.num) / rowSums(exp(log.w.num))
  w <- cbind( exp(-log( 1+exp(log.w.num[,2]-log.w.num[,1]) )),exp(-log( 1+exp(log.w.num[,1]-log.w.num[,2]) )) ) ##log-sum-exp trick!
  if(is.nan(sum(w))){
    print("post.fn: NaNs produced")
  }
  return(w)
}
#dim(p)
#dim(p.true)
#alpha.start <- rep(0.01,n.part)
#sigma.sq.true.matrix <- matrix(0,length(sigma.sq.true),length(sigma.sq.true))
#diag(sigma.sq.true.matrix) <- sigma.sq.true
#post <- post.fn12(kappa.mat=kappa.mat.true,alpha=alpha.start,p=p,beta=beta.true,v=v.true,des=cbind(1,x),X=X4,z=z,delta=delta,nclass=2,alpha.method="Cox",ind.mat=ind.mat,t.mat=t.mat,X2=X2,mu=mu.true,sigma.sq=sigma.sq.true.matrix )
#post <- post.fn12(kappa.mat=kappa.mat.true,alpha=alpha.start,p=p,beta=beta.true,v=v.true,des=cbind(1,x),X=X4,z=z,delta=delta,nclass=2,alpha.method="PEM",ind.mat=ind.mat,t.mat=t.mat,X2=X2,mu=mu.true,sigma.sq=sigma.sq.true.matrix )
#table(y2,round(post[,2]))
#post1 <- post.fn12(kappa.mat=kappa.mat.true,alpha=alpha.start,p=p,beta=beta.true,v=v.true,des=cbind(1,x),X=X4,z=z,delta=delta,nclass=2,alpha.method="Cox",ind.mat=NULL,t.mat=NULL,X2=X2,mu=mu.true,sigma.sq=sigma.sq.true.matrix,U2=0,v2=NULL)
#post2 <- post.fn11(kappa.mat=kappa.mat.true,alpha=alpha.start,p=p,beta=beta.true,v=v.true,des=cbind(1,x),X=X4,z=z,delta=delta,nclass=2,alpha.method="Cox",ind.mat=NULL,t.mat=NULL,X2=X2,mu=mu.true,sigma.sq=sigma.sq.true.matrix)
#max(post1-post2)
#post <- post.fn12(kappa.mat=kappa.mat.true,alpha=alpha.start,p=p,beta=beta.true,v=v.true,des=cbind(1,x),X=X4,z=z,delta=delta,nclass=2,alpha.method="Cox",ind.mat=NULL,t.mat=NULL,X2=X2,mu=mu.true,sigma.sq=sigma.sq.true.matrix)
#post <- post.fn12(kappa.mat=kappa.mat.true,alpha=alpha.start,p=p,beta=beta.true,v=v.true,des=cbind(1,x),X=X4,z=z,delta=delta,nclass=2,alpha.method="Cox",ind.mat=NULL,t.mat=NULL,X2=X2,mu=mu.true,sigma.sq=sigma.sq.true.matrix,GAMMA=NULL,z2=NULL)
#post <- post.fn12(kappa.mat=kappa.mat.true,alpha=alpha.PEM.start,p=p,beta=beta.true,v=v.true,des=cbind(1,x),X=X4,z=z,delta=delta,nclass=2,alpha.method="PEM",ind.mat=ind.mat,t.mat=t.mat,X2=X2,mu=mu.true,sigma.sq=sigma.sq.true.matrix,GAMMA=NULL,z2=NULL)
#post <- post.fn12(kappa.mat=kappa.mat.true,alpha=alpha.start,p=p,beta=beta.true,v=v.true,des=cbind(1,x),X=X4,z=z,delta=delta,nclass=2,alpha.method="Cox",ind.mat=NULL,t.mat=NULL,X2=X2,mu=mu.true,sigma.sq=sigma.sq.true.matrix,GAMMA=GAMMA.true,z2=z)
#post <- post.fn12(kappa.mat=kappa.mat.true,alpha=alpha.PEM.start,p=p,beta=beta.true,v=v.true,des=cbind(1,x),X=X4,z=z,delta=delta,nclass=2,alpha.method="PEM",ind.mat=ind.mat,t.mat=t.mat,X2=X2,mu=mu.true,sigma.sq=sigma.sq.true.matrix,GAMMA=GAMMA.true,z2=z)

#test.post <- post.fn12(kappa.mat=kappa.mat.true,alpha=alpha.start[1:6],p=p,beta=beta.true,v=v.true,des=cbind(1,x[1:6]),X=X4[1:6,],z=z[1:6],delta=delta[1:6],nclass=2,alpha.method="Cox",ind.mat=NULL,t.mat=NULL,X2=NULL,mu=NULL,sigma.sq=NULL,GAMMA=GAMMA.true,z2=z[1:6])
#test.post


###########################################################################
loglik.fn12 <- function(w,kappa.mat,alpha,p,beta,v,des,X,z,delta,nclass,X2=NULL,mu=NULL,sigma.sq=NULL,alpha.method="Cox",ind.mat=NULL,t.mat=NULL,GAMMA=NULL,z2=NULL){
  v.vec <- c(0,v)

  if(!is.null(GAMMA)){
    print("loglik.fn: z2-by-class interaction included")
    GAMMA.mat <- as.matrix(cbind(0,GAMMA))
  }

  J <- nclass
  P <- ncol(des)
  fullkappa.mat <- rbind(rep(0,P),kappa.mat)

  if(is.null(X)){
    print("loglik.fn: No categorical manifest variables")
  }
  if(is.null(X2) | is.null(mu) | is.null(sigma.sq)){
    print("loglik.fn: No continuous manifest variables")
  }

  if(!is.null(X)){
    n <- nrow(X)
    if(min(p,(1-p))==0){
      print("LOG-LIKELIHOOD NOT CALCULABLE SINCE phat = 0/1") ##since we can't take log of zero
    }
  } else {
    n <- nrow(X2)
  }

  loglik.i2 <- matrix(NA,n,J)
 
  for(j in 1:J){
    if(!is.null(X)){
      loglik.i2.j.X  <- apply(X,1,function(row) sum(row*log(p[j,]) ))
    } else {
      loglik.i2.j.X <- 0
    } 
    if(!is.null(X2)){
      loglik.i2.j.X2 <- mvtnorm::dmvnorm(as.matrix(X2),mu[j,],sigma=as.matrix(sigma.sq),log=TRUE)
    } else {
      loglik.i2.j.X2 <- 0
    }

    if(alpha.method=="Cox"){
      print("loglik.fn: alpha.method is Cox")
      A <- cumsum(alpha)

      loglik.i2[,j] = des %*% fullkappa.mat[j,] - log( rowSums(apply(fullkappa.mat,1,function(kappa.j.vec) exp(des %*% kappa.j.vec))) ) + 
        loglik.i2.j.X + loglik.i2.j.X2 + ifelse(delta==1,delta*(log(alpha) + as.matrix(z) %*% beta + v.vec[j] + if(!is.null(GAMMA)){as.matrix(z2) %*% GAMMA.mat[,j]} else {0}) ,0) -
        A*exp(as.matrix(z) %*% beta + v.vec[j] + if(!is.null(GAMMA)){as.matrix(z2) %*% GAMMA.mat[,j]} else {0})  ##ifelse removes the possibility of log(0)=-Inf
    } else if(alpha.method=="PEM"){
      print("loglik.fn: alpha.method is PEM")
	  
	S <- ncol(t.mat)
      PEM.loglik.contribution.j <- PEM.loglik.contribution.i.j <- vector("list",J)
      PEM.loglik.contribution.j[[j]] <- matrix(NA,n,S)
      for(s in 1:S){
        PEM.loglik.contribution.j[[j]][,s] <- ind.mat[,s]*delta * log(alpha[s]*exp(z %*% as.matrix(beta)+ v.vec[j] + if(!is.null(GAMMA)){as.matrix(z2) %*% GAMMA.mat[,j]} else {0})) -
          ind.mat[,s]*( alpha[s]*t.mat[,s] + cbind(0,t(apply(t.mat,1,function(x) cumsum(x*alpha))))[,s] ) * exp(z %*% as.matrix(beta) + v.vec[j] + if(!is.null(GAMMA)){as.matrix(z2) %*% GAMMA.mat[,j]} else {0})
      }
      PEM.loglik.contribution.i.j[[j]] <- rowSums(PEM.loglik.contribution.j[[j]])	  
    
      loglik.i2[,j] <- des %*% fullkappa.mat[j,] - log( rowSums(apply(fullkappa.mat,1,function(kappa.j.vec) exp(des %*% kappa.j.vec)))) +
        loglik.i2.j.X + loglik.i2.j.X2 + PEM.loglik.contribution.i.j[[j]]
      }
  }

  L.exp <- sum(w*loglik.i2) 
  L.obs <- sum(log((rowSums(exp(loglik.i2))))) 

  return.list <- list(loglik.i2,L.exp,L.obs)
  names(return.list) <- c("loglik.i2","L.exp","L.obs")
  return(return.list)
}
#alpha.start <- rep(0.1,n)
#dim(p.true)
#dim(p)
#loglik.fn12(w=post,kappa.mat=kappa.mat.true,alpha=alpha.start,p=p,beta=beta.true,v=v.true,des=cbind(1,x),X=X4,z=z,delta=delta,nclass=2,X2=X2,mu=mu.true,sigma.sq=sigma.sq.true.matrix,alpha.method="Cox")
#loglik.fn12(w=post,kappa.mat=kappa.mat.true,alpha=alpha.start,p=p,beta=beta.true,v=v.true,des=cbind(1,x),X=X4,z=z,delta=delta,nclass=2,X2=X2,mu=mu.true,sigma.sq=sigma.sq.true.matrix,alpha.method="PEM",ind.mat=ind.mat,t.mat=t.mat)
#loglik.fn12(w=post,kappa.mat=kappa.mat.true,alpha=alpha.start,p=p,beta=beta.true,v=v.true,des=cbind(1,x),X=X4,z=z,delta=delta,nclass=2,X2=X2,mu=mu.true,sigma.sq=sigma.sq.true.matrix,alpha.method="Cox",U2=0,v2=NULL)
#loglik.fn11(w=post,kappa.mat=kappa.mat.true,alpha=alpha.start,p=p,beta=beta.true,v=v.true,des=cbind(1,x),X=X4,z=z,delta=delta,nclass=2,X2=X2,mu=mu.true,sigma.sq=sigma.sq.true.matrix,alpha.method="Cox")
#loglik.fn12(w=post,kappa.mat=kappa.mat.true,alpha=alpha.start,p=p,beta=beta.true,v=v.true,des=cbind(1,x),X=X4,z=z,delta=delta,nclass=2,X2=X2,mu=mu.true,sigma.sq=sigma.sq.true.matrix,alpha.method="Cox",U2=U,v2=v2.test)
#loglik.fn12(w=post,kappa.mat=kappa.mat.true,alpha=alpha.start,p=p,beta=beta.true,v=v.true,des=cbind(1,x),X=X4,z=z,delta=delta,nclass=2,X2=X2,mu=mu.true,sigma.sq=sigma.sq.true.matrix,alpha.method="Cox",GAMMA=GAMMA.true,z2=z)
#loglik.fn12(w=post,kappa.mat=kappa.mat.true,alpha=alpha.PEM.start,p=p,beta=beta.true,v=v.true,des=cbind(1,x),X=X4,z=z,delta=delta,nclass=2,X2=X2,mu=mu.true,sigma.sq=sigma.sq.true.matrix,alpha.method="PEM",t.mat=t.mat,ind.mat=ind.mat,GAMMA=GAMMA.true,z2=z)



###########################################################################
mu.sigma.sq.fn12 <- function(post,X2,nclass,cov.type="homo",custom.zeros=NULL){

  ##cov.type="homo" - unstructured and homogenous across classes
  ##cov.type="localind" - homogenous across classes and diagonial (zero covariance within class)
  ##cov.type="custom" - need to specify which elements of the lower.tri of the cov matrix should equal zero

  X2 <- as.matrix(X2)
  J <- nclass
  n <- nrow(X2)

  mu <- matrix(NA,nclass,ncol(X2))
  L <- ncol(X2)

  for(j in 1:J){
    mu[j,] <- colSums(post[,j]*X2) / sum(post[,j])
  }

  ##Update cov (see Mclachlan and Peel - pages 82-82, equations 3.4 & 3.10)
  cov.list <- sweep.list <- vector("list",J)
  for(j in 1:J){
    sweep.list[[j]] <- sweep(as.matrix(X2),2,mu[j,],"-")
    cov.list[[j]] <- Reduce("+",  lapply( 1:n,function(x) post[x,j]*(as.matrix(sweep.list[[j]][x,]) %*% t(as.matrix(sweep.list[[j]][x,]))) )  ) / sum(post[,j])
  }
  cov <- as.matrix(( sum(post[,1])*cov.list[[1]] + sum(post[,2])*cov.list[[2]] ) / n)

  if(cov.type=="localind"){
    zeros <- matrix(0,ncol(as.matrix(X2)),ncol(as.matrix(X2)))
    diag(zeros) <- diag(cov)
    cov <- zeros
  }
  if(cov.type=="custom"){
    cov[lower.tri(cov)][custom.zeros] <- 0
    cov[upper.tri(cov)] <- t(cov)[upper.tri(cov)]
  }
  sigma.sq <- cov
  
  return.list <- list(mu,sigma.sq)
  names(return.list) <- c("mu","sigma.sq")
  return(return.list)
}
#mu.sigma.sq.fn12(post=post,X2=X2,nclass=2)


###########################################################################
phat.fn12 <- function(w,X,nclass){
  J <- nclass
  M <- ncol(X)
  phat2 <- matrix(NA,J,M)

  for(j in 1:J){
    phat2[j,] <- apply(X,2,function(col) sum(col*w[,j]) / sum(w[,j]))
  }
  return(phat2)
}
#round(phat.fn12(w=post,X=X4,nclass=2),2)
#p



###########################################################################
Skappa.fn12 <- function(w,kappa.mat,des,nclass){
  J <- nclass
  P <- ncol(des)
  fullkappa.mat <-  rbind(rep(0,P),kappa.mat)
  Sk <- matrix(NA,(J-1),P)
  
  for(j in 2:J){
    for(p in 1:P){
      Sk[(j-1),p] <- colSums(des[,p]*(w[,j] - exp(des %*% fullkappa.mat[j,]) / (rowSums(apply(fullkappa.mat,1,function(kappa.j.vec) exp(des %*% kappa.j.vec))))))
    }
  }

  return(Sk)
}


###########################################################################
kappaNR.fn12 <- function(w,kappa.mat,des,nclass){
  J <- nclass
  P <- ncol(des)
  fullkappa.mat <-  rbind(rep(0,P),kappa.mat)

  #FIRST DERIVATIVES
  Skappa <- Skappa.fn12(w,kappa.mat,des,nclass)

  #SECOND DERIVATIVES
  Dkappa.vec <- matrix(NA,(P*P*(J-1)*(J-1)),1)
 
  counter = 0
  denom <- rowSums(apply(fullkappa.mat,1,function(row) exp(des %*% row)))

  for(j1 in 1:(J-1)){
    for(j2 in 1:(J-1)){
      for(q in 1:P){
        for(p in 1:P){
          counter <- counter + 1
          
          if(j1 == j2){
            ind = 1
          } else {
            ind = 0
          }

          Dkappa.vec[counter] <- colSums( des[,p]*des[,q] * (ind - exp(des %*% kappa.mat[j1,]) / denom) *
            exp(des %*% kappa.mat[j2,]) / denom )
        }
      }
    }
  }
  ######################
  print("Dkappa.vec is:")
  print(Dkappa.vec)
  ######################
  Dkappa2 <- matrix(Dkappa.vec,sqrt(length(Dkappa.vec)),sqrt(length(Dkappa.vec)))
  ######################
  print("Dkappa2 is:")
  print(Dkappa2)
  ######################

  ##Newton-Raphson with error trap for singular matrix
  if(is.nan(sum(Dkappa2)) || matrixcalc::is.singular.matrix(Dkappa2)){
    Dkappa2.break = 1
    newkappa = NA
    print("Dkappa2 is singular (kappaNR.fn12)")
  } else {
    Dkappa2.break = 0
    newkappa <- t(kappa.mat) + solve(Dkappa2) %*% t(Skappa)
  }

  return.list <- list(newkappa,Dkappa2,Skappa,Dkappa2.break)
  names(return.list) <- c("newkappa","Dkappa2","Skappa","Dkappa2.break")
  return(return.list)
}


###########################################################################
Sbeta.fn12 <- function(w,alpha,beta,v,delta,z,alpha.method="Cox",ind.mat=NULL,t.mat=NULL,GAMMA=NULL,z2=NULL){
  v.vec <- c(0,v)
  Q <- length(beta)

  if(!is.null(GAMMA)){
    GAMMA.mat <- as.matrix(cbind(0,GAMMA))
    Q3 <- nrow(GAMMA.mat)
    SGAMMA <- matrix(NA,Q3,1)
  } else {
    SGAMMA <- NULL
  }

  if(alpha.method=="Cox"){
    Sbeta <- matrix(NA,Q,1)
    A <- cumsum(alpha)
    for(q in 1:Q){
      #Sbeta[q] <- sum( as.matrix(z)[,q]*(delta - A*exp(as.matrix(z) %*% beta)*apply(w,1,function(row) sum(row*exp(v.vec + if(!is.null(GAMMA)){as.matrix(z2) %*% GAMMA.mat} else {0})))) )
      Sbeta[q] <- sum( as.matrix(z)[,q]*(delta - A*exp(as.matrix(z) %*% beta)*
          ((w[,1] * exp(v.vec[1] + if(!is.null(GAMMA)){as.matrix(z2) %*% GAMMA.mat[,1]} else {0}) +  
          w[,2] * exp(v.vec[2] + if(!is.null(GAMMA)){as.matrix(z2) %*% GAMMA.mat[,2]} else {0})
	    ))) )
    }
    if(!is.null(GAMMA)){
      for(Q3 in 1:Q3){
        SGAMMA[Q3] <- sum( w[,2]*as.matrix(z2)[,Q3]*(delta - A*exp(as.matrix(z) %*% beta + v.vec[2] + as.matrix(z2) %*% GAMMA.mat[,2])) )
      }
    }
  } else if(alpha.method=="PEM"){
    
    S <- ncol(t.mat)
    Sbeta.q.s <- matrix(NA,Q,S)

    for(q in 1:Q){
      for(s in 1:S){
        #Sbeta.q.s[q,s] <- sum( as.matrix(z)[,q]*(ind.mat[,s]*delta - 
        #  (ind.mat[,s]*(alpha[s]*t.mat[,s] + cbind(0,t(apply(t.mat,1,function(x) cumsum(x*alpha))))[,s])*
        #  exp(as.matrix(z) %*% beta)*apply(w,1,function(row) sum(row*exp(v.vec + if(!is.null(GAMMA)){as.matrix(z2) %*% GAMMA.mat[,2]} else {0}))))) )
        Sbeta.q.s[q,s] <- sum( as.matrix(z)[,q]*(ind.mat[,s]*delta - 
          (ind.mat[,s]*(alpha[s]*t.mat[,s] + cbind(0,t(apply(t.mat,1,function(x) cumsum(x*alpha))))[,s])*
          exp(as.matrix(z) %*% beta)*(w[,1] * exp(v.vec[1] + if(!is.null(GAMMA)){as.matrix(z2) %*% GAMMA.mat[,1]} else {0}) +  
          w[,2] * exp(v.vec[2] + if(!is.null(GAMMA)){as.matrix(z2) %*% GAMMA.mat[,2]} else {0})
	    ))) )
      }
      Sbeta <- as.matrix(rowSums(Sbeta.q.s))
    }
    if(!is.null(GAMMA)){
      SGAMMA.Q3.s <- matrix(NA,Q3,S)
      for(Q3 in 1:Q3){
        for(s in 1:S){
          SGAMMA.Q3.s[Q3,s] <- sum( w[,2]*as.matrix(z2)[,Q3]*(ind.mat[,s]*delta - 
            (ind.mat[,s]*(alpha[s]*t.mat[,s] + cbind(0,t(apply(t.mat,1,function(x) cumsum(x*alpha))))[,s])*
            exp(as.matrix(z) %*% beta + v.vec[2] + if(!is.null(GAMMA)){as.matrix(z2) %*% GAMMA.mat[,2]} else {0}))) )
        }
      }
      SGAMMA <- as.matrix(rowSums(SGAMMA.Q3.s))
    }
  }

  names(Sbeta) <- rep("Sbeta",Q)

  if(!is.null(GAMMA)){
    names(SGAMMA) <- rep("SGAMMA",Q3) 
  }

  return(c(Sbeta,SGAMMA))
}
#Sbeta.fn12(w=post,alpha=rep(0.01,n),beta=beta.true,v=v.true,delta=delta,z=z,alpha.method="Cox")
#Sbeta.fn12(w=post,alpha=alpha.start,beta=beta.true,v=v.true,delta=delta,z=z,alpha.method="PEM",ind.mat=ind.mat,t.mat=t.mat)
#Sbeta.fn12(w=post,alpha=rep(0.01,n),beta=beta.true,v=v.true,delta=delta,z=z,alpha.method="Cox",U2=0,v2=NULL)
#Sbeta.fn11(w=post,alpha=rep(0.01,n),beta=beta.true,v=v.true,delta=delta,z=z,alpha.method="Cox")
#Sbeta.fn12(w=post,alpha=rep(0.01,n),beta=beta.true,v=v.true,delta=delta,z=z,alpha.method="Cox",GAMMA=GAMMA.true,z2=z)
#Sbeta.fn12(w=post,alpha=alpha.PEM.start,beta=beta.true,v=v.true,delta=delta,z=z,alpha.method="PEM",ind.mat=ind.mat,t.mat=t.mat,GAMMA=GAMMA.true,z2=z)
#test.Sbeta <- Sbeta.fn12(w=test.post,alpha=alpha.start[1:6],beta=beta.true,v=v.true,delta=delta[1:6],z=z[1:6],alpha.method="Cox",GAMMA=GAMMA.true,z2=z[1:6])
#test.Sbeta
#test.Sbeta <- Sbeta.fn12(w=test.post,alpha=alpha.start[1:3],beta=beta.true,v=v.true,delta=delta[1:6],z=z[1:6],alpha.method="PEM",GAMMA=GAMMA.true,z2=z[1:6],ind.mat=ind.mat[1:6,],t.mat=t.mat[1:6,])
#test.Sbeta


###########################################################################
Sv.fn12 <- function(w,alpha,beta,v,delta,z,nclass,alpha.method="Cox",ind.mat=NULL,t.mat=NULL,GAMMA=NULL,z2=NULL){
  J <- nclass

  if(!is.null(GAMMA)){
    GAMMA.mat <- as.matrix(cbind(0,GAMMA))
  }

  if(alpha.method=="Cox"){
    Sv <- matrix(NA,length(v),1) 
    A <- cumsum(alpha)
    for(j in 2:J){
      Sv[j-1] <- sum( w[,j]*(delta - A*exp(as.matrix(z) %*% beta + v[j-1] + if(!is.null(GAMMA)){as.matrix(z2) %*% GAMMA.mat[,j]} else {0})) )
    }
  } else if(alpha.method=="PEM"){
    S <- ncol(t.mat)
    Sv.s <- matrix(NA,length(v),S)

    for(j in 2:J){
      for(s in 1:S){
        Sv.s[(j-1),s] <- sum( w[,j]*ind.mat[,s]*(delta - (alpha[s]*t.mat[,s] + cbind(0,t(apply(t.mat,1,function(x) cumsum(x*alpha))))[,s])*exp(as.matrix(z) %*% beta + v[j-1] + if(!is.null(GAMMA)){as.matrix(z2) %*% GAMMA.mat[,2]} else {0})) )  
      }
    }
    Sv <- sum(Sv.s)
  }

  names(Sv) <- c("Sv")
  return(Sv)
}
#Sv.fn12(w=post,alpha=alpha.start,beta=beta.true,v=v.true,delta=delta,z=z,nclass=2,alpha.method="Cox")
#Sv.fn12(w=post,alpha=alpha.start,beta=beta.true,v=v.true,delta=delta,z=z,nclass=2,alpha.method="PEM",ind.mat=ind.mat,t.mat=t.mat)
#Sv.fn12(w=post,alpha=alpha.start,beta=beta.true,v=v.true,delta=delta,z=z,nclass=2,alpha.method="Cox",U2=0,v2=NULL)
#Sv.fn11(w=post,alpha=alpha.start,beta=beta.true,v=v.true,delta=delta,z=z,nclass=2,alpha.method="Cox")
#Sv.fn12(w=post,alpha=alpha.start,beta=beta.true,v=v.true,delta=delta,z=z,nclass=2,alpha.method="Cox",U2=U,v2=v2.test)
#Sv.fn12(w=post,alpha=alpha.start,beta=beta.true,v=v.true,delta=delta,z=z,nclass=2,alpha.method="Cox",GAMMA=GAMMA.true,z2=z)
#Sv.fn12(w=post,alpha=alpha.PEM.start,beta=beta.true,v=v.true,delta=delta,z=z,nclass=2,alpha.method="PEM",t.mat=t.mat,ind.mat=ind.mat,GAMMA=GAMMA.true,z2=z)

#test.Sv <- Sv.fn12(w=test.post,alpha=alpha.start[1:6],beta=beta.true,v=v.true,delta=delta[1:6],z=z[1:6],nclass=2,alpha.method="Cox",GAMMA=GAMMA.true,z2=z[1:6])
#test.Sv

###########################################################################
betaNR.fn12 <- function(w,alpha,beta,v,delta,z,nclass,alpha.method="Cox",ind.mat=NULL,t.mat=NULL,GAMMA=NULL,z2=NULL){
  J <- nclass
  v.vec <- c(0,v)
  Q <- Q1 <- Q2 <- length(beta)

  if(!is.null(GAMMA)){
    GAMMA.mat <- as.matrix(cbind(0,GAMMA))
    Q3 <- nrow(GAMMA.mat)
    SGAMMA <- matrix(NA,Q3,1)
  } else {
    SGAMMA <- NULL
  }

  if(alpha.method=="Cox"){
    ##COX
    print("betaNR.fn: alpha.method is Cox")
    A <- cumsum(alpha)

    ##FIRST DERIVATIVES
    Sbeta.SGAMMA <- Sbeta.fn12(w=w,alpha=alpha,beta=beta,v=v,delta=delta,z=z,GAMMA=GAMMA,z2=z2)
    Sv <- Sv.fn12(w=w,alpha=alpha,beta=beta,v=v,delta=delta,z=z,nclass=nclass,GAMMA=GAMMA,z2=z2)  

    Sbeta <- Sbeta.SGAMMA[1:Q1]

    if(is.null(GAMMA)){
      S. <- c(Sbeta,Sv)
    } else {
      SGAMMA <- Sbeta.SGAMMA[(Q1+1):(Q1+Q3)]
      S. <- c(Sbeta,Sv,SGAMMA)     
    }
    print(S.)


    ####NEGATIVE SECOND DERIVATIVES

    ##Second derivative for beta
    S2beta <- matrix(NA,Q,Q)

    for(q1 in 1:Q){
      for(q2 in 1:Q){
        #S2beta[q1,q2] <- sum( as.matrix(z)[,q1]*as.matrix(z)[,q2] * A*exp(as.matrix(z) %*% beta) * 
        #  apply(w,1,function(row) sum(row*exp(v.vec + if(!is.null(GAMMA)){as.matrix(z2) %*% GAMMA.mat} else {0}))) )
        S2beta[q1,q2] <- sum( as.matrix(z)[,q1]*as.matrix(z)[,q2] * A*exp(as.matrix(z) %*% beta) * 
          ((w[,1] * exp(v.vec[1] + if(!is.null(GAMMA)){as.matrix(z2) %*% GAMMA.mat[,1]} else {0}) +  
          w[,2] * exp(v.vec[2] + if(!is.null(GAMMA)){as.matrix(z2) %*% GAMMA.mat[,2]} else {0})
	    )) )
      }
    }

    ##Second derivative for v - NB: class 1 effect is zero.
    S2v <- matrix(NA,length(v),1)

    for(j in 2:J){
      S2v[j-1] <- sum( w[,j]*A*exp(as.matrix(z) %*% beta + v[j-1] + if(!is.null(GAMMA)){as.matrix(z2) %*% GAMMA.mat[,2]} else {0}) )
    }

    ##Differentiate Sbeta w.r.t. v (or vice-versa)
    S2beta.v <- matrix(NA,length(v),Q)

    for(j in 2:J){
      for(q in 1:Q){
        S2beta.v[(j-1),q] <- sum( A*as.matrix(z)[,q]*exp(as.matrix(z) %*% beta) * 
          w[,j]*exp(v[j-1] + if(!is.null(GAMMA)){as.matrix(z2) %*% GAMMA.mat[,2]} else {0}) )
      }
    }

    if(!is.null(GAMMA)){
      S2GAMMA <- matrix(NA,Q3,Q3)
      S2GAMMA.v <- matrix(NA,Q3,1)
      S2GAMMA.beta <- matrix(NA,Q3,Q1)

      ##Second derivative for GAMMA
      for(q3 in 1:Q3){
        for(q4 in 1:Q3){
          S2GAMMA[q3,q4] <- sum( w[,2]*as.matrix(z2)[,q3]*as.matrix(z2)[,q4]*(A*exp(as.matrix(z) %*% beta + v.vec[2] + as.matrix(z2) %*% GAMMA.mat[,2])) )
        }
      }
      
      ##Differentiate SGAMMA w.r.t. v
      for(q3 in 1:Q3){
        S2GAMMA.v[q3] <- sum( w[,2]*as.matrix(z2)[,q3]*(A*exp(as.matrix(z) %*% beta + v.vec[2] + as.matrix(z2) %*% GAMMA.mat[,2])) )
      }

      ##Differentiate SGAMMA w.r.t. beta
      for(q3 in 1:Q3){  ##length of GAMMA = Q3
        for(q1 in 1:Q1){  ##length of beta = Q1
          S2GAMMA.beta[q3,q1] <- sum( w[,2]*as.matrix(z2)[,q3]*as.matrix(z)[,q1]*(A*exp(as.matrix(z) %*% beta + v.vec[2] + as.matrix(z2) %*% GAMMA.mat[,2])) )
        }
      }
 
    }


  } else if(alpha.method=="PEM"){
    
    ##PEM
    print("betaNR.fn: alpha.method is PEM")

    ##FIRST DERIVATIVES
    Sbeta.SGAMMA <- Sbeta.fn12(w=w,alpha=alpha,beta=beta,v=v,delta=delta,z=z,alpha.method="PEM",ind.mat=ind.mat,t.mat=t.mat,GAMMA=GAMMA,z2=z2)
    Sv <- Sv.fn12(w=w,alpha=alpha,beta=beta,v=v,delta=delta,z=z,nclass=nclass,alpha.method="PEM",ind.mat=ind.mat,t.mat=t.mat,GAMMA=GAMMA,z2=z2)

    Sbeta <- Sbeta.SGAMMA[1:Q1]

    if(is.null(GAMMA)){
      S. <- c(Sbeta,Sv)
    } else {
      SGAMMA <- Sbeta.SGAMMA[(Q1+1):(Q1+Q3)]
      S. <- c(Sbeta,Sv,SGAMMA)     
    }
    print(S.)


    ##NEGATIVE SECOND DERIVATIVES
    S <- ncol(t.mat)

    ##Second derivative for beta
    S2beta.s <- vector("list",S)

    for(s in 1:S){
      S2beta.s[[s]] <-  matrix(NA,Q,Q)
      for(q1 in 1:Q){
        for(q2 in 1:Q){
          #S2beta.s[[s]][q1,q2] <- sum( as.matrix(z)[,q1]*as.matrix(z)[,q2] * 
          #  ind.mat[,s]*(alpha[s]*t.mat[,s] + cbind(0,t(apply(t.mat,1,function(x) cumsum(x*alpha))))[,s])*exp(as.matrix(z) %*% beta) * 
          #  apply(w,1,function(row) sum(row*exp(v.vec + if(!is.null(GAMMA)){as.matrix(z2) %*% GAMMA.mat} else {0}))) )
          S2beta.s[[s]][q1,q2] <- sum( as.matrix(z)[,q1]*as.matrix(z)[,q2] * 
            ind.mat[,s]*(alpha[s]*t.mat[,s] + cbind(0,t(apply(t.mat,1,function(x) cumsum(x*alpha))))[,s])*exp(as.matrix(z) %*% beta) * 
            (w[,1] * exp(v.vec[1] + if(!is.null(GAMMA)){as.matrix(z2) %*% GAMMA.mat[,1]} else {0}) +  
            w[,2] * exp(v.vec[2] + if(!is.null(GAMMA)){as.matrix(z2) %*% GAMMA.mat[,2]} else {0})
	      ) )
        }
      }
    }
    S2beta <- Reduce("+",S2beta.s)

    ##Second derivative for v - NB: class 1 effect is zero.
    S2v.s <- vector("list",S)
    
    for(s in 1:S){
      S2v.s[[s]] <- matrix(NA,length(v),1)

      for(j in 2:J){
       S2v.s[[s]][j-1] <- sum( w[,j]*ind.mat[,s]*(alpha[s]*t.mat[,s] + cbind(0,t(apply(t.mat,1,function(x) cumsum(x*alpha))))[,s]) *
         exp(as.matrix(z) %*% beta + v[j-1] + if(!is.null(GAMMA)){as.matrix(z2) %*% GAMMA.mat[,2]} else {0}) )
      }
    }
    S2v <- Reduce("+",S2v.s)

    ##Differentiate Sbeta w.r.t. v (or vice-versa)
    S2beta.v.s <- vector("list",S)
    
    for(s in 1:S){
      S2beta.v.s[[s]] <- matrix(NA,length(v),Q)

      for(j in 2:J){
        for(q in 1:Q){
          S2beta.v.s[[s]][(j-1),q] <- sum( as.matrix(z)[,q]*ind.mat[,s]*(alpha[s]*t.mat[,s] + cbind(0,t(apply(t.mat,1,function(x) cumsum(x*alpha))))[,s]) *
            exp(as.matrix(z) %*% beta)* w[,j]*exp(v[j-1] + if(!is.null(GAMMA)){as.matrix(z2) %*% GAMMA.mat[,2]} else {0}) )
        }
      }
    }
    S2beta.v <- Reduce("+",S2beta.v.s)


    if(!is.null(GAMMA)){
      S2GAMMA.s <- vector("list",S)
      S2GAMMA.v.s <- vector("list",S)
      S2GAMMA.beta <- vector("list",S)

      ##Second derivative for GAMMA
      for(s in 1:S){
        S2GAMMA.s[[s]] <- matrix(NA,Q3,Q3)
        for(q3 in 1:Q3){
          for(q4 in 1:Q3){
            S2GAMMA.s[[s]][q3,q4] <- sum( w[,2]*as.matrix(z2)[,q3]*as.matrix(z2)[,q4]*ind.mat[,s]*(alpha[s]*t.mat[,s] + cbind(0,t(apply(t.mat,1,function(x) cumsum(x*alpha))))[,s]) *
              exp(as.matrix(z) %*% beta + v.vec[2] + as.matrix(z2) %*% GAMMA.mat[,2]) )
          }
        }
      }
      S2GAMMA <- Reduce("+",S2GAMMA.s)
   
      ##Differentiate SGAMMA w.r.t. v
      S2GAMMA.v.s <- vector("list",S)
    
      for(s in 1:S){
        S2GAMMA.v.s[[s]] <- matrix(NA,Q3,1)

        for(q3 in 1:Q3){
          S2GAMMA.v.s[[s]][q3] <- sum( w[,2]*as.matrix(z2)[,q3]*ind.mat[,s]*(alpha[s]*t.mat[,s] + cbind(0,t(apply(t.mat,1,function(x) cumsum(x*alpha))))[,s]) *
            exp(as.matrix(z) %*% beta + v.vec[2] + as.matrix(z2) %*% GAMMA.mat[,2]) )
        }
      }
      S2GAMMA.v <- Reduce("+",S2GAMMA.v.s)


      ##Differentiate SGAMMA w.r.t. beta
      S2GAMMA.beta.s <- vector("list",S)

      for(s in 1:S){
        S2GAMMA.beta.s[[s]] <- matrix(NA,Q3,Q1)
        for(q3 in 1:Q3){  ##length of GAMMA = Q3
          for(q1 in 1:Q1){  ##length of beta = Q1
            S2GAMMA.beta.s[[s]][q3,q1] <- sum( w[,2]*as.matrix(z2)[,q3]*as.matrix(z)[,q1]*ind.mat[,s]*(alpha[s]*t.mat[,s] + cbind(0,t(apply(t.mat,1,function(x) cumsum(x*alpha))))[,s]) *
              exp(as.matrix(z) %*% beta + v.vec[2] + as.matrix(z2) %*% GAMMA.mat[,2]) )
          }
        }
      }
      S2GAMMA.beta <- Reduce("+",S2GAMMA.beta.s)


    } ##end !is.null(GAMMA)

  } ##end PEM loop

  ##Newton-Raphson
  if(is.null(GAMMA)){
    D <- cbind(rbind(S2beta,S2beta.v),rbind(t(S2beta.v),S2v))
    newbeta.v <- as.vector(c(beta,v)) + solve(D) %*% S.
  } else {
    D <- cbind(rbind(S2beta,S2beta.v,S2GAMMA.beta),rbind(t(S2beta.v),S2v,S2GAMMA.v),rbind(t(S2GAMMA.beta),t(S2GAMMA.v),S2GAMMA))
    newbeta.v <- as.vector(c(beta,v,GAMMA)) + solve(D) %*% S.    
  }

  return.list <- list(newbeta.v,D,S.)
  names(return.list) <- c("newbeta.v","D","S.")
  return(return.list)
}
#a$betaNR
#betaNR.fn12(w=a$post,alpha=rep(mean(a$alpha),4),beta=a$betaNR$newbeta.v[1],v=a$betaNR$newbeta.v[2],delta=delta,z=z,nclass=2,alpha.method="PEM",ind.mat=ind.mat,t.mat=t.mat)
#betaNR.fn12(w=post,alpha=alpha.start,beta=beta.true,v=v.true,delta=delta,z=z,nclass=2,alpha.method="Cox",ind.mat=NULL,t.mat=NULL,U2=0,v2=NULL)
#betaNR.fn11(w=post,alpha=alpha.start,beta=beta.true,v=v.true,delta=delta,z=z,nclass=2,alpha.method="Cox",ind.mat=NULL,t.mat=NULL)
#betaNR.fn12(w=post,alpha=alpha.start,beta=beta.true,v=v.true,delta=delta,z=z,nclass=2,alpha.method="Cox",ind.mat=NULL,t.mat=NULL,U2=U,v2=v2.test)
#betaNR.fn12(w=post,alpha=alpha.start,beta=beta.true,v=v.true,delta=delta,z=z,nclass=2,alpha.method="Cox",ind.mat=NULL,t.mat=NULL,U2=log(U)-mean(log(U)),v2=v2.test)
#betaNR.fn12(w=post,alpha=alpha.start,beta=beta.true,v=v.true,delta=delta,z=z,nclass=2,alpha.method="Cox",ind.mat=NULL,t.mat=NULL,GAMMA=NULL,z2=NULL)
#betaNR.fn12(w=post,alpha=alpha.start,beta=beta.true,v=v.true,delta=delta,z=z,nclass=2,alpha.method="Cox",ind.mat=NULL,t.mat=NULL,GAMMA=GAMMA.true,z2=z)
#betaNR.fn12(w=post,alpha=alpha.PEM.start,beta=beta.true,v=v.true,delta=delta,z=z,nclass=2,alpha.method="PEM",ind.mat=ind.mat,t.mat=t.mat,GAMMA=GAMMA.true,z2=z)

#test.betaNR <- betaNR.fn12(w=test.post,alpha=alpha.start[1:6],beta=beta.true,v=v.true,delta=delta[1:6],z=z[1:6],nclass=2,alpha.method="Cox",ind.mat=NULL,t.mat=NULL,GAMMA=GAMMA.true,z2=z[1:6])
#test.betaNR


###########################################################################
alpha.fn12 <- function(w,beta,v,z,U,delta,alpha.method="Cox",ind.mat=NULL,t.mat=NULL,GAMMA=NULL,z2=NULL){
  v.vec <- c(0,v)

  if(!is.null(GAMMA)){
    GAMMA.mat <- as.matrix(cbind(0,GAMMA))
  }

  times <- factor(U)
  n <- length(U)

  if(alpha.method=="Cox"){
    alpha.denom.i <- matrix(NA,n,1)
    alpha0 <- matrix(NA,length(times),1)
    rownames(alpha0) <- times

    #for(i in 1:n){
    #  alpha.denom.i[i] <- exp(as.matrix(z)[i,] %*% beta)*sum(w[i,]*exp(v.vec + if(!is.null(GAMMA)){as.matrix(z2)[i,] %*% GAMMA.mat} else {0}))
    #  #alpha.denom.i[i] <- exp(as.matrix(z)[i,] %*% beta)*(w[i,1] + w[i,2]*exp(v.vec[2] + if(!is.null(GAMMA)){as.matrix(z2)[i,] %*% GAMMA.mat[,2]} else {0}))
    #}

    alpha.denom.i <- exp(as.matrix(z) %*% beta)*(w[,1] + w[,2]*exp(v.vec[2] + if(!is.null(GAMMA)){as.matrix(z2) %*% GAMMA.mat[,2]} else {0}))

    for(t in 1:n){
      alpha0[t] <- delta[t] / sum(alpha.denom.i[t:n])
    }

  } else if(alpha.method=="PEM"){
    S <- ncol(t.mat)
    alpha0.s <- matrix(NA,1,S)
      for(s in 1:S){
        #alpha0.s[s] <- sum(  ind.mat[,s]*delta / sum( t.mat[,s] * exp(as.matrix(z) %*% beta) * rowSums(t(apply(w,1,function(x) x*exp(v.vec)))) )  )
        alpha0.s[s] <- sum(  ind.mat[,s]*delta / sum( t.mat[,s] * exp(as.matrix(z) %*% beta) * (w[,1] + w[,2]*exp(v.vec[2] + if(!is.null(GAMMA)){as.matrix(z2) %*% GAMMA.mat[,2]} else {0})) )  )
      }
    alpha0 <- colSums(alpha0.s)
    }

  return(alpha0)
}
#alpha.fn12(w=post,beta=beta.true,v=v.true,z=z,U=U,delta=delta,alpha.method="Cox")
#alpha.fn12(w=post,beta=beta.true,v=v.true,z=z,U=U,delta=delta,alpha.method="PEM",ind.mat=ind.mat,t.mat=t.mat)
#alpha.fn12(w=post,beta=beta.true,v=v.true,z=z,U=U,delta=delta,alpha.method="Cox",U2=0,v2=NULL)
#alpha.fn11(w=post,beta=beta.true,v=v.true,z=z,U=U,delta=delta,alpha.method="Cox")
#alpha.fn12(w=post,beta=beta.true,v=v.true,z=z,U=U,delta=delta,alpha.method="Cox")
#alpha.fn12(w=post,beta=beta.true,v=v.true,z=z,U=U,delta=delta,alpha.method="Cox",GAMMA=GAMMA.true,z2=z)
#alpha.fn12(w=post,beta=beta.true,v=v.true,z=z,U=U,delta=delta,alpha.method="PEM",t.mat=t.mat,ind.mat=ind.mat,GAMMA=GAMMA.true,z2=z)



#alpha.fn12(w=test.post,beta=beta.true,v=v.true,z=z[1:6],U=U[1:6],delta=delta[1:6],alpha.method="Cox",GAMMA=GAMMA.true,z2=z[1:6])


###########################################################################
#lcsm.fn12 <- function(X,x=NULL,z,U,delta,nclass=2,start.list=NULL,hist=TRUE,tol=1e-5,nreps=1,maxiter=1000,printiter=TRUE,boundary.crit=0.01,maxlik.method="L.obs",lambda1=1,lambda2=1,X2=X2,alpha.method="Cox",ind.mat=NULL,t.mat=NULL,cov.type="homo",X.var.vec=NULL,custom.zeros=NULL,z2=NULL){
LCSM <- function(y=NULL,x=NULL,z,t,delta,nclass=2,start.list=NULL,hist=TRUE,tol=1e-5,nreps=1,maxiter=1000,printiter=TRUE,boundary.crit=0.01,maxlik.method="L.obs",lambda1=1,lambda2=1,w=NULL,alpha.method="Cox",ind.mat=NULL,t.mat=NULL,cov.type="homo",y.var.vec=NULL,custom.zeros=NULL,z2=NULL){



##Initial step to switch argument names to match those listed in the thesis
X <- y
X.var.vec <- y.var.vec
X2 <- w
U <- t
rm(y,y.var.vec,w,t)
####


ptm <- proc.time()

if(!is.null(X)){
  if(ncol(X)!=length(X.var.vec)){
    print("ncol(X) and X.var.vec are different lengths!")
    break
  }
}

if(is.null(X) & is.null(X2)){
  print("NEED TO SPECIFY ONE OR BOTH OF X OR X2")
  break
}

if(is.null(z2)){
  GAMMA <- NULL
}

#Setup
n <- max(nrow(X),nrow(X2))
J <- nclass
Q <- ncol(as.matrix(z))
if(is.null(x)){
  des <- as.matrix(rep(1,n))
} else {
  des <- cbind(1,x)
}
P <- ncol(des)

if(!is.null(X)){
  M <- ncol(X)
}
if(!is.null(X2)){
  X2 <- as.matrix(X2)
  L <- ncol(X2)
}

if(alpha.method=="Cox"){
  S <- n
} else if(alpha.method=="PEM"){
  S <- ncol(t.mat)
}

if(!is.null(z2)){
  Q3 <- ncol(as.matrix(z2))
}


#Repetitions
loglik.best <- -Inf
for(rep in 1:nreps){
print(paste("Rep =",rep))


if(hist==TRUE){
  m.post1 <- m.post2 <- matrix(NA,n,maxiter)
  if(is.null(z2)){
    m.beta.v <- m.Sbeta.v <- matrix(NA,maxiter,Q+(J-1))
  } else {
    m.beta.v <- m.Sbeta.v <- matrix(NA,maxiter,Q+1+Q3)
  }
  m.alpha <- matrix(NA,S,maxiter)
  m.L.exp <- m.L.obs <- m.L.exp.cum <- matrix(NA,maxiter,1)
  m.kappa <- m.Skappa <- matrix(NA,maxiter,P)

  if(!is.null(X)){  
    m.phat1 <- m.phat2 <- matrix(NA,maxiter,M)
  }

  if(!is.null(X2)){
    m.mu1 <- m.mu2 <- matrix(NA,maxiter,L)
    m.sigma.sq <- vector("list",maxiter)
  }
}

#Start values
if(is.null(start.list)){
  beta.start <- matrix(runif(Q)*3,Q,1)
  v.start <- matrix(runif(J-1)*3,1,(J-1))
  alpha.start <- rep(0.01,S)                                 ####!!!!alpha.start
  kappa.mat.start <- matrix(runif((J-1)*P)*3,(J-1),P)
  if(!is.null(X)){
    p.start <- matrix(runif(J*M),J,M)
  } else {
    p.start <- NA
  }
  if(!is.null(X2)){
    mu.start <- rmvnorm(J,apply(as.matrix(X2),2,mean),sigma=cov(as.matrix(X2)*2))    ####!!!! 
    sigma.sq.start <- (runif(1)+0.001)*cov(as.matrix(X2)) / 2    ####!!!! START VALUES UNAFFECTED BY cov.type!
  } else {
    mu.start <- sigma.sq.start <- NA  
  }
  if(!is.null(z2)){
    GAMMA.start <- as.matrix(runif(Q3)/2) ####Maybe too large?
  } else {
    GAMMA.start <- NULL
  }

} else {
  beta.start <- start.list$beta.start
  v.start <- start.list$v.start
  alpha.start <- start.list$alpha.start
  kappa.mat.start <- start.list$kappa.mat.start
  if(!is.null(X)){ 
    p.start <- start.list$p.start
  } else {
    p.start <- NA
  }
  if(!is.null(X2)){
    mu.start <- start.list$mu.start
    sigma.sq.start <- start.list$sigma.sq.start
  } else {
    mu.start <- sigma.sq.start <- NA
  }
  if(!is.null(z2)){
    GAMMA.start <- start.list$GAMMA.start
  } else {
    GAMMA.start <- NULL
  }
}
start.values <- list(p.start,kappa.mat.start,beta.start,v.start,alpha.start,mu.start,sigma.sq.start,GAMMA.start)
names(start.values) <- c("p.start","kappa.mat.start","beta.start","v.start","alpha.start","mu.start","sigma.sq.start","GAMMA.start")
print(paste("GAMMA.start = ",GAMMA.start))


#Initial step
post <- post.fn12(kappa.mat=kappa.mat.start,z=z,alpha=alpha.start,p=p.start,beta=beta.start,v=v.start,des=des,X=X,delta=delta,nclass=nclass,lambda1=lambda1,lambda2=lambda2,X2=X2,mu=mu.start,sigma.sq=sigma.sq.start,alpha.method=alpha.method,ind.mat=ind.mat,t.mat=t.mat,GAMMA=GAMMA.start,z2=z2)
if(is.nan(sum(post))){
  print("lcsm.fn: initial post.fn creates NaNs")
  break
}
loglik <- loglik.fn12(w=post,kappa.mat=kappa.mat.start,alpha=alpha.start,p=p.start,beta=beta.start,v=v.start,des=des,X=X,z=z,delta=delta,nclass=nclass,X2=X2,mu=mu.start,sigma.sq=sigma.sq.start,alpha.method=alpha.method,ind.mat=ind.mat,t.mat=t.mat,GAMMA=GAMMA.start,z2=z2)

if(!is.null(X)){
  phat <- phat.fn12(w=post,X=X,nclass=nclass)
  } else {
     phat <- NA
  }
kappaNR <- kappaNR.fn12(w=post,kappa.mat=kappa.mat.start,des=des,nclass=nclass)
if(kappaNR$Dkappa2.break==1) break

if(!is.null(X2)){
  mu.sigma.sq <- mu.sigma.sq.fn12(post=post,X2=X2,nclass=nclass,cov.type=cov.type,custom.zeros=custom.zeros)
} else {
  mu.sigma.sq <- as.data.frame(NA)
  mu.sigma.sq$mu <- mu.sigma.sq$sigma.sq <- NA
}
betaNR <- betaNR.fn12(w=post,alpha=alpha.start,beta=beta.start,v=v.start,delta=delta,z=z,nclass=nclass,alpha.method=alpha.method,ind.mat=ind.mat,t.mat=t.mat,GAMMA=GAMMA.start,z2=z2)
##############
print("lcsm.fn: initial beta.v = ")
print(betaNR$newbeta.v)
##############

if(!is.null(z2)){
  GAMMA <- betaNR$newbeta.v[(Q+1+1):((Q+1+Q3))]
}

alpha <- alpha.fn12(w=post,beta=betaNR$newbeta.v[1:Q],v=betaNR$newbeta.v[Q+1],z=z,U=U,delta=delta,alpha.method=alpha.method,ind.mat=ind.mat,t.mat=t.mat,GAMMA=GAMMA,z2=z2)


#While loop
diff <- latentgolddiff <- Inf
step = 0

print("Start while loop")
while(diff > tol & step < maxiter){
#for(qq in 1:maxiter){
  #print(paste("(start of while loop)Dkappa2.break=",kappaNR$Dkappa2.break))

  step <- step + 1
  if(printiter==TRUE){
    print(paste("Iteration",step))
  }
  class.switch <- 0 ##Don't think this works

  newpost <- post.fn12(kappa.mat=t(kappaNR$newkappa),z=z,alpha=alpha,p=phat,beta=betaNR$newbeta.v[1:Q],v=betaNR$newbeta.v[Q+1],des=des,X=X,delta=delta,nclass=nclass,lambda1=lambda1,lambda2=lambda2,X2=X2,mu=mu.sigma.sq$mu,sigma.sq=mu.sigma.sq$sigma.sq,alpha.method=alpha.method,ind.mat=ind.mat,t.mat=t.mat,GAMMA=GAMMA,z2=z2)
  if(is.nan(sum(newpost))){
    print("lcsm.fn: post.fn (newpost) creates NaNs")
    break
  }
  ##########################
  print(paste("newpost: min =",min(newpost)," and max = ",max(newpost)))
  ##########################

  newloglik <- loglik.fn12(w=newpost,kappa.mat=t(kappaNR$newkappa),alpha=alpha,p=phat,beta=betaNR$newbeta.v[1:Q],v=betaNR$newbeta.v[Q+1],des=des,X=X,z=z,delta=delta,nclass=nclass,X2=X2,mu=mu.sigma.sq$mu,sigma.sq=mu.sigma.sq$sigma.sq,alpha.method=alpha.method,ind.mat=ind.mat,t.mat=t.mat,GAMMA=GAMMA,z2=z2)

  if(is.nan(newloglik$L.obs)){
      post <- newpost
      loglik <- newloglik
      print(loglik)
      betaNR <- newbetaNR
      alpha <- newalpha
      if(!is.null(X)){
        phat <- newphat
      }
      kappaNR <- newkappaNR
      if(!is.null(X2)){
        mu.sigma.sq <- newmu.sigma.sq
      } 
      break
  }
  if(is.nan(newloglik$L.obs)) break ##Don't think this works properly because break is in the wrong place isn't it?


  if(maxlik.method=="L.obs"){
    diff <- newloglik$L.obs - loglik$L.obs
  }
  if(maxlik.method=="latentgold"){
    diff <- latentgolddiff
  }
  if(maxlik.method=="L.exp"){
    diff <- L.exp.change
  }


  if(!is.null(X)){
    newphat <- phat.fn12(w=newpost,X=X,nclass=nclass)
    #############################
    print("phat = ")
    print(phat)
    #############################
  }
  newkappaNR <- kappaNR.fn12(w=newpost,kappa.mat=t(kappaNR$newkappa),des,nclass)
    #print(paste("(midle of while loop)Dkappa2.break=",newkappaNR$Dkappa2.break))
    #print(newkappaNR)
    if(newkappaNR$Dkappa2.break==1){
      break
    }
  if(!is.null(X2)){
    newmu.sigma.sq <- mu.sigma.sq.fn12(post=newpost,X2=X2,nclass=nclass,cov.type=cov.type,custom.zeros=custom.zeros)
    #############################
    print("newmu.sigma.sq = ")
    print(newmu.sigma.sq)
    #############################
  }

  newbetaNR <- betaNR.fn12(w=newpost,alpha=alpha,beta=betaNR$newbeta.v[1:Q],v=betaNR$newbeta.v[Q+1],delta=delta,z=z,nclass=nclass,alpha.method=alpha.method,ind.mat=ind.mat,t.mat=t.mat,GAMMA=GAMMA,z2=z2)
  ##############
  print("lcsm.fn: while loop newbeta.v = ")
  print(newbetaNR$newbeta.v)
  ##############

  newalpha <- alpha.fn12(w=newpost,beta=newbetaNR$newbeta.v[1:Q],v=newbetaNR$newbeta.v[Q+1],z=z,U=U,delta=delta,alpha.method=alpha.method,ind.mat=ind.mat,t.mat=t.mat,GAMMA=GAMMA,z2=z2)

  #store results
  if(hist==TRUE){
    m.beta.v[step,] <- t(newbetaNR$newbeta.v)
    m.Sbeta.v[step,] <- t(newbetaNR$S)
    m.alpha[,step] <- newalpha
    m.post1[,step] <- newpost[,1]
    m.post2[,step] <- newpost[,2]
    m.L.exp[step] <- newloglik$L.exp
    m.L.obs[step] <- newloglik$L.obs
    if(!is.null(X)){
      m.phat1[step,] <- newphat[1,]
      m.phat2[step,] <- newphat[2,]
    }
    m.kappa[step,] <- t(newkappaNR$newkappa)
    m.Skappa[step,] <- t(newkappaNR$S)
    if(!is.null(X2)){
      m.mu1[step,] <- newmu.sigma.sq$mu[1,]
      m.mu2[step,] <- newmu.sigma.sq$mu[2,]
      m.sigma.sq[[step]] <- newmu.sigma.sq$sigma.sq
    }
  }

  ##latentgold: Calculate differences between old and new
  if(maxlik.method=="latentgold"){
    ##page 62: https://www.statisticalinnovations.com/wp-content/uploads/LGtecnical.pdf
    postdiff <- abs(newpost - post)
    phatdiff <- abs(newphat - phat)
    kappadiff <- abs(newkappaNR$newkappa - kappaNR$newkappa)
    betadiff <- abs(newbetaNR$newbeta.v - betaNR$newbeta.v)
    alphadiff <- abs(newalpha - alpha)
    mudiff <- abs(newmu.sigma.sq$mu - mu.sigma.sq$mu)
    sigma.sqdiff <- abs(newmu.sigma.sq$sigma.sq - mu.sigma.sq$sigma.sq)
    latentgolddiff <- sum(postdiff,phatdiff,kappadiff,betadiff,alphadiff,mudiff,sigma.sqdiff)
  }

  ##replace old with new for re-start
  post <- newpost
  loglik <- newloglik
  if(!is.null(X)){
    phat <- newphat
  }
  kappaNR <- newkappaNR
  if(!is.null(X2)){
    mu.sigma.sq <- newmu.sigma.sq
  }
  betaNR <- newbetaNR
  alpha <- newalpha

  if(!is.null(z2)){
    GAMMA <- betaNR$newbeta.v[(Q+1+1):((Q+1+Q3))]
  } else {
    GAMMA <- NULL
  }
} ##While loop end
print("while loop end")

time <- proc.time() - ptm

if(!is.null(X)){
  if(min(phat,(1-phat)) < boundary.crit){
    boundary <- 1
  } else {
    boundary <- 0
  }
} else {
  boundary <- 0
}


  ##Store params
  params.switch <- NA

  if(!is.null(X)){
    phat.remove.ref <- rep(NA,max(X.var.vec))
    for(gg in 1:max(X.var.vec)){
      phat.remove.ref[gg] <- min(which(X.var.vec==gg))
    }

  phat.noref <- as.matrix(phat[,-c(phat.remove.ref)])

    params.X <- c(phat.noref[1,],phat.noref[2,])
    names(params.X) <- c(rep("phat1",ncol(phat.noref)),rep("phat2",ncol(phat.noref)))
  } else {
    params.X <- NULL
  }

  if(!is.null(X2)){
    L <- ncol(X2)
    if(cov.type=="homo"){
      params.X2 <- c(mu.sigma.sq$mu[1,],mu.sigma.sq$mu[2,],diag(mu.sigma.sq$sigma.sq),
        as.vector(mu.sigma.sq$sigma.sq[lower.tri(mu.sigma.sq$sigma.sq)]))
      names(params.X2) <- c(rep("mu1",L),rep("mu2",L),rep("var",L),rep("covar",(L^2-L)/2))
    }
    if(cov.type=="localind"){
      params.X2 <- c(mu.sigma.sq$mu[1,],mu.sigma.sq$mu[2,],diag(mu.sigma.sq$sigma.sq))
      names(params.X2) <- c(rep("mu1",L),rep("mu2",L),rep("var",L))
    }
    if(cov.type=="custom"){
      params.X2 <- c(mu.sigma.sq$mu[1,],mu.sigma.sq$mu[2,],diag(mu.sigma.sq$sigma.sq),
        as.vector(mu.sigma.sq$sigma.sq[lower.tri(mu.sigma.sq$sigma.sq)][mu.sigma.sq$sigma.sq[lower.tri(mu.sigma.sq$sigma.sq)]!=0]))
      names(params.X2) <- c(rep("mu1",L),rep("mu2",L),rep("var",L),rep("covar",length(as.vector(mu.sigma.sq$sigma.sq[lower.tri(mu.sigma.sq$sigma.sq)][mu.sigma.sq$sigma.sq[lower.tri(mu.sigma.sq$sigma.sq)]!=0]))))
    }
  } else {
    L <- 0
    params.X2 <- NULL
  }

  names(kappaNR$newkappa) <- rep("kappa",P)
  names(betaNR$newbeta.v) <- c(rep("beta",Q),rep("v",(length(betaNR$newbeta.v)-Q)))
  
  params <- as.matrix(c(params.X,kappaNR$newkappa,betaNR$newbeta.v,params.X2))


##Remove NAs from matrices
if(hist==TRUE){
  m.beta.v <- m.beta.v[complete.cases(m.beta.v),]
  m.Sbeta.v <- m.Sbeta.v[complete.cases(m.Sbeta.v),]
  m.alpha <- m.alpha[,complete.cases(t(m.alpha))]
  m.post1 <- m.post1[,complete.cases(t(m.post1))]
  m.post2 <- m.post2[,complete.cases(t(m.post2))]
  m.L.exp <- m.L.exp[complete.cases(m.L.exp),]
  m.L.obs <- m.L.obs[complete.cases(m.L.obs),]

  if(!is.null(X)){
    m.phat1 <- m.phat1[complete.cases(m.phat1),]
    m.phat2 <- m.phat2[complete.cases(m.phat2),]
  } else {
    m.phat1 <- m.phat2 <- NA
  }
  m.kappa <- m.kappa[complete.cases(m.kappa),]
  m.Skappa <- m.Skappa[complete.cases(m.Skappa),]
  if(!is.null(X2)){
    m.mu1 <- m.mu1[complete.cases(m.mu1),]
    m.mu2 <- m.mu2[complete.cases(m.mu2),]
    m.sigma.sq <- m.sigma.sq[!sapply(m.sigma.sq,is.null)] 
  }

  niter <- nrow(m.beta.v)
  nparams <- length(params) + length(alpha)

  AIC <- -2*loglik$L.exp+2*nparams
  BIC <- -2*loglik$L.exp+nparams*log(n)

  return.list <- list(start.values,X,x,z,U,delta,m.beta.v,m.Sbeta.v,m.kappa,m.Skappa,
    m.alpha,m.phat1,m.phat2,m.post1,m.post2,m.L.exp,m.L.obs,m.L.exp.cum,
    post,alpha,phat,loglik,kappaNR,betaNR,step,time,boundary,params,params.switch,lambda1,lambda2,kappaNR$Dkappa2.break,mu.sigma.sq,alpha.method,niter,X2,cov.type,X.var.vec,nparams,AIC,BIC,custom.zeros,z2)
  names(return.list) <- c("start.values","X","x","z","U","delta","m.beta.v","m.Sbeta.v","m.kappa","m.Skappa",
    "m.alpha","m.phat1","m.phat2","m.post1","m.post2","m.L.exp","m.L.obs","m.L.exp.cum",
    "post","alpha","phat","loglik","kappaNR","betaNR","iterations","time","boundary","params","params.switch","lambda1","lambda2","Dkappa2.break","mu.sigma.sq","alpha.method","niter","X2","cov.type","X.var.vec","nparams","AIC","BIC","custom.zeros","z2")  
} else {
  return.list <- list(start.values,X,x,z,U,delta,post,alpha,phat,loglik,kappaNR,betaNR,step,time,boundary,params,lambda1,lambda2,Dkappa2.break,mu.sigma.sq,alpha.method,X2,cov.type,X.var.vec,nparams,AIC,BIC,custom.zeros,z2)
  names(return.list) <- c("start.values","X","x","z","U","delta","post","alpha","phat","loglik","kappaNR","betaNR","iterations","time","boundary","params","params.switch","lambda1","lambda2","Dkappa2.break","mu.sigma.sq","alpha.method","X2","cov.type","X.var.vec","nparams","AIC","BIC","custom.zeros","z2")  
}


print(paste("Current best loglik =",loglik.best))
print(paste("Current rep loglik=",return.list$loglik$L.obs))

  if(is.nan(return.list$loglik$L.obs)){
    print("LOGLIK IS NaN")
    best <- return.list ## but is actually the one that returns NaN, not the best
  } else {
    if(return.list$loglik$L.obs > loglik.best){
      loglik.best <- return.list$loglik$L.obs
      best <- return.list
      print("Best model updated")
      }
  }

print("End of rep")
if(boundary==1){
  print(paste("BOUNDARY SOLUTION OBTAINED [crit: min(phat) <",boundary.crit,"]"))
}
} ##end nreps

return(best)
}

#rm(a)
#a <- lcsm.fn12(X=X4,x=NULL,z=z,U=U,delta=delta,nclass=2,start.list=NULL,hist=TRUE,tol=1e-7,nreps=3,maxiter=1000,printiter=TRUE,boundary.crit=0.001,maxlik.method="L.obs",lambda1=1,lambda2=1,X2=X2,alpha.method="Cox",ind.mat=NULL,t.mat=NULL,cov.type="localind",X.var.vec=X4.var.vec,custom.zeros=NULL,z2=NULL)
#a$params

#rm(a$PEM)
#aPEM <- lcsm.fn12(X=X4,x=NULL,z=z,U=U,delta=delta,nclass=2,start.list=NULL,hist=TRUE,tol=1e-7,nreps=3,maxiter=1000,printiter=TRUE,boundary.crit=0.001,maxlik.method="L.obs",lambda1=1,lambda2=1,X2=X2,alpha.method="PEM",ind.mat=ind.mat,t.mat=t.mat,cov.type="localind",X.var.vec=X4.var.vec,custom.zeros=NULL,z2=NULL)
#aPEM$params
#cbind(a$params,aPEM$params)





#######################
## Interaction model ## ##fix seeds so no class switch between Cox and PEM
####################### ##might need checking is diff data spec

##Random start interaction models
#rm(a3)
#set.seed(98765)
#a3 <- lcsm.fn12(X=X4,x=NULL,z=z,U=U,delta=delta,nclass=2,start.list=NULL,hist=TRUE,tol=1e-7,nreps=3,maxiter=1000,printiter=TRUE,boundary.crit=0.001,maxlik.method="L.obs",lambda1=1,lambda2=1,X2=X2,alpha.method="Cox",ind.mat=NULL,t.mat=NULL,cov.type="localind",X.var.vec=X4.var.vec,custom.zeros=NULL,z2=z)
#a3$params
#a3$betaNR$newbeta.v
#a3$loglik$L.obs
#cbind(a3$betaNR$newbeta.v,c(a$betaNR$newbeta.v,NA),c(beta.true,v.true,GAMMA.true))

#a.se <- lcsm.se.fn11(lcsm.fit=a,inner.method="EM")
#a3.se <- lcsm.se.fn12(lcsm.fit=a3,inner.method="EM")


#rm(a3PEM)
#set.seed(9876)
#a3PEM <- lcsm.fn12(X=X4,x=NULL,z=z,U=U,delta=delta,nclass=2,start.list=NULL,hist=TRUE,tol=1e-7,nreps=3,maxiter=1000,printiter=TRUE,boundary.crit=0.001,maxlik.method="L.obs",lambda1=1,lambda2=1,X2=X2,alpha.method="PEM",ind.mat=ind.mat,t.mat=t.mat,cov.type="localind",X.var.vec=X4.var.vec,custom.zeros=NULL,z2=z)
#a3PEM$params
#a3PEM$betaNR$newbeta.v
#a3PEM$loglik$L.obs

#cbind(a3$betaNR$newbeta.v,a3PEM$betaNR$newbeta.v)


##########################################
## try using maxlik.method="latentgold" ##
##########################################

#rm(a3PEM2)
#a3PEM2 <- lcsm.fn12(X=X4,x=NULL,z=z,U=U,delta=delta,nclass=2,start.list=NULL,hist=TRUE,tol=1e-7,nreps=3,maxiter=1000,printiter=TRUE,boundary.crit=0.001,maxlik.method="latentgold",lambda1=1,lambda2=1,X2=X2,alpha.method="PEM",ind.mat=ind.mat,t.mat=t.mat,cov.type="localind",X.var.vec=X4.var.vec,custom.zeros=NULL,z2=z)
#max(abs(a3PEM$params - a3PEM2$params))


#######################################
## try PEM with perfect start values ##
#######################################

#rm(a3.PEM.start.list)
#a3.PEM.start.list <- list(a3PEM$betaNR$newbeta.v[1],
#  a3PEM$betaNR$newbeta.v[2],
#  a3PEM$alpha,
#  a3PEM$kappaNR$newkappa,
#  a3PEM$phat,
#  a3PEM$mu.sigma.sq$mu,
#  a3PEM$mu.sigma.sq$sigma.sq,
#  a3PEM$betaNR$newbeta.v[3]
#  )
#names(a3.PEM.start.list) <- c("beta.start",
#  "v.start",
#  "alpha.start",
#  "kappa.mat.start",
#  "p.start",
#  "mu.start",
#  "sigma.sq.start",
#  "GAMMA.start"
#  )
#a3.PEM.start.list


#rm(a3PEM3)
#a3PEM3 <- lcsm.fn12(X=X4,x=NULL,z=z,U=U,delta=delta,nclass=2,start.list=a3.PEM.start.list,hist=TRUE,tol=1e-7,nreps=1,maxiter=1000,printiter=TRUE,boundary.crit=0.001,maxlik.method="L.obs",lambda1=1,lambda2=1,X2=X2,alpha.method="PEM",ind.mat=ind.mat,t.mat=t.mat,cov.type="localind",X.var.vec=X4.var.vec,custom.zeros=NULL,z2=z)
#cbind(a3PEM$params,a3PEM3$params)
#max(abs(a3PEM$params - a3PEM3$params))


######################################
##try Cox with perfect start values ##
######################################

#rm(a3.start.list)
#a3.start.list <- list(a3$betaNR$newbeta.v[1],
#  a3$betaNR$newbeta.v[2],
#  a3$alpha,
#  a3$kappaNR$newkappa,
#  a3$phat,
#  a3$mu.sigma.sq$mu,
#  a3$mu.sigma.sq$sigma.sq,
#  a3$betaNR$newbeta.v[3]
#  )
#names(a3.start.list) <- c("beta.start",
#  "v.start",
#  "alpha.start",
#  "kappa.mat.start",
#  "p.start",
#  "mu.start",
#  "sigma.sq.start",
#  "GAMMA.start"
#  )
#a3.start.list


#rm(a33)
#a33 <- lcsm.fn12(X=X4,x=NULL,z=z,U=U,delta=delta,nclass=2,start.list=a3.start.list,hist=TRUE,tol=1e-7,nreps=1,maxiter=1000,printiter=TRUE,boundary.crit=0.001,maxlik.method="L.obs",lambda1=1,lambda2=1,X2=X2,alpha.method="Cox",ind.mat=ind.mat,t.mat=t.mat,cov.type="localind",X.var.vec=X4.var.vec,custom.zeros=NULL,z2=z)
#cbind(a3$params,a33$params)

#max(abs(a3$params - a33$params))


###################################
## adding variables to z and z2  ##
###################################

#set.seed(111)
#newvar <- rnorm(n)/100


#rm(d)
#d <- lcsm.fn12(X=X4,x=NULL,z=cbind(z,newvar),U=U,delta=delta,nclass=2,start.list=NULL,hist=TRUE,tol=1e-7,nreps=10,maxiter=1000,printiter=TRUE,boundary.crit=0.001,maxlik.method="L.obs",lambda1=1,lambda2=1,X2=X2,alpha.method="Cox",ind.mat=NULL,t.mat=NULL,cov.type="localind",X.var.vec=X4.var.vec,custom.zeros=NULL,z2=cbind(z,newvar))
#d$betaNR$newbeta.v

#rm(dPEM)
#dPEM <- lcsm.fn12(X=X4,x=NULL,z=cbind(z,newvar),U=U,delta=delta,nclass=2,start.list=NULL,hist=TRUE,tol=1e-7,nreps=10,maxiter=1000,printiter=TRUE,boundary.crit=0.001,maxlik.method="L.obs",lambda1=1,lambda2=1,X2=X2,alpha.method="PEM",ind.mat=ind.mat,t.mat=t.mat,cov.type="localind",X.var.vec=X4.var.vec,custom.zeros=NULL,z2=cbind(z,newvar))
#dPEM$betaNR$newbeta.v

#cbind(d$betaNR$newbeta.v,dPEM$betaNR$newbeta.v)


###########################################################################
plot.lcsm.fn <- function(obj,param=c("loglik","kappa","beta","v","p","A"),U,truth=NULL,switch=FALSE,start){
  test <- match.arg(param)
  nsteps <- obj$iterations

  if(switch==TRUE){
    truth$p = rbind(truth$p[2,],truth$p[1,])
    truth$v = -truth$v
    truth$lambda = truth$lambda*exp(beta.true*z+v.true)
  }


  if(param=="loglik"){
    par(mfrow=c(1,2))
    plot(seq(1:nsteps),obj$m.L.exp,type="l",main="Expected log likelihood",xlab="Iteration",ylab="expected log-likelihood")
    plot(seq(1:nsteps),obj$m.L.obs,type="l",main="Observed log likelihood",xlab="Iteration",ylab="observed log-likelihood")
  }

  if(param=="kappa"){
    par(mfrow=c(1,2))
    plot(seq(1:nsteps),as.matrix(obj$m.Skappa)[,1],type="l",ylim=c(min(obj$m.Skappa),max(obj$m.Skappa)),main="Score/gradient for kappa parameter(s)",ylab="Score/gradient",xlab="Iteration")
    if(ncol(as.matrix(obj$m.kappa)) == 2){
      points(seq(1:nsteps),as.matrix(obj$m.Skappa)[,2],type="l",lty=2,ylim=c(min(obj$m.Skappa),max(obj$m.Skappa)))
    }
    abline(h=0,lty=1,lwd=2)
    legend("right",legend=c("Skappa0","Skappa1"),lty=c(1,2))

    plot(seq(1:nsteps),as.matrix(obj$m.kappa)[,1],type="l",ylim=c(min(obj$m.kappa,truth$kappa),max(obj$m.kappa,truth$kappa)),main="Estimates of kappa parameter(s)",ylab="Score/gradient",xlab="Iteration")
    abline(h=truth$kappa[1],lty=1,lwd=2)
    if(ncol(as.matrix(obj$m.kappa)) == 2){
      points(seq(1:nsteps),as.matrix(obj$m.kappa)[,2],type="l",lty=2)
      abline(h=truth$kappa[2],lty=2,lwd=2)
    }
    legend("right",legend=c("kappa0","kappa1","True kappa0","True kappa1"),lty=c(1,2,1,2),lwd=c(1,1,2,2))
  }

  if(param=="beta"){
    par(mfrow=c(1,2))
    plot(seq(1:nsteps),obj$m.Sbeta.v[,1],type="l",ylim=c(min(obj$m.Sbeta.v[,1]),max(obj$m.Sbeta.v[,1])),main="Score/gradient for Beta",ylab="Score/gradient",xlab="Iteration")
    abline(h=0,lty=1,lwd=2)
    legend("right",legend=c("Sbeta"),lty=1)

    plot(seq(1:nsteps),obj$m.beta.v[,1],type="l",ylim=c(min(obj$m.beta.v[,1],truth$beta),max(obj$m.beta.v[,1],truth$beta)),main="Beta estimate",ylab="Score/gradient",xlab="Iteration")
    abline(h=truth$beta,lty=1,lwd=2)
    legend("right",legend=c("beta"),lty=1)
  }

  if(param=="v"){
    par(mfrow=c(1,2))
    plot(seq(1:nsteps),obj$m.Sbeta.v[,2],type="l",ylim=c(min(obj$m.Sbeta.v[,2]),max(obj$m.Sbeta.v[,2])),main="Score/gradient for v",ylab="Score/gradient",xlab="Iteration")
    abline(h=0,lty=1,lwd=2)
    legend("right",legend=c("Sv"),lty=1)

    plot(seq(1:nsteps),obj$m.beta.v[,2],type="l",ylim=c(min(obj$m.beta.v[,2],truth$v),max(obj$m.beta.v[,2],truth$v)),main="v estimate",ylab="Score/gradient",xlab="Iteration")
    abline(h=truth$v,lty=1,lwd=2)
    legend("right",legend=c("v"),lty=1)
  }

  if(param=="p"){
    par(mfrow=c(1,ncol(obj$phat)))
    for(h in 1:ncol(obj$phat)){
      plot(seq(1:nsteps),obj$m.phat1[,h],lty=1,type="l",ylim=c(0,1),main=paste("phat",h),xlab="Iteration",ylab="est.")
      points(seq(1:nsteps),obj$m.phat2[,h],lty=2,type="l")
      abline(h=truth$p[1,h],lty=1,lwd=2,col="blue")
      abline(h=truth$p[2,h],lty=2,lwd=2,col="blue")      
    }
    legend("bottomright",legend=c("phat, LC=1","phat, LC=2","True phat for LC=1","True phat for LC=2"),lty=c(1,2,1,2),lwd=c(1,1,2,2),col=c("black","black","blue","blue"))

  }

  if(param=="A"){
    if(switch==TRUE){ print("not displayed correctly when switch=TRUE") }
    m.A <- cumsum(obj$m.alpha[,nsteps])
    par(mfrow=c(1,1))
    plot(U,m.A,type="l",xlab="Time",ylab="Estimated cumulative hazard",main="Final iter: Estimated versus true cumulative hazard")
    points(U,truth$lambda*U,type="l",col="blue",lwd=2)
    legend("topleft",legend=c("Estimated H","True H"),lwd=c(1,2),col=c("black","blue"))
  }

}

###########################################################################

fit.lcsm.fn12 <- function(obj=NULL){

  J <- nclass <- 2
  M <- ncol(obj$X)
  n <- nrow(obj$X)
  observed.df <- as.data.frame.table(table(as.data.frame(obj$X)))

  if(is.null(obj$x)){
    des <- 1
  } else {
    #des <- cbind(rep(1,n),x)
    print("Doesn't work for LCR models")
    break
  }

  class2p <- exp(des %*% obj$kappaNR$newkappa) / (1 + exp(des %*% obj$kappaNR$newkappa)) 
  class1p <- 1 - class2p

  observed <- observed.df

  for(m in 1:M){
    observed[,m] <- as.numeric(observed.df[,m])-1
  }

  observed2 <- observed[,1:M]

  for(var in 1:max(obj$X.var.vec)){
    observed2 <- observed2[rowSums(observed2[,obj$X.var.vec==var])==1,]
  }

  expected <- matrix(NA,nrow(observed2),1)
  expected <- n*(class1p*apply(as.matrix(observed2),1,function(row) prod(obj$phat[1,]^row))+
				 class2p*apply(as.matrix(observed2),1,function(row) prod(obj$phat[2,]^row)))

  freq <- as.matrix(observed[rownames(observed) %in% rownames(observed2),M+1])
  colnames(freq) <- c("Freq")

  pearson.resid <- (freq-expected) / sqrt(expected)
  colnames(pearson.resid) <- c("pearson.resid")

  fit.table <- cbind(observed2,freq,expected,pearson.resid)
  fit.chisq <- sum((fit.table$Freq-fit.table$expected)^2 / fit.table$expected)
  fit.df <- length(obj$kappaNR$newkappa) + length(obj$phat) - J*length(unique(obj$X.var.vec))  ##Is this right? Should we only count parameters in this part of the model?
  fit.chisq.pval <- pchisq(fit.chisq,fit.df,lower.tail=FALSE)
  fit.Gsq <- 2*sum(fit.table$Freq*log(fit.table$Freq/fit.table$expected))
  fit.Gsq.pval <- pchisq(fit.Gsq,fit.df,lower.tail=FALSE)

  ##sub-tables to test for local independence
  factor.mat <- matrix(NA,n,length(unique(obj$X.var.vec)))
  for(var in 1:length(unique(obj$X.var.vec))){
    factor.mat[,var] <- factor(obj$X[,obj$X.var.vec==var] %*% 1:length(obj$X.var.vec[obj$X.var.vec==var]))
  }
   
  if(length(unique(obj$X.var.vec))==1){ ##Skip subtables when only 1 variable in X
    combinations <- subtable.list <- subtable.chisq <- subtable.Gsq <- subtable.chisq.pval <- subtable.Gsq.pval <- NA
  } else {
    combinations <- combn(1:length(unique(obj$X.var.vec)),2)
    subtable.list <- subtable.chisq <- subtable.Gsq <- subtable.chisq.pval <- subtable.Gsq.pval <- vector("list",ncol(combinations))
    for(col in 1:ncol(combinations)){
      subtable.list[[col]] <- as.data.frame.table(table(factor.mat[,combinations[1,col]],factor.mat[,combinations[2,col]]))

      subtable.list[[col]][,c(4,5)] <- t(apply(subtable.list[[col]],1,function(row) colSums(
        fit.table[observed2[,obj$X.var.vec==combinations[1,col]][,as.numeric(row[1])]==1 & observed2[,obj$X.var.vec==combinations[2,col]][,as.numeric(row[2])]==1,(M+1):(M+2)]
      )))

      names(subtable.list[[col]])[4:5] <- c("freq2","expected")

      subtable.chisq[[col]] <- sum((subtable.list[[col]]$Freq - subtable.list[[col]]$expected)^2/subtable.list[[col]]$expected) 
      subtable.Gsq[[col]] <- 2*sum(subtable.list[[col]]$Freq*log(subtable.list[[col]]$Freq/subtable.list[[col]]$expected))
      subtable.chisq.pval[[col]] <- pchisq(subtable.chisq[[col]],(max(as.numeric(subtable.list[[col]][,1]))-1)*(max(as.numeric(subtable.list[[col]][,2]))-1),lower.tail=FALSE)
      subtable.Gsq.pval[[col]] <- pchisq(subtable.Gsq[[col]],(max(as.numeric(subtable.list[[col]][,1]))-1)*(max(as.numeric(subtable.list[[col]][,2]))-1),lower.tail=FALSE)
    }
  }

  return.list <- list(fit.table,fit.chisq,fit.chisq.pval,fit.Gsq,fit.Gsq.pval,combinations,subtable.list,subtable.chisq,subtable.Gsq,subtable.chisq.pval,subtable.Gsq.pval)
  names(return.list) <- c("fit.table","fit.chisq","fit.chisq.pval","fit.Gsq","fit.Gsq.pval","combinations","subtable.list","subtable.chisq","subtable.Gsq","subtable.chisq.pval","subtable.Gsq.pval")
  return(return.list)
}

#################################################################################################


####################################################################
## ind.mat function - creates ind.mat and t.mat for PEM models    ##
## NB: these need to be assigned to the lcsm.fit object for SEs!! ##
####################################################################
 
#part <- c(0,30,60,max(U+1))

ind.mat.fn12  <- function(U=NULL,delta=NULL,part=NULL){

  n <- length(U)

  ## part is a vector a time partitions including 0
  n.part <- length(part)-1
  t.mat <- ind.mat <- matrix(0,n,n.part)

  for(i in 1:n){   ##Jacko had 1:(n-1)?!?
    #print(i)
    tt <- U[i]-part
    tt[which(tt<0)]<-0
    t.mat[i,]<-tt[1:n.part]
    ind.mat[i,max(which(tt>0))]<-1
  }
  #rm(tt)

  for (i in 1:n.part){
    t.mat[which(t.mat[,i]>diff(part)[i]),i] <- diff(part)[i] 
    #ind.mat <- ind.mat*delta
  }

  events <- colSums(delta*ind.mat)
  if(min(colSums(delta*ind.mat))==0){
    print("WARNING: NO EVENTS IN AT LEAST ONE TIME PERIOD. CHANGE PARTITIONS.")
  }

  return.list <- list(ind.mat,t.mat,part,events)
  names(return.list) <- c("ind.mat","t.mat","part","events")
  return(return.list)
}

#ind.mat.list <- ind.mat.fn12(U=byar2$dtime.rand,delta=byar2$died,part=part)
#names(ind.mat.list)
#ind.mat.list$ind.mat
#ind.mat.list$t.mat
#ind.mat.list$part
#ind.mat.list$events


#################################################################################################


####################################################################
## survival residuals function - creates Cox-Snell, Martingale and##
##    Deviance residuals for lcsm models                          ##
####################################################################

lcsm.survresid.fn12 <- function(lcsm.fit){
  ##NB: this only works for alpha.method=Cox and exponential LCSM model (i.e. length(alpha) = 1))

  cox.snell.lc <- xb.lc <- matrix(NA,nrow(lcsm.fit$z),2)

  ##Cox-Snell
  if(lcsm.fit$alpha.method=="Cox"){
    cox.snell.lc[,1] <- cumsum(lcsm.fit$alpha)*exp( cbind(lcsm.fit$z,0,lcsm.fit$z2*0) %*% as.matrix(lcsm.fit$betaNR$newbeta.v) )
    cox.snell.lc[,2] <- cumsum(lcsm.fit$alpha)*exp( cbind(lcsm.fit$z,1,lcsm.fit$z2*1) %*% as.matrix(lcsm.fit$betaNR$newbeta.v) )
  }
  if(lcsm.fit$alpha.method=="PEM" & length(lcsm.fit$alpha)==1){
    cox.snell.lc[,1] <- lcsm.fit$alpha*lcsm.fit$U*exp( cbind(lcsm.fit$z,0,lcsm.fit$z2*0) %*% as.matrix(lcsm.fit$betaNR$newbeta.v) )
    cox.snell.lc[,2] <- lcsm.fit$alpha*lcsm.fit$U*exp( cbind(lcsm.fit$z,1,lcsm.fit$z2*1) %*% as.matrix(lcsm.fit$betaNR$newbeta.v) )
  } 
  cox.snell.lc.comb <- apply(as.matrix(1:nrow(cox.snell.lc)),1,function(gg) cox.snell.lc[gg,] %*% lcsm.fit$post[gg,])

  fitres.cox.snell.lc <- survfit(coxph(Surv(cox.snell.lc.comb,lcsm.fit$delta)~1),type='aalen')
  #Note we're using the residual as the 'time' variable and we estimate the 
  # cumulative hazard of the residuals using a Nelson-Aalen estimate

  #Martingales
  xb.lc[,1] <- cbind(lcsm.fit$z,0,lcsm.fit$z2*0) %*% as.matrix(lcsm.fit$betaNR$newbeta.v)
  xb.lc[,2] <- cbind(lcsm.fit$z,1,lcsm.fit$z2*1) %*% as.matrix(lcsm.fit$betaNR$newbeta.v)
  xb.lc.comb <- apply(as.matrix(1:nrow(xb.lc)),1,function(gg) xb.lc[gg,] %*% lcsm.fit$post[gg,])

  martingale.lc.comb <- lcsm.fit$delta - cox.snell.lc.comb


  ##Deviance
  deviance.lc.comb <- sign(martingale.lc.comb)*( -1*(martingale.lc.comb+lcsm.fit$delta*log(lcsm.fit$delta-martingale.lc.comb)) )^0.5


  ##Plots
  par(mfrow=c(2,2))
  plot(fitres.cox.snell.lc$time,-log(fitres.cox.snell.lc$surv))
  abline(0,1,lty=2)
  #plot(seq(1:length(lcsm.fit$U)),martingale.lc)
  plot(xb.lc.comb,martingale.lc.comb)
  abline(h=0,lty=2)  
  #hist(deviance.lc)
  #plot(1:n,deviance.lc.comb)
  plot(xb.lc.comb,deviance.lc.comb)
  abline(h=0,lty=2)
  par(mfrow=c(1,1))

  return.mat <- cbind(fitres.cox.snell.lc$time,fitres.cox.snell.lc$surv,
                  cox.snell.lc.comb,xb.lc.comb,martingale.lc.comb,deviance.lc.comb)
  colnames(return.mat) <- c("cox.snell.lc.time","cox.snell.lc.surv",
                  "cox.snell.lc.comb","xb.lc.comb","martingale.lc.comb","deviance.lc.comb")

  return(as.data.frame(return.mat))
}


##########################################################################################################
##########################################################################################################

#############
## Entropy ##
#############


entropy.fn12 <- function(lcsm.fit){
  1 - sum( -lcsm.fit$post[,1]*log(lcsm.fit$post[,1]) + -lcsm.fit$post[,2]*log(lcsm.fit$post[,2]) ) / (nrow(lcsm.fit$post)*log(2))
}
#entropy.fn12(LCSM.hunt.custom4.GAMMA)
#entropy.fn12(LCSM.hunt.custom4)
#entropy.fn12(LCSM.hunt.localind)


############################
## LCSM X2 resid function ##
############################
#lcsm.fit <- LCSM.hunt.custom4.PEM.GAMMA

lcsm.X2resid.fn12 <- function(lcsm.fit){

  lcsm.fit$X2
  lcsm.fit$mu.sigma.sq$mu
  lcsm.fit$post

  yhat <- lcsm.fit$post %*% lcsm.fit$mu.sigma.sq$mu
  resid <- lcsm.fit$X2 - yhat

  par(mfrow=c(3,2))
  for(i in 1:ncol(lcsm.fit$X2)){
    hist(resid[,i])
  }

  return.list <- list(yhat,resid)
  names(return.list) <- c("yhat","resid")
  return(return.list)
}



##########################################################################################################
##########################################################################################################


####################
## Jacko PEM code ##
####################

PEM.fn <- function(T,CEN,part,form){

###Setting structures
n<-length(T)
n.part<-length(part)-1
t.mat<-matrix(0,n,n.part)
w.mat<-t.mat


### Defining time , tt , and observation
### i n d i c a t o r , w. mat , in v e c t o r form
for(i in 1:(n-1)){
  tt<-T[i]-part
  tt[which(tt<0)]<-0
  t.mat[i,]<-tt[1:n.part]
  w.mat[i,max(which(tt>0))]<-1
}
for (i in 1:n.part){
  t.mat[which(t.mat[,i]>diff(part)[i]),i]<-diff(part)[i] 
  w.mat<-w.mat*CEN
}

  ### d e f i n i n g s t r u c t u r e s f o r l o g l i n e a r model
  resp<-c(w.mat)
  time<-log(c(t.mat)+0.000001)
  int<-gl(n.part,n)


### S e t t i n g formula
  if(form==~1){
  modmat<-matrix(1,n,1)
  } else {
  modmat<-model.matrix(form)
  }
  n.co<-ncol(modmat)
  coef.mat<-matrix(NA,n*n.part,n.co)
  
  for(i in 1:n.co){
    coef.mat[,i]<-rep(modmat[,i],n.part)
  }
  colnames(coef.mat)<-colnames(modmat)
 
  ### F i t t i n g models
  if(n.part==1){
    loglin<-glm(resp~coef.mat,family="poisson"(link=log),offset=time)
  } else {
    loglin<-glm(resp~-1+int+coef.mat,family="poisson"(link=log),offset=time)
  }
  return(loglin)
}

