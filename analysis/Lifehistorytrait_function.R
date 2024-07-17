lifeTimeRepEvents <- function(matU, matF, startLife = 1) {
  
  uDim = dim(matU)[1]
  surv = colSums(matU)
  repLifeStages = colSums(matF)
  repLifeStages[which(repLifeStages>0)] <- 1
  
  if(missing(matF) | missing(matU)) {stop('matU or matF missing')}
  if(sum(matF, na.rm=T)==0) {stop('matF contains only 0 values')}
  
  # Probability of survival to first reprod event
  Uprime <- matU
  Uprime[, which(repLifeStages==1)] <- 0
  Mprime = matrix(0,2,uDim)
  for (p in 1:uDim[1]) {
    if (repLifeStages[p]==1) {Mprime[2,p] = 1} else {
      Mprime[1,p] = 1-surv[p]
    }
  }
  Bprime = Mprime%*%(MASS::ginv(diag(uDim)-Uprime))
  pRep = Bprime[2,startLife]
  
  out = data.frame(pRep = pRep)
  
  # first derive age-trajectories of survivorship (lx) and fecundity (mx)
  lx <- Rage::mpm_to_lx(matU, start = startLife)
  mx <- Rage::mpm_to_mx(matU, matF, start = startLife)
  # out$lx<-Rage::mpm_to_lx(matU, start = startLife)
  # out$mx <- Rage::mpm_to_mx(matU, matF, start = startLife)
  
  # Age at first reproduction (La; Caswell 2001, p 124)
  D = diag(c(Bprime[2,]))
  Uprimecond = D%*%Uprime%*%MASS::ginv(D)
  expTimeReprod = colSums(MASS::ginv(diag(uDim)-Uprimecond))
  La = expTimeReprod[startLife]
  
  out$La = La
  
  # Mean life expectancy conditional on 
  # entering the life cycle in the first reproductive stage
  firstRepLifeStage = min(which(repLifeStages==1))
  N = solve(diag(uDim[1])-matU)
  meanRepLifeExpectancy = colSums(N)[firstRepLifeStage]
  
  out$meanRepLifeExpectancy = meanRepLifeExpectancy
  
  # Life expectancy from mean maturity
  remainingMatureLifeExpectancy = colSums(N)[startLife]-La
  
  out$remainingMatureLifeExpectancy = remainingMatureLifeExpectancy
  
  # mean life expectancy
  out$meanelexp<-Rage::life_expect_mean(matU, start = startLife)
  #  Variance of life expectancy 
  out$variancelexp<-Rage::life_expect_var(matU, start = startLife)
  
  ##longevity
  out$longevity<-Rage::longevity(matU, start = startLife)
  
  # Generation time 
  matA=matU + matF
  N <- solve(diag(dim(matU)[1])-matU)
  R <- matF %*% N
  out$Ro <- lambda(R)
  out$lambda <- lambda(matA)
  out$generation.time <- log(out$Ro)/log(out$lambda)
  
  # Demetrius' entropy
  out$entropyd<-Rage::entropy_d(lx, mx)
  # calculate Keyfitz' entropy
  out$entropyk<-Rage::entropy_k(lx)
  # shape of survival/mortality trajectory
  out$ssmt<-Rage::shape_surv(lx) 
  # shape of fecundity trajectory
  out$sft<-Rage::shape_rep(lx)
  #Spread of reproduction:  Gini index function `Gini` from the `ineq` package.
  ssmt<-Rage::shape_surv(lx) 
  sft<-Rage::shape_rep(lx)
  out$gini<-ineq::Gini(ssmt*sft,corr = F)
  return(out)
}

