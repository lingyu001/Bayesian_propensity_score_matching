rm(list=ls())
gc()

## Locate and read data
# setwd("C:/Users/lul165/Documents/lly/Rfiles")

rsdata<-read.csv("rs_notre_sites_state.csv",na.strings=c(""))


rsdata<-subset(rsdata,shortseg==0,na.strings=c(""))


rsdata<-rsdata[order(rsdata$no),]
rsdata$ints=ifelse(rsdata$no_ints>0,1,0)


## Create simulation sample and result container
## Results including ASMD(absolute standardized mean difference) before and after matching
## and coefficient of treatment, using 1000 samples in the simulation
sampleno<-1000

simul_data<-vector("list",sampleno) 
asmd_mean_bf<-vector("list",sampleno) 
asmd_ub_bf<-vector("list",sampleno)
asmd_mean_psm<-vector("list",sampleno)
asmd_mean_bpsm<-vector("list",sampleno)
asmd_ub_psm<-vector("list",sampleno)
asmd_ub_bpsm<-vector("list",sampleno)
coeff_tre_psm<-vector("list",sampleno)
coeff_tre_bpsm<-vector("list",sampleno)
mcmc_time<-vector("list",sampleno)


## Data and variables used for simulation
  data=list(
    
    lnaadt0308=rsdata$lnaadt0308,
    lnlength=rsdata$lnlength,
    n=length(rsdata$total_crash),
    rhr567=rsdata$rhr567,
    cd=rsdata$curve_density,
    ad=rsdata$access_density,
    dsm=rsdata$d_seg_mi,
    width=rsdata$width,
    hspeed=rsdata$highspeed,
    rspave=rsdata$r_s_pave,
    ints=ifelse(rsdata$no_ints>0,1,0)
  )

## Data used for matching analysis
  mdata<-data.frame(
    lnaadt0308=rsdata$lnaadt0308,
    lnlength=rsdata$lnlength,
    rhr567=rsdata$rhr567,
    curve_density=rsdata$curve_density,
    access_density=rsdata$access_density,
    d_seg_mi=rsdata$d_seg_mi,
    width=rsdata$width,
    highspeed=rsdata$highspeed,
    r_s_pave=rsdata$r_s_pave,
    ints=ifelse(rsdata$no_ints>0,1,0)
  )

## Data simulation
ptm<-proc.time()
  
for(i in 1:1000){
  
  ## Define simulation model
  simul_model="
  
  data{
  
  
  
  for(i in 1:n){
  
  ## Create unobserved covariate, defined as aggressive driving
  ## Unobserved variable is correlated with observed covariates
  aggre_dr[i] ~ dbern(prob[i])
  
  logit(prob[i])<- -1.500-0.1*rhr567[i]-0.01*ad[i]+0.05*width[i]+0.05*hspeed[i]+0.1*ints[i]
  
  
  ## Treatment variable simulation
  ## Treatment variable is correlated with observed and unobserved covariates
  rs_sm[i] ~ dbern(pr[i])
  
  logit(pr[i])<- -1.420-0.316*lnaadt0308[i]+0.120*lnlength[i]-0.041*rhr567[i]-0.007*cd[i]-0.022*ad[i]-0.014*dsm[i]
  +0.003*width[i]-1.085*hspeed[i]+0.346*rspave[i]-0.353*ints[i]+0.5*aggre_dr[i]
  
  ## Outcome variable simulation
  ## Outcome variable is defined as crash frequency, correlated with treatment, observed and unobserved covariates
  ## Outcome variable following negative binomial distribution
  mu[i]<- -5.732+0.709*lnaadt0308[i]+1*lnlength[i]+0.049*rhr567[i]+0.033*cd[i]+0.002*dsm[i]
  -0.149*hspeed[i]-0.024*rspave[i]+0.064*ints[i]+0.1*aggre_dr[i]-0.1*rs_sm[i]
  
  lambda[i]<-exp(mu[i])
  p[i]<-r/(r+lambda[i])
  y[i] ~ dnegbin(p[i],r)
  }
  r ~ dunif(0,50)
  }
  
  model{
  fake <- 0
  }
  
  "
  
  writeLines( simul_model, con="simul_model.txt" )
  
  
  ## Use JAGS for the sampling
  ## Required JAGS installation: http://mcmc-jags.sourceforge.net/
  # install.packages("rjags")
  
  library(rjags)
  
  ## Store outcome, treatment and unobserved variables
  parameters = c("y","rs_sm","aggre_dr")
  
  ## Use JAGS model for to generate one sample
  jagsModel = jags.model( "simul_model.txt" , data=data, n.chains=1 , n.adapt=1000)
  simulated_codasamples = coda.samples(jagsModel,variable.names=parameters,thin=1,n.iter=1)
  
  ## Write out the sample for further analysis
  sd.sp<-as.matrix(simulated_codasamples)
  sd.sp<-t(sd.sp)
  
  aggre_dr<-data.frame(sd.sp[1:length(rsdata$total_crash)])
  colnames(aggre_dr)<-c("aggre_dr")
  
  rs_sm<-data.frame(sd.sp[(length(rsdata$total_crash)+1):(length(rsdata$total_crash)*2)])
  colnames(rs_sm) <-c("rs_sm")
  
  y<-data.frame(sd.sp[(length(rsdata$total_crash)*2+1):(length(rsdata$total_crash)*3)])
  colnames(y) <-c("y")
  
  ## Check the simulated variables
  sum(aggre_dr$aggre_dr)
  sum(rs_sm$rs_sm)
  mean(y$y)
  
  
  ## Add the simulated varaible to analytic data
  mdata$aggre_dr<-aggre_dr$aggre_dr
  mdata$rs<-rs_sm$rs_sm
  mdata$y<-y$y
  
  
  # ## Verify unobserved variable simulation
  # library(MASS)
  # logit_ub<-glm(aggre_dr~rhr567+access_density+width+highspeed+ints,family=binomial(link='logit'),data=mdata)
  # summary(logit_ub)
  # 
  # ## Verify outcome simulation
  # 
  # nbmodel_bf=glm.nb(y~lnaadt0308+lnlength+rhr567+curve_density+d_seg_mi+r_s_pave+highspeed+ints+rs+aggre_dr,data=mdata)
  # summary(nbmodel_bf)
  
  ## Propensity score matching analysis
  
  ## Function to compute matching covariate balance metric ASMD
  
  ## Before matching ASMD
  matching_asmd_bf<-function(var_con_bf,var_tre_bf){
    asmd<-abs(mean(var_tre_bf)-mean(var_con_bf))/((sd(var_tre_bf)^2+sd(var_con_bf)^2)/2)^0.5
    return(asmd)
  }
  
  ## After matching ASMD
  matching_asmd_af<-function(var_con_bf,var_tre_bf,var_con_af,var_tre_af){
    asmd<-abs(mean(var_tre_af)-mean(var_con_af))/((sd(var_tre_bf)^2+sd(var_con_bf)^2)/2)^0.5
    return(asmd)
  }
  
  ## Covariate balance before-matching acording to ASMD values
  
  asmd_lnaadt0308<-matching_asmd_bf(mdata$lnaadt0308[mdata$rs==0],mdata$lnaadt0308[mdata$rs==1])
  asmd_lnlength<-matching_asmd_bf(mdata$lnlength[mdata$rs==0],mdata$lnlength[mdata$rs==1])
  asmd_rhr567<-matching_asmd_bf(mdata$rhr567[mdata$rs==0],mdata$rhr567[mdata$rs==1])
  asmd_cd<-matching_asmd_bf(mdata$curve_density[mdata$rs==0],mdata$curve_density[mdata$rs==1])
  asmd_ad<-matching_asmd_bf(mdata$access_density[mdata$rs==0],mdata$access_density[mdata$rs==1])
  
  asmd_d<-matching_asmd_bf(mdata$d_seg_mi[mdata$rs==0],mdata$d_seg_mi[mdata$rs==1])
  asmd_width<-matching_asmd_bf(mdata$width[mdata$rs==0],mdata$width[mdata$rs==1])
  asmd_hspeed<-matching_asmd_bf(mdata$highspeed[mdata$rs==0],mdata$highspeed[mdata$rs==1])
  asmd_rspave<-matching_asmd_bf(mdata$r_s_pave[mdata$rs==0],mdata$r_s_pave[mdata$rs==1])
  asmd_ints<-matching_asmd_bf(mdata$ints[mdata$rs==0],mdata$ints[mdata$rs==1])
  
  asmd.table.bf<-data.frame(
    asmd_lnaadt0308,asmd_lnlength,asmd_rhr567,asmd_cd,asmd_ad,
    asmd_d,asmd_width,asmd_hspeed,asmd_rspave,asmd_ints
  )
  ## Taking average of ASMD value of all variables
  asmd.mean.bf<-rowMeans(asmd.table.bf)
  asmd.table.bf
  asmd.mean.bf
  
  ##Unobserved covariate ASMD before
  asmd_ub.bf<-matching_asmd_bf(mdata$aggre_dr[mdata$rs==0],mdata$aggre_dr[mdata$rs==1])
  asmd_ub.bf
  
  ## Traditional propensity score matching (PSM)
  ## Omit unobserved variable through the whole process
  
  ## Binary logit model estimating propensit scores
  logit_ps<-glm(rs~lnaadt0308+lnlength+rhr567+curve_density+access_density+d_seg_mi+width+highspeed+r_s_pave+ints
                ,family=binomial(link='logit'),data=mdata)
  summary(logit_ps)
  ## Estimate propensity scores
  mdata$psm_ps<-predict(logit_ps,type="response",data=mdata)
  intercept<- as.numeric(coef(logit_ps)["(Intercept)"])
  
  ## Conduct Nearest-neighbour matching using MatchIt Package
  # install.packages("MatchIt")
  library(MatchIt)
  
  ## Matching specification: 1:2 caliper matching (caliper width = 0.2*sd(propensity scores))
  match.nearest1 <- matchit(rs ~ lnaadt0308+lnlength+rhr567+curve_density+access_density+d_seg_mi+width+highspeed+r_s_pave+ints,
                            #formula for propensity score
                            data = mdata,
                            method="nearest", #matching method
                            distance = "logit", #distance measure; propensity score
                            # distance = rsdata$ps_psm,
                            discard = "both", #observations to discard beyond support
                            replace=TRUE, #matching controls with or without replacement
                            ratio=2, #number of controls per treatment
                            reestimate = FALSE, #re-estimating distance after discarding observations
                            m.order = "largest", #order in which to match treated units
                            caliper = 0.20*sd(mdata$psm_ps), # matching caliper width
                            # mahvars=c("psm_ps") # if two units has same matching distance, how to choose
                            
  )
  summary(match.nearest1)
  # summary(match.nearest1,standardize = TRUE)
  # plot(match.nearest1, type="hist", col="grey")
  
  afpsm <- match.data(match.nearest1)
  summary(afpsm$distance)
  summary(afpsm$psm_ps)
  
  ## Covariable balance after PSM
  
  asmd_lnaadt0308<-matching_asmd_af(mdata$lnaadt0308[mdata$rs==0],mdata$lnaadt0308[mdata$rs==1],afpsm$lnaadt0308[afpsm$rs==0],afpsm$lnaadt0308[afpsm$rs==1])
  asmd_lnlength<-matching_asmd_af(mdata$lnlength[mdata$rs==0],mdata$lnlength[mdata$rs==1],afpsm$lnlength[afpsm$rs==0],afpsm$lnlength[afpsm$rs==1])
  asmd_rhr567<-matching_asmd_af(mdata$rhr567[mdata$rs==0],mdata$rhr567[mdata$rs==1],afpsm$rhr567[afpsm$rs==0],afpsm$rhr567[afpsm$rs==1])
  asmd_cd<-matching_asmd_af(mdata$curve_density[mdata$rs==0],mdata$curve_density[mdata$rs==1],afpsm$curve_density[afpsm$rs==0],afpsm$curve_density[afpsm$rs==1])
  asmd_ad<-matching_asmd_af(mdata$access_density[mdata$rs==0],mdata$access_density[mdata$rs==1],afpsm$access_density[afpsm$rs==0],afpsm$access_density[afpsm$rs==1])
  
  asmd_d<-matching_asmd_af(mdata$d_seg_mi[mdata$rs==0],mdata$d_seg_mi[mdata$rs==1],afpsm$d_seg_mi[afpsm$rs==0],afpsm$d_seg_mi[afpsm$rs==1])
  asmd_width<-matching_asmd_af(mdata$width[mdata$rs==0],mdata$width[mdata$rs==1],afpsm$width[afpsm$rs==0],afpsm$width[afpsm$rs==1])
  asmd_hspeed<-matching_asmd_af(mdata$highspeed[mdata$rs==0],mdata$highspeed[mdata$rs==1],afpsm$highspeed[afpsm$rs==0],afpsm$highspeed[afpsm$rs==1])
  asmd_rspave<-matching_asmd_af(mdata$r_s_pave[mdata$rs==0],mdata$r_s_pave[mdata$rs==1],afpsm$r_s_pave[afpsm$rs==0],afpsm$r_s_pave[afpsm$rs==1])
  asmd_ints<-matching_asmd_af(mdata$ints[mdata$rs==0],mdata$ints[mdata$rs==1],afpsm$ints[afpsm$rs==0],afpsm$ints[afpsm$rs==1])
  
  asmd.table.psm<-data.frame(
    asmd_lnaadt0308,asmd_lnlength,asmd_rhr567,asmd_cd,asmd_ad,
    asmd_d,asmd_width,asmd_hspeed,asmd_rspave,asmd_ints
  )
  
  asmd.mean.psm<-rowMeans(asmd.table.psm)
  asmd.table.psm
  asmd.mean.psm
  
  ## Unobserved covariate ASMD after PSM
  ## Check to see if this variable also got balanced by the matching
  asmd_ub.psm<-matching_asmd_af(mdata$aggre_dr[mdata$rs==0],mdata$aggre_dr[mdata$rs==1],afpsm$aggre_dr[afpsm$rs==0],afpsm$aggre_dr[afpsm$rs==1])
  asmd_ub.psm
  asmd_ub.bf
  
  ## Outcome model - Negative binomial regression
  
  # install.packages("MASS")
  library(MASS)
  nbmodel_psm=glm.nb(y~lnaadt0308+lnlength+rhr567+curve_density+d_seg_mi+r_s_pave+highspeed+ints+rs,data=afpsm)
  # summary(nbmodel_psm)
  
  coefficient<-coef(nbmodel_psm)["rs"]
  
  
  ## Store all the simulation and PSM results
  
  simul_data[[i]]<-cbind(aggre_dr,rs_sm,y) # simulated sample
  asmd_mean_bf[[i]]<-asmd.mean.bf # average ASMD before matching
  asmd_ub_bf[[i]]<-asmd_ub.bf # unobserved variable before matching
  asmd_mean_psm[[i]]<-asmd.mean.psm # average ASMD after PSM
  asmd_ub_psm[[i]]<-asmd_ub.psm # unobserved variable after PSM
  coeff_tre_psm[[i]]<- coef(nbmodel_psm)["rs"] # coefficient of treatment estimated after PSM
  # write.csv(simul_data,"simul_data_test.csv")
  
  ## Bayesian propensity score matching (PSM)
  ## Using JAGS, but package R2jags in R for parallel computation.
  # install.packages("R2jags")
  library(R2jags)
  
  ## Data used for Bayesian modeling
  data=list(
    y=mdata$rs,
    lnaadt0308=mdata$lnaadt0308,
    lnlength=mdata$lnlength,
    n=length(mdata$total_crash),
    rhr567=mdata$rhr567,
    cd=mdata$curve_density,
    ad=mdata$access_density,
    dsm=mdata$d_seg_mi,
    width=mdata$width,
    hspeed=mdata$highspeed,
    rspave=mdata$r_s_pave,
    ints=mdata$ints,
    intercept=intercept
    
  )
  
  ## Define Bayesian propensity score model
  modelString = "

  model{
  
  ## Model specification
  
  for (i in 1:n){
  y[i] ~ dbern(p[i])
  
  logit(p[i]) <- b0+b1*lnaadt0308[i]+b2*lnlength[i]+b3*rhr567[i]+b4*cd[i]+b5*ad[i]+b6*dsm[i]+b7*width[i]+b8*hspeed[i]+b9*rspave[i]+b10*ints[i]
  
   
  }
  
  
  ## Prior setting for each parameter
  ## Use non-informative priors (normal distribution)
  
  b0 ~ dnorm(intercept,1/2^2)
  b1 ~ dnorm(0,1/1^2)
  b2 ~ dnorm(0,1/1^2)
  b3 ~ dnorm(0,1/1^2)
  b4 ~ dnorm(0,1/1^2)
  b5 ~ dnorm(0,1/1^2)
  b6 ~ dnorm(0,1/1^2)
  b7 ~ dnorm(0,1/1^2)
  b8 ~ dnorm(-1,1/1^2)
  b9 ~ dnorm(0,1/1^2)
  b10 ~ dnorm(0,1/1^2)
  
  #examine the posterior
  #lambda_mean<-mean(lambda)
  p_mean <- mean(p)
  
  #examine the likelihood of data
  #y_mean <- mean(y)
  
  }
  "
  writeLines( modelString , con="psmodel.txt" )
  
  ## Parameters to save, estimates of Bayesian propensity scores
  parameters = c("p")
  
  ## Parallel sampling
  ptm_mcmc<-proc.time()
  
  jagsmodel.p=
    jags.parallel(
      data=data, 
      parameters.to.save=parameters, 
      model.file="psmodel.txt",
      n.chains=4, # number of Markov Chains
      n.iter=800, # 80000 number samples in each chain
      n.burnin=50, # 5000 burning steps before save the sample
      n.thin=3, # Thinning steps remove autocorrelation
      n.cluster = 4, # use 4 cores for parallel computation
      DIC=TRUE, working.directory=NULL, 
      jags.seed = 123,
      # refresh = n.iter/50, 
      # progress.bar = "text", digits=5,
      RNGname = c("Wichmann-Hill", "Marsaglia-Multicarry",
                  "Super-Duper", "Mersenne-Twister"),
      # jags.module = c("glm","dic")
      export_obj_names=NULL,
      envir = .GlobalEnv
    )
  nchain=4
  niter=80000
  nburnin=5000
  nthin=3
  samplenum<-(niter-nburnin)/nthin*nchain # print sample number
  durtime_mcmc <- proc.time()-ptm_mcmc
  
  jagsmodel.p.mcmc<-as.mcmc(jagsmodel.p)
  # source("posteriorSummaryStats.R")
  # print(summarizePost(jagsmodel.p.mcmc, filters = c("p")) )
  
  ## Output the sample from mcmc object and write them into analytic data
  ps.sp<-as.matrix(jagsmodel.p.mcmc)
  ps<-as.matrix(colMeans(ps.sp))
  ps<-as.data.frame(ps[-1,])
  colnames(ps)<-c("ps")
  ps$X <- rownames(ps)
  ps$no<-gsub("[^0-9]", "", ps$X) 
  ps<-ps[order(as.numeric(ps$no)),]
  mdata$ps<-as.numeric(ps$ps)
  
  ## Conduct Nearest-neighbour matching for BPSM
  library(MatchIt)
  
  match.nearest2 <- matchit(rs ~ lnaadt0308+lnlength+rhr567+curve_density+access_density+d_seg_mi+width+highspeed+r_s_pave+ints,
                            
                            #formula for propensity score
                            data = mdata,
                            method="nearest", #matching method
                            distance = mdata$ps, #distance measure: Bayesian propensity score
                            # distance = rsdata$ps_psm,
                            discard = "both", #observations to discard beyond support
                            replace=TRUE, #matching controls with or without replacement
                            ratio=2, #number of controls per treatment
                            reestimate = FALSE, #re-estimating distance after discarding observations
                            m.order = "largest", #order in which to match treated units
                            caliper = 0.20*sd(mdata$ps),
                            mahvars=c("ps")
  )
  summary(match.nearest2)
  # summary(match.nearest1,standardize = TRUE)
  # plot(match.nearest1, type="hist", col="grey")
  afbpsm <- match.data(match.nearest2)
  
  summary(afbpsm$distance)
  summary(afbpsm$ps)
  
  ## Covariable balance after BPSM
  
  asmd_lnaadt0308<-matching_asmd_af(mdata$lnaadt0308[mdata$rs==0],mdata$lnaadt0308[mdata$rs==1],afbpsm$lnaadt0308[afbpsm$rs==0],afbpsm$lnaadt0308[afbpsm$rs==1])
  asmd_lnlength<-matching_asmd_af(mdata$lnlength[mdata$rs==0],mdata$lnlength[mdata$rs==1],afbpsm$lnlength[afbpsm$rs==0],afbpsm$lnlength[afbpsm$rs==1])
  asmd_rhr567<-matching_asmd_af(mdata$rhr567[mdata$rs==0],mdata$rhr567[mdata$rs==1],afbpsm$rhr567[afbpsm$rs==0],afbpsm$rhr567[afbpsm$rs==1])
  asmd_cd<-matching_asmd_af(mdata$curve_density[mdata$rs==0],mdata$curve_density[mdata$rs==1],afbpsm$curve_density[afbpsm$rs==0],afbpsm$curve_density[afbpsm$rs==1])
  asmd_ad<-matching_asmd_af(mdata$access_density[mdata$rs==0],mdata$access_density[mdata$rs==1],afbpsm$access_density[afbpsm$rs==0],afbpsm$access_density[afbpsm$rs==1])
  
  asmd_d<-matching_asmd_af(mdata$d_seg_mi[mdata$rs==0],mdata$d_seg_mi[mdata$rs==1],afbpsm$d_seg_mi[afbpsm$rs==0],afbpsm$d_seg_mi[afbpsm$rs==1])
  asmd_width<-matching_asmd_af(mdata$width[mdata$rs==0],mdata$width[mdata$rs==1],afbpsm$width[afbpsm$rs==0],afbpsm$width[afbpsm$rs==1])
  asmd_hspeed<-matching_asmd_af(mdata$highspeed[mdata$rs==0],mdata$highspeed[mdata$rs==1],afbpsm$highspeed[afbpsm$rs==0],afbpsm$highspeed[afbpsm$rs==1])
  asmd_rspave<-matching_asmd_af(mdata$r_s_pave[mdata$rs==0],mdata$r_s_pave[mdata$rs==1],afbpsm$r_s_pave[afbpsm$rs==0],afbpsm$r_s_pave[afbpsm$rs==1])
  asmd_ints<-matching_asmd_af(mdata$ints[mdata$rs==0],mdata$ints[mdata$rs==1],afbpsm$ints[afbpsm$rs==0],afbpsm$ints[afbpsm$rs==1])
  
  asmd.table.bpsm<-data.frame(
    asmd_lnaadt0308,asmd_lnlength,asmd_rhr567,asmd_cd,asmd_ad,
    asmd_d,asmd_width,asmd_hspeed,asmd_rspave,asmd_ints
  )
  
  asmd.mean.bpsm<-rowMeans(asmd.table.bpsm)
  asmd.table.bpsm
  asmd.mean.bpsm
  
  ## Unobserved covariate ASMD after PSM
  ## Check to see if this variable also got balanced by the matching
  asmd_ub.bpsm<-matching_asmd_af(mdata$aggre_dr[mdata$rs==0],mdata$aggre_dr[mdata$rs==1],afbpsm$aggre_dr[afbpsm$rs==0],afbpsm$aggre_dr[afbpsm$rs==1])
  asmd_ub.bpsm
  asmd_ub.bf
  
  ## Outcome model using Negarive binomial regression
  library(MASS)
  nbmodel_bpsm=glm.nb(y~lnaadt0308+lnlength+rhr567+curve_density+d_seg_mi+r_s_pave+highspeed+ints+rs,data=afbpsm)
  # summary(nbmodel_bpsm)
  
  coef(nbmodel_bpsm)["rs"]
  
  
  ## Store BPSM results
  
  asmd_mean_bpsm[[i]]<-asmd.mean.bpsm
  asmd_ub_bpsm[[i]]<-asmd_ub.bpsm
  coeff_tre_bpsm[[i]]<- coef(nbmodel_bpsm)["rs"]
  mcmc_time[[i]]<-durtime_mcmc

}

durtime <- proc.time()-ptm
durtime


## Write out all the saved results
coeff_tre_psm_results<-data.matrix(coeff_tre_psm)
coeff_tre_bpsm_results<-data.matrix(coeff_tre_bpsm)
write.csv(coeff_tre_psm_results,"coeff_tre_psm.csv")
write.csv(coeff_tre_bpsm_results,"coeff_tre_bpsm.csv")

asmd_mean_psm_results<-data.matrix(asmd_mean_psm)
asmd_mean_bpsm_results<-data.matrix(asmd_mean_bpsm)
write.csv(asmd_mean_psm_results,"asmd_mean_psm_results.csv")
write.csv(asmd_mean_bpsm_results,"asmd_mean_bpsm_results.csv")

asmd_mean_bf_results<-data.matrix(asmd_mean_bf)
write.csv(asmd_mean_bf_results,"asmd_mean_bf_results.csv")

asmd_ub_psm_results<-data.matrix(asmd_ub_psm)
asmd_ub_bpsm_results<-data.matrix(asmd_ub_bpsm)
write.csv(asmd_ub_psm_results,"asmd_ub_psm_results.csv")
write.csv(asmd_ub_bpsm_results,"asmd_ub_bpsm_results.csv")


asmd_ub_bf_results<-data.matrix(asmd_ub_bf)
write.csv(asmd_ub_bf_results,"asmd_ub_bf_results.csv")

dir.create("C:/Users/lul165/Documents/lly/Rfiles/simual_data")
for (i in 1:1000){
  data_name <- paste(paste("data_", i, ".csv",sep = ""))
  write.csv(simul_data[i],data_name,dir="C:/Users/lul165/Documents/lly/Rfiles/simual_data")
}


## Plot distribution of treatment and control group for comparing PSM and BPSM

library(ggplot2)

## Distribution of ASMD for unobserved variables
aggre_psm<-ggplot(add,aes(x=add$aggre_psm))+
  geom_density(color="darkblue", fill="lightblue")+
  scale_x_continuous(name="ASMD")+
  scale_y_continuous(name="Density")+
  ggtitle("ASMD of unobserved (PSM)")+
  theme_bw()+theme(plot.title=element_text(size=24,face="bold",hjust=0.5),text=element_text(size=20),legend.position = "bottom")
aggre_psm

aggre_bpsm<-ggplot(add,aes(x=add$aggre_bpsm))+
  geom_density(color="darkred", fill="pink")+
  scale_x_continuous(name="ASMD")+
  scale_y_continuous(name="Density")+
  ggtitle("ASMD of unobserved (BPSM)")+
  theme_bw()+theme(plot.title=element_text(size=24,face="bold",hjust=0.5),text=element_text(size=20),legend.position = "bottom")
aggre_bpsm

## ASMD observed
asmd_mean_psm<-ggplot(add,aes(x=add$asmd_mean_psm))+
  geom_density(color="darkblue", fill="lightblue")+
  scale_x_continuous(name="ASMD")+
  scale_y_continuous(name="Density")+
  ggtitle("ASMD of observed (PSM)")+
  theme_bw()+theme(plot.title=element_text(size=24,face="bold",hjust=0.5),text=element_text(size=20),legend.position = "bottom")
asmd_mean_psm

asmd_mean_bpsm<-ggplot(add,aes(x=add$asmd_mean_bpsm))+
  geom_density(color="darkred", fill="pink")+
  scale_x_continuous(name="ASMD")+
  scale_y_continuous(name="Density")+
  ggtitle("ASMD of observed (BPSM)")+
  theme_bw()+theme(plot.title=element_text(size=24,face="bold",hjust=0.5),text=element_text(size=20),legend.position = "bottom")
asmd_mean_bpsm


## CMF Bias

bias_psm<-ggplot(add,aes(x=add$bias_psm))+
  geom_density(color="darkblue", fill="lightblue")+
  scale_x_continuous(name="CMF estimates - True CMF")+
  scale_y_continuous(name="Density")+
  ggtitle("Bias of CMF(PSM)")+
  theme_bw()+theme(plot.title=element_text(size=24,face="bold",hjust=0.5),text=element_text(size=20),legend.position = "bottom")
bias_psm

bias_bpsm<-ggplot(add,aes(x=add$bias_bpsm))+
  geom_density(color="darkred", fill="pink")+
  scale_x_continuous(name="CMF estimates - True CMF")+
  scale_y_continuous(name="Density")+
  ggtitle("Bias of CMF(BPSM)")+
  theme_bw()+theme(plot.title=element_text(size=24,face="bold",hjust=0.5),text=element_text(size=20),legend.position = "bottom")
bias_bpsm

