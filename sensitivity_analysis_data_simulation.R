rm(list=ls())
gc()

## Locate and read data
setwd("C:/Users/liuxi/Google Drive/Dissertation/Github project/sensitivity_analysis_data_simulation")

rsdata<-read.csv("rs_notre_sites_state.csv",na.strings=c(""))
rsdata<-subset(rsdata,shortseg==0,na.strings=c(""))
rsdata<-rsdata[order(rsdata$no),]
rsdata$ints=ifelse(rsdata$no_ints>0,1,0)

## Data used for analysis
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
num_cov<-length(mdata)

## Create simulation sample and result container
## Results including ASMD(absolute standardized mean difference) before and after matching
## and coefficient of treatment, using 1000 samples in the simulation
sampleno<-1000

simul_data<- # Simulated data container
asmd_mean_bf<-asmd_ub_bf<-asmd_mean_psm<-asmd_mean_bpsm<-asmd_ub_psm<-asmd_ub_bpsm<- # Covariate balance metric container
coeff_tre_psm<-coeff_tre_bpsm<- # Treatment effect estimates container
mcmc_time<-vector("list",sampleno) # MCMC time container


## Data used for simulation
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
  rs[i] ~ dbern(pr[i])
  
  logit(pr[i])<- -1.420-0.316*lnaadt0308[i]+0.120*lnlength[i]-0.041*rhr567[i]-0.007*cd[i]-0.022*ad[i]-0.014*dsm[i]
  +0.003*width[i]-1.085*hspeed[i]+0.346*rspave[i]-0.353*ints[i]+0.5*aggre_dr[i]
  
  ## Outcome variable simulation
  ## Outcome variable is defined as crash frequency, correlated with treatment, observed and unobserved covariates
  ## Outcome variable following negative binomial distribution
  mu[i]<- -5.732+0.709*lnaadt0308[i]+1*lnlength[i]+0.049*rhr567[i]+0.033*cd[i]+0.002*dsm[i]
  -0.149*hspeed[i]-0.024*rspave[i]+0.064*ints[i]+0.1*aggre_dr[i]-0.1*rs[i]
  
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

## Define Bayesian propensity score model
modelString = "

model{

## Model specification

for (i in 1:n){
rs[i] ~ dbern(p[i])

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
#y_mean <- mean(rs)

}
"
writeLines( modelString , con="psmodel.txt" )


## Data simulation
ptm<-proc.time()
  
for(i in 1:1000){
  
 
  
  ## Define parameters of interest in the simulations
  parameters = c("y","rs","aggre_dr")
  
  ## Simulation
  source("jags_simulation.R")
  simul_sample<-simul_data("simul_model.txt",data,parameters)
  
  summary(simul_sample)
  
  ## Add the simulated varaible to analytic data
  mdata[names(simul_sample)]<-simul_sample
 
  summary(mdata)
  
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
  
  ## Produce covariate balance metric ASMD
  
  ## Before matching ASMD
  
  asmd.table.bf<-matrix(nrow=1,ncol=num_cov)
  tre_data_bf<-subset(mdata,mdata$rs==1)
  con_data_bf<-subset(mdata,mdata$rs==0)
  
  source("matching_asmd.R")
  for (i in 1:num_cov){
    
    asmd.table.bf[,i]<-matching_asmd_bf(tre_data_bf[[i]],con_data_bf[[i]])
  }
  asmd.mean.bf<-mean(asmd.table.bf)

  ##Unobserved covariate ASMD before matching
  asmd_ub.bf<-matching_asmd_bf(tre_data_bf$aggre_dr,con_data_bf$aggre_dr)
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
  match.nearest.psm <- matchit(rs ~ lnaadt0308+lnlength+rhr567+curve_density+access_density+d_seg_mi+width+highspeed+r_s_pave+ints,
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
  summary(match.nearest.psm)
  # summary(match.nearest.psm,standardize = TRUE)
  # plot(match.nearest.psm, type="hist", col="grey")
  
  afpsm <- match.data(match.nearest.psm)
  summary(afpsm$distance)
  summary(afpsm$psm_ps)
  
  ## Covariable balance after PSM
  asmd.table.psm<-matrix(nrow=1,ncol=num_cov)
  tre_data_psm<-subset(afpsm,afpsm$rs==1)
  con_data_psm<-subset(afpsm,afpsm$rs==0)
  
  for (i in 1:num_cov){
    
    asmd.table.psm[,i]<-matching_asmd_af(tre_data_bf[[i]],con_data_bf[[i]],
                                         tre_data_psm[[i]],con_data_psm[[i]])
  }
  asmd.mean.psm<-mean(asmd.table.psm)
  
  ##Unobserved covariate ASMD after PSM
  asmd_ub.psm<-matching_asmd_af(tre_data_bf$aggre_dr,con_data_bf$aggre_dr,
                                tre_data_psm$aggre_dr,con_data_psm$aggre_dr)
  asmd_ub.psm
  
  
  ## Outcome model - Negative binomial regression
  
  # install.packages("MASS")
  library(MASS)
  nbmodel_psm=glm.nb(y~lnaadt0308+lnlength+rhr567+curve_density+d_seg_mi+r_s_pave+highspeed+ints+rs,data=afpsm)
  # summary(nbmodel_psm)
  
  ## Store all the simulation and PSM results
  
  simul_data[[i]]<-cbind(aggre_dr,rs_sm,y) # simulated sample
  asmd_mean_bf[[i]]<-asmd.mean.bf # average ASMD before matching
  asmd_ub_bf[[i]]<-asmd_ub.bf # unobserved variable before matching
  asmd_mean_psm[[i]]<-asmd.mean.psm # average ASMD after PSM
  asmd_ub_psm[[i]]<-asmd_ub.psm # unobserved variable after PSM
  coeff_tre_psm[[i]]<- coef(nbmodel_psm)["rs"] # coefficient of treatment estimated after PSM
  # write.csv(simul_data,"simul_data_test.csv")
  
  ## Bayesian propensity score matching (PSM)


  ## Data used for Bayesian modeling
  bpsmdata=list(
    rs=mdata$rs,
    lnaadt0308=mdata$lnaadt0308,
    lnlength=mdata$lnlength,
    n=length(mdata$lnaadt0308),
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
  
    ## Parameters to save, estimates of Bayesian propensity scores
  parameters.bpsm = c("p")
  
  ## Bayesian modeling using MCMC simulating
  nchain=ncore=4
  niter=800
  nburnin=50
  nthin=3
  samplenum<-(niter-nburnin)/nthin*nchain # print sample number
  
  source("jags_simulation.R")
  
  ptm_mcmc<-proc.time()
  
  ps<-jags_modeling("psmodel.txt",bpsmdata,parameters.bpsm,nchain,niter,nburnin,nthin,ncore)
  
  durtime_mcmc <- proc.time()-ptm_mcmc
  
  ## Add Bayesian PS to the analytic data
  mdata$ps<-as.numeric(ps$ps)
  
  ## Conduct Nearest-neighbour matching for BPSM
  library(MatchIt)
  
  match.nearest.bpsm <- matchit(rs ~ lnaadt0308+lnlength+rhr567+curve_density+access_density+d_seg_mi+width+highspeed+r_s_pave+ints,
                            
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
                            # mahvars=c("ps")
  )
  summary(match.nearest.bpsm)
  # summary(match.nearest.bpsm,standardize = TRUE)
  # plot(match.nearest.bpsm, type="hist", col="grey")
  afbpsm <- match.data(match.nearest.bpsm)
  
  summary(afbpsm$distance)
  summary(afbpsm$ps)
  
  ## Covariable balance after BPSM
  source("matching_asmd.R")
  
  asmd.table.bpsm<-matrix(nrow=1,ncol=num_cov)
  tre_data_bpsm<-subset(afbpsm,afbpsm$rs==1)
  con_data_bpsm<-subset(afbpsm,afbpsm$rs==0)
  
  for (i in 1:num_cov){
    
    asmd.table.bpsm[,i]<-matching_asmd_af(tre_data_bf[[i]],con_data_bf[[i]],
                                         tre_data_bpsm[[i]],con_data_bpsm[[i]])
  }
  asmd.mean.bpsm<-mean(asmd.table.bpsm)
  
  ##Unobserved covariate ASMD after bpsm
  asmd_ub.bpsm<-matching_asmd_af(tre_data_bf$aggre_dr,con_data_bf$aggre_dr,
                                tre_data_bpsm$aggre_dr,con_data_bpsm$aggre_dr)
  asmd_ub.bpsm
  
  
  ## Outcome model using Negarive binomial regression
  library(MASS)
  nbmodel_bpsm=glm.nb(y~lnaadt0308+lnlength+rhr567+curve_density+d_seg_mi+r_s_pave+highspeed+ints+rs,data=afbpsm)
  # summary(nbmodel_bpsm)
  
  # coef(nbmodel_bpsm)["rs"]
  
  
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

