
## Use JAGS for the simulate samples
## Required JAGS installation: http://mcmc-jags.sourceforge.net/

# install.packages("rjags")
require(rjags)

simul_data<-function(simul_model,data,parameters){
  
parameters<-sort(parameters) # sort parameter for output simulation sample purpose

## Use JAGS model for to generate one sample
# jagsModel = jags.model( "simul_model.txt" , data=data, n.chains=1 , n.adapt=1000)
jagsModel = jags.model( simul_model, data=data, n.chains=1 , n.adapt=100)

codasamples = coda.samples(jagsModel,variable.names=parameters,thin=1,n.iter=1)

## Comvert MCMC object to data frame to output the simulated sample

simul_matrix<-as.matrix(codasamples) # convert mcmc object to a matrix
simul_matrix<-t(simul_matrix) 

num_param<-length(parameters)
length_data<-length(data[[1]])
simul_sample<-matrix(nrow = length_data,ncol = num_param)

for (i in 1:num_param){
  
  simul_sample[,i]=simul_matrix[((i-1)*length_data+1):(i*length_data)]

}

colnames(simul_sample)<-parameters
simul_sample<-data.frame(simul_sample)

## Return a data frame contains simulated variables
return(simul_sample)
}


# Using JAGS, but package R2jags in R for Bayesian modeling (parallel computation).


jags_modeling<-function(jagsmodel,data,parameters,nchain,niter,nburnin,nthin,ncores){
  
  nburnin<-as.numeric(nburnin)
  nchain<-as.numeric(nchain)
  niter<-as.numeric(niter)
  nthin<-as.numeric(nthin)
  ncores<-as.numeric(ncores)
  
  ## Parallel sampling
  jagsmodel.p=
    jags.parallel(
      data=data, 
      parameters.to.save=parameters, 
      model.file=jagsmodel,
      n.chains=nchain, # number of Markov Chains
      n.iter=niter, # 80000 number samples in each chain
      n.burnin=nburnin, # 5000 burning steps before save the sample
      n.thin=nthin, # Thinning steps remove autocorrelation
      n.cluster = ncores, # use 4 cores for parallel computation
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
 
  samplenum<-(niter-nburnin)/nthin*nchain # print sample number
  
  ## Summary statistics of Bayesian posterior
  # source("posteriorSummaryStats.R")
  # print(summarizePost(jagsmodel.p.mcmc, filters = parameters) )
  
  ## Output the sample from mcmc object and write them into analytic data
  jagsmodel.p.mcmc<-as.mcmc(jagsmodel.p)
  ps.sp<-as.matrix(jagsmodel.p.mcmc)
  ps<-as.matrix(colMeans(ps.sp))
  ps<-as.data.frame(ps[-1,])
  colnames(ps)<-c("ps")
  ps$X <- rownames(ps)
  ps$no<-gsub("[^0-9]", "", ps$X) 
  ps<-ps[order(as.numeric(ps$no)),]
  
  ## Return a data frame contained the PS estimates
  return(ps)
}