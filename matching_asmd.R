
## Function to compute matching covariate balance metric ASMD

## Before matching ASMD
matching_asmd_bf<-function(var_tre_bf,var_con_bf){
  
  asmd<-abs(mean(var_tre_bf)-mean(var_con_bf))/((sd(var_tre_bf)^2+sd(var_con_bf)^2)/2)^0.5
  return(asmd)
}

## After matching ASMD
matching_asmd_af<-function(var_tre_bf,var_con_bf,var_tre_af,var_con_af){
  asmd<-abs(mean(var_tre_af)-mean(var_con_af))/((sd(var_tre_bf)^2+sd(var_con_bf)^2)/2)^0.5
  return(asmd)
}