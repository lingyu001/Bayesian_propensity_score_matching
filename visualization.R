rm(list=ls())
gc()


## Plot distribution of treatment and control group for comparing PSM and BPSM

library(ggplot2)
library(reshape)

simul_results<-read.csv("simul_results.csv")
head(simul_results)

## Average ASMD density plot
asmd_mean<-simul_results[c("asmd_mean_bf","asmd_mean_psm","asmd_mean_bpsm")]
colnames(asmd_mean)<-c("Before matching","PSM","BPSM")
asmd_mean<-melt(asmd_mean,v.names="asmd")
colnames(asmd_mean)<-c("method","average_asmd")

avg_asmd_plot<-ggplot(asmd_mean,aes(x=asmd_mean$average_asmd,fill=asmd_mean$method))+
  geom_density(position="identity",alpha=0.6)+
  scale_x_continuous(name="Average ASMD")+
  scale_y_continuous(name="Density")+
  scale_fill_manual(values=c("dodgerblue3","green3","lightpink2"))+
  ggtitle("ASMD of observed covariates")+guides(fill=guide_legend(title="Method"))+
  theme_bw()+theme(plot.title=element_text(size=24,face="bold",hjust=0.5),text=element_text(size=20),legend.position = "bottom")
avg_asmd_plot

## Unobserved covariate density plot
asmd_ub<-simul_results[c("asmd_ub_bf","asmd_ub_psm","asmd_ub_bpsm")]
colnames(asmd_ub)<-c("Before matching","PSM","BPSM")
asmd_ub<-melt(asmd_ub,v.names="asmd")
colnames(asmd_ub)<-c("method","aggre_asmd")

asmd_ub_plot<-ggplot(asmd_ub,aes(x=asmd_ub$aggre_asmd,fill=asmd_ub$method))+
  geom_density(position="identity",alpha=0.6)+
  scale_x_continuous(name="ASMD")+
  scale_y_continuous(name="Density")+
  scale_fill_manual(values=c("dodgerblue3","green3","lightpink2"))+
  ggtitle("ASMD of unobserved covariate")+guides(fill=guide_legend(title="Method"))+
  theme_bw()+theme(plot.title=element_text(size=24,face="bold",hjust=0.5),text=element_text(size=20),legend.position = "bottom")
asmd_ub_plot



## CMF Bias
cmf_bias<-simul_results[c("coeff_tre_psm","coeff_tre_bpsm")]
colnames(cmf_bias)<-c("PSM","BPSM")
cmf_bias<-melt(cmf_bias)
colnames(cmf_bias)<-c("method","coeff")
cmf_bias$bias<-exp(cmf_bias$coeff)-exp(-0.1)


bias_plot<-ggplot(cmf_bias,aes(x=cmf_bias$bias,fill=cmf_bias$method))+
  geom_density(position="identity",alpha=0.6)+
  scale_x_continuous(name="CMF estimates - True CMF")+
  scale_y_continuous(name="Density")+
  scale_fill_manual(values=c("green3","lightpink2"))+
  ggtitle("Bias of CMF")+guides(fill=guide_legend(title="Method"))+
  theme_bw()+theme(plot.title=element_text(size=24,face="bold",hjust=0.5),text=element_text(size=20),legend.position = "bottom")
bias_plot

