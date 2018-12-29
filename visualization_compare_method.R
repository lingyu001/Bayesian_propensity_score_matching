rm(list=ls())
gc()





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
