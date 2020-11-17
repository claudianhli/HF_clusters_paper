library(poLCA)
library(dplyr)
library(ggplot2)
library(reshape2)
library(Hmisc)

### Read in data file

df<-readRDS("data.rds")

### Format dataset for poLCA

df_lca <- df %>%
  dplyr::select(CM_CAD, CM_AF, CM_LIVERDISEASE, CM_CANCER, CM_PAD, #comorbidities used to derive the comorbidity clusters
                CM_RENALFAIL, CM_DIABETES, CM_HYPERTENSION,  
                CM_CVA, CM_PEPTIC,  CM_DEPRESSION, CM_ALCOHOL,
                CM_ANEMIA,  CM_OBESITY, CM_COPD, CM_DEMENTIA)

df_lca$CM_AF<-ifelse(df_lca$CM_AF==0 , 1, 2) #1=no comorbidity present; 2= comorbidity present
df_lca$CM_CAD<-ifelse(df_lca$CM_CAD==0 , 1, 2)
df_lca$CM_LIVERDISEASE<-ifelse(df_lca$CM_LIVERDISEASE==0 , 1, 2)
df_lca$CM_CANCER<-ifelse(df_lca$CM_CANCER==0 , 1, 2)
df_lca$CM_PAD<-ifelse(df_lca$CM_PAD==0 , 1, 2)
df_lca$CM_RENALFAIL<-ifelse(df_lca$CM_RENALFAIL==0 , 1, 2)
df_lca$CM_DIABETES<-ifelse(df_lca$CM_DIABETES==0 , 1, 2)
df_lca$CM_HYPERTENSION<-ifelse(df_lca$CM_HYPERTENSION==0 , 1, 2)
df_lca$CM_CVA<-ifelse(df_lca$CM_CVA==0 , 1, 2)
df_lca$CM_PEPTIC<-ifelse(df_lca$CM_PEPTIC==0 , 1, 2)
df_lca$CM_DEPRESSION<-ifelse(df_lca$CM_DEPRESSION==0 , 1, 2)
df_lca$CM_ALCOHOL<-ifelse(df_lca$CM_ALCOHOL==0 , 1, 2)
df_lca$CM_ANEMIA<-ifelse(df_lca$CM_ANEMIA==0 , 1, 2)
df_lca$CM_OBESITY<-ifelse(df_lca$CM_OBESITY==0 , 1, 2)
df_lca$CM_COPD<-ifelse(df_lca$CM_COPD==0, 1, 2)
df_lca$CM_DEMENTIA<-ifelse(df_lca$CM_DEMENTIA==0, 1, 2)

### Run 8 models- from 2 to 9 clusters 

f=cbind(CM_CAD, CM_AF, CM_LIVERDISEASE, CM_CANCER, CM_PAD,
        CM_RENALFAIL, CM_DIABETES, 
        CM_CVA,  CM_DEPRESSION,
        CM_ANEMIA,  CM_OBESITY, CM_COPD) ~ 1

two <- poLCA(f, df_lca, nclass=2, maxiter = 5000)

three <- poLCA(f, df_lca, nclass=3, maxiter = 5000)

four <- poLCA(f, df_lca, nclass=4, maxiter = 5000) 

five <- poLCA(f, df_lca, nclass=5, maxiter = 5000)

six <- poLCA(f, df_lca, nclass=6, maxiter = 5000)

seven <- poLCA(f, df_lca, nclass=7, maxiter = 10000) #Increase no. interations to achieve convergence

eight <- poLCA(f, df_lca, nclass=8, maxiter = 10000)

nine <- poLCA(f, df_lca, nclass=9, maxiter = 12000)

### Increase no. rep for chosen class model to reach global rather than local maximum of the log-likelihood function

five_max<-poLCA(f, df_lca, nclass=5, maxiter = 5000, nrep=100)

five<-five_max

################################################
##### Create statistics table for scree plot ###
################################################

results <- data.frame(Model=c("Model 2"),
                      log_likelihood=two$llik,
                      df = two$resid.df,
                      BIC=two$bic,
                      ABIC=  (-2*two$llik) + ((log((two$N + 2)/24)) * two$npar),
                      CAIC = (-2*two$llik) + two$npar * (1 + log(two$N)), 
                      likelihood_ratio=two$Gsq)

results$Model<-as.integer(results$Model)

results[1,1]<-c("Model 2")
results[2,1]<-c("Model 3")
results[3,1]<-c("Model 4")
results[4,1]<-c("Model 5")
results[5,1]<-c("Model 6")
results[6,1]<-c("Model 7")
results[7,1]<-c("Model 8")
results[8,1]<-c("Model 9")

results[1,2]<-two$llik
results[2,2]<-three$llik
results[3,2]<-four$llik
results[4,2]<-five$llik
results[5,2]<-six$llik
results[6,2]<-seven$llik
results[7,2]<-eight$llik
results[8,2]<-nine$llik

results[1,3]<-two$resid.df
results[2,3]<-three$resid.df
results[3,3]<-four$resid.df
results[4,3]<-five$resid.df
results[5,3]<-six$resid.df
results[6,3]<-seven$resid.df
results[7,3]<-eight$resid.df
results[8,3]<-nine$resid.df

results[1,4]<-two$bic
results[2,4]<-three$bic
results[3,4]<-four$bic
results[4,4]<-five$bic
results[5,4]<-six$bic
results[6,4]<-seven$bic
results[7,4]<-eight$bic
results[8,4]<-nine$bic

results[1,5]<-(-2*two$llik) + ((log((two$N + 2)/24)) * two$npar) #abic
results[2,5]<-(-2*three$llik) + ((log((three$N + 2)/24)) * three$npar)
results[3,5]<-(-2*four$llik) + ((log((four$N + 2)/24)) * four$npar)
results[4,5]<-(-2*five$llik) + ((log((five$N + 2)/24)) * five$npar)
results[5,5]<-(-2*six$llik) + ((log((six$N + 2)/24)) * six$npar)
results[6,5]<-(-2*seven$llik) + ((log((seven$N + 2)/24)) * seven$npar)
results[7,5]<-(-2*eight$llik) + ((log((eight$N + 2)/24)) * eight$npar)
results[8,5]<-(-2*nine$llik) + ((log((nine$N + 2)/24)) * nine$npar)

results[1,6]<- (-2*two$llik) + two$npar * (1 + log(two$N)) #caic
results[2,6]<- (-2*three$llik) + three$npar * (1 + log(three$N))
results[3,6]<- (-2*four$llik) + four$npar * (1 + log(four$N))
results[4,6]<- (-2*five$llik) + five$npar * (1 + log(five$N))
results[5,6]<- (-2*six$llik) + six$npar * (1 + log(six$N))
results[6,6]<- (-2*seven$llik) + seven$npar * (1 + log(seven$N))
results[7,6]<- (-2*eight$llik) + eight$npar * (1 + log(eight$N))
results[8,6]<- (-2*nine$llik) + nine$npar * (1 + log(nine$N))

results[1,7]<-two$Gsq
results[2,7]<-three$Gsq
results[3,7]<-four$Gsq
results[4,7]<-five$Gsq
results[5,7]<-six$Gsq
results[6,7]<-seven$Gsq
results[7,7]<-eight$Gsq
results[8,7]<-nine$Gsq

results$model < - as.factor(results$model) 

### Convert to long format to prepare for plot ###

results2<-tidyr::gather(results,Criterion,Value,4:7)
results2
results2$Criterion[results2$Criterion=="likelihood_ratio"] <- "Likelihood ratio"

### Create scree plot ###

tiff("screeplot.tiff", units="in", width = 10.2, height=7, res=300)
ggplot(results2) + 
  geom_point(aes(x=Model,y=Value),size=3) +
  geom_line(aes(Model, Value, group = 1)) +
  theme_bw()+
  labs(x = "", y="", title = "") + 
  facet_grid(Criterion ~. ,scales = "free") +
  theme_bw(base_size = 16, base_family = "") +   
  theme(panel.grid.major.x = element_blank() ,
        panel.grid.major.y = element_line(colour="grey", size=0.5),
        legend.title = element_text(size = 16, face = 'bold'),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text=  element_text(size=16),
        axis.line = element_line(colour = "black")) # 
dev.off()


#######################################################
### Plot best class solution -   five classes model ###
######################################################

plotmodel5<-reshape2::melt(five$probs)
plotmodel5$L1[plotmodel5$L1=="CM_OBESITY"]<-"Obesity"
plotmodel5$L1[plotmodel5$L1=="CM_DEPRESSION"]<-"Depression"
plotmodel5$L1[plotmodel5$L1=="CM_RENALFAIL"]<-"Renal failure"
plotmodel5$L1[plotmodel5$L1=="CM_LIVERDISEASE"]<-"Liver disease"
plotmodel5$L1[plotmodel5$L1=="CM_CANCER"]<-"Cancer"
plotmodel5$L1[plotmodel5$L1=="CM_ANEMIA"]<-"Anemia"
plotmodel5$L1[plotmodel5$L1=="CM_COPD"]<-"COPD"
plotmodel5$L1[plotmodel5$L1=="CM_PAD"]<-"PAD"
plotmodel5$L1[plotmodel5$L1=="CM_CVA"]<-"CVA"
plotmodel5$L1[plotmodel5$L1=="CM_AF"]<-"AF"
plotmodel5$L1[plotmodel5$L1=="CM_DIABETES"]<-"Diabetes"
plotmodel5$L1[plotmodel5$L1=="CM_CAD"]<-"CAD"
plotmodel5$model<-"5-class model"

###  Plot 5 class model ###

tiff(file="5class_bar.tiff", units="in", width=10, height=6, res=300)
ggplot(plotmodel5, aes(x=L1, y=value, fill=Var2))+
  geom_bar(stat="identity", position="stack") +
  facet_wrap(Var1 ~ .) +    
  scale_fill_brewer(type="seq", palette="Greys")+
  theme_bw() +
  theme(axis.text.x=element_text(angle=90))+
  labs(x="Comorbidity", y="Probability")
dev.off()

#########################################################################
### Save posterior probabilities to data frame of individual patients ###
#########################################################################

### Population shares of classes
round(colMeans(five$posterior)*100, 2)

posteriors <- data.frame(five$posterior, #matrix of probabilities for each respondent and each class
                         predclass=five$predclass) #Vector of predicted class membership for each patient done by modal assignment

df<-cbind(df, posteriors)

saveRDS(df, "dflca.rds")

########################################################################
### Run model with "hypertension" added #################################
########################################################################

f_sens=cbind(CM_CAD, CM_AF, CM_LIVERDISEASE, CM_CANCER, CM_PAD,
        CM_RENALFAIL, CM_DIABETES, 
        CM_CVA,  CM_DEPRESSION,
        CM_ANEMIA,  CM_OBESITY, CM_COPD, CM_HYPERTENSION) ~ 1

five_sens<-poLCA(f, df_lca, nclass=5, maxiter = 5000, nrep=100)
#BIC: 4374723
#AIC: 4373987

five
#BIC: 4268117
#AIC: 4267434


######### END LCA #####################


