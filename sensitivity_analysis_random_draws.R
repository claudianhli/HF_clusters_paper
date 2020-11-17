####################################################################
########  Sensitivity analysis - Random draws ######################
####################################################################
df<-readRDS("dflca.rds")

library(nnet)
library(mice)

#### Create 20 simulated datasets, sampling from the matrix of probabilities for each respondent, per each class 

sink('lc.sens.txt')
simul.classes <- raply(20, # nr simulations
                       aaply(df[, paste0('X', 1:5)], 1, 
                             sample, x=1:5, size=1, replace=FALSE, # sampling parameter
                             .expand = FALSE))

### Run 20 Cox regression models with the imputed datasets
simul.res <- alply(simul.classes, 1, function(assign.class) {
  df$assign.class <- as.factor(assign.class)
  
  sens.model<-coxph(Surv(time_to_event, admitted) ~  AGE + SEX + 
                      + CM_DEMENTIA +  CM_PEPTIC + CM_ALCOHOL + CM_HYPERTENSION
                      + RACE + EDUCATION  
                      + COM  # Insurance status (Commercial vs. Medicare Advantage)
                      + IP_DIAGNOSIS #Place of diagnosis (Inpatient vs. Outpatient)
                      + CARDIOSELECTIVE_BASELINE + NONCARDIOSELECTIVE_BASELINE + ACEI_ARB_BASELINE #medications at baseline
                      + MRA_BASELINE + THIAZIDE_BASELINE + POTASSIUMSPARING_BASELINE + LOOP_BASELINE #medications at baseline
                      + assign.class 
                    , data=df)
  sens.model
})

simul.se <- laply(simul.res, function(x) {sqrt(diag(x$var))}, .drop = TRUE)
simul.var <- laply(simul.res, function(x) {diag(x$var)}, .drop = TRUE)
simul.mean <- laply(simul.res, function(x) {x$coefficients}, .drop = TRUE)

## Rubin's rules  to combine estimates from the 20 Cox models
D <- nrow(simul.mean)
theta.bars <- colMeans(simul.mean)
V.bars <- colMeans(simul.var) # within-imputation var
B <- colMeans(aaply(simul.mean, 1, '-', theta.bars)^2) # between-imputation var
T <- V.bars + (1 - (1/D))*B
names(V.bars) <- names(simul.res[[1]]$coefficients)
names(B) <- names(V.bars)
names(T) <- names(B)
V.bars <- as.matrix(V.bars)
B <- as.matrix(B)
T <- as.matrix(T)
print('within imputation var:')
print(V.bars)
print('between imputation var:')
print(B)
sqrt(B)

simul.se
sink()


### Save results from combined results 

res<-as.data.frame(summary(res.pooled))

res$lrCI<-res$estimate-(1.96*res$std.error)
res$upCI<-res$estimate+(1.96*res$std.error)
res$HR<-exp(res$estimate)
res$lowerCI_HR<-exp(res$lrCI)
res$upperCI_HR <- exp(res$upCI)

sink('imputed_results_combined.txt')
res
sink()


##########################################
###### Save imputed datasets if needed ###
##########################################

sink('imputed_results20.txt')
print('1st imputed datset:')
simul.res[[1]]
exp(coef(simul.res[[1]]))
exp(confint(simul.res[[1]]))

print('2nd imputed datset:')
simul.res[[2]]
exp(coef(simul.res[[2]]))
exp(confint(simul.res[[2]]))


print('3rd imputed datset:')
simul.res[[3]]
exp(coef(simul.res[[3]]))
exp(confint(simul.res[[3]]))

print('4th imputed datset:')
simul.res[[4]]
exp(coef(simul.res[[4]]))
exp(confint(simul.res[[4]]))

print('5th imputed datset:')
simul.res[[5]]
exp(coef(simul.res[[5]]))
exp(confint(simul.res[[5]]))

print('6th imputed datset:')
simul.res[[6]]
exp(coef(simul.res[[6]]))
exp(confint(simul.res[[6]]))

print('7th imputed datset:')
simul.res[[7]]
exp(coef(simul.res[[7]]))
exp(confint(simul.res[[7]]))

print('8th imputed datset:')
simul.res[[8]]
exp(coef(simul.res[[8]]))
exp(confint(simul.res[[8]]))

print('9th imputed datset:')
simul.res[[9]]
exp(coef(simul.res[[9]]))
exp(confint(simul.res[[9]]))

print('10th imputed datset:')
simul.res[[10]]
exp(coef(simul.res[[10]]))
exp(confint(simul.res[[10]]))

print('11th imputed datset:')
simul.res[[11]]
exp(coef(simul.res[[11]]))
exp(confint(simul.res[[11]]))

print('12th imputed datset:')
simul.res[[12]]
exp(coef(simul.res[[12]]))
exp(confint(simul.res[[12]]))

print('13th imputed datset:')
simul.res[[13]]
exp(coef(simul.res[[13]]))
exp(confint(simul.res[[13]]))

print('14th imputed datset:')
simul.res[[14]]
exp(coef(simul.res[[14]]))
exp(confint(simul.res[[14]]))

print('15th imputed datset:')
simul.res[[15]]
exp(coef(simul.res[[15]]))
exp(confint(simul.res[[15]]))

print('16th imputed datset:')
simul.res[[16]]
exp(coef(simul.res[[16]]))
exp(confint(simul.res[[16]]))

print('17th imputed datset:')
simul.res[[17]]
exp(coef(simul.res[[17]]))
exp(confint(simul.res[[17]]))

print('18th imputed datset:')
simul.res[[18]]
exp(coef(simul.res[[18]]))
exp(confint(simul.res[[18]]))

print('19th imputed datset:')
simul.res[[19]]
exp(coef(simul.res[[19]]))
exp(confint(simul.res[[19]]))

print('20th imputed datset:')
simul.res[[20]]
exp(coef(simul.res[[20]]))
exp(confint(simul.res[[20]]))

sink()

