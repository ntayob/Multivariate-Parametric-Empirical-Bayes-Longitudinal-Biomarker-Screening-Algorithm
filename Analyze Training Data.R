########################################################################
#Load Packages
########################################################################
library(nlme)
library(reshape2)

T_int=1 #Duration of pre-clinical period of interest (Tau in Manuscript)

########################################################################
#Load Training Data (K=3 biomarkers)
########################################################################
train_data=read.csv("data_training.csv",stringsAsFactors = F)
names(train_data)
#ID: patient identifier
#Cancer: disease status of patient (0 if never diagnosed with cancer, 1 if diagnosed with cancer during follow-up)
#d: time from first screening visit to the clinical detection of disease if the patient is diseased and end of study time if the patient is disease-free
#t: time since first screening visit for each patient visit
#J: total number of visits for the patient
#obs_number: index of screening visits
#Y1: levels of biomarker 1
#Y2: levels of biomarker 2
#Y3: levels of biomarker 3

########################################################################
#Get parameter estimates for multivariate models
########################################################################
#Change data shape from wide to long
sdata=melt(train_data[,c("ID","Cancer","d","t","Y1","Y2","Y3")],measure.vars=c("Y1","Y2","Y3"))

#Estimate parameters in patients that are disease-free during the study
control_data=subset(sdata,sdata$Cancer==0)
model_controls <- lme(value ~ variable - 1, random = ~variable - 1 | ID, weights = varIdent(form = ~1 | variable), correlation = corSymm(form = ~1 | ID/t), control=lmeControl(msMaxIter = 200, msMaxEval = 500, opt="optim"), data = control_data)
#summary(model_controls)

#Mean vector of multivariate normal distribution
mu_theta_hat=as.numeric(model_controls$coefficients$fixed)
write.csv(mu_theta_hat,"mu.csv",row.names=FALSE)

#Between-subject covariance matrix
temp=suppressWarnings(as.numeric(VarCorr(model_controls)))
sigma_theta_hat=diag(temp[1:3])
sigma_theta_hat[1,2]=sigma_theta_hat[2,1]=temp[10]*temp[5]*temp[6]
sigma_theta_hat[1,3]=sigma_theta_hat[3,1]=temp[11]*temp[5]*temp[7]
sigma_theta_hat[2,3]=sigma_theta_hat[3,2]=temp[15]*temp[6]*temp[7]
write.csv(sigma_theta_hat,"sigma_theta.csv",row.names=FALSE)

#Within-subject covariance matrix
temp2=coef(model_controls$modelStruct$corStruct,unconstrained=FALSE)
temp3=coef(model_controls$modelStruct$varStruct,unconstrained=FALSE)

sigma_hat=diag(c(temp[8]*temp[8],temp[8]*temp[8]*temp3[1]*temp3[1],temp[8]*temp[8]*temp3[2]*temp3[2]))
sigma_hat[1,2]=sigma_hat[2,1]=temp[8]*temp[8]*temp3[1]*temp2[1]
sigma_hat[1,3]=sigma_hat[3,1]=temp[8]*temp[8]*temp3[2]*temp2[2]
sigma_hat[2,3]=sigma_hat[3,2]=temp[8]*temp3[1]*temp[8]*temp3[2]*temp2[3]
write.csv(sigma_hat,"sigma.csv",row.names=FALSE)

#Estimate parameters in patients that develop disease during the study
#Exclude cases with flat trajectory for all biomarkers to increase signal
temp_case_data=subset(train_data,Cancer==1)
case_ID=unique(temp_case_data$ID)
n_1=length(case_ID)
increasing_1=rep(NA,n_1)
increasing_2=rep(NA,n_1)
increasing_3=rep(NA,n_1)
for(j in 1:n_1)
{
  if(unique(temp_case_data$J[temp_case_data$ID==case_ID[j]])>1)
  {
    model_1=lm(Y1~t,data=subset(temp_case_data,ID==case_ID[j]))
    increasing_1[j]=as.numeric(coef(model_1)[2]>0)
    model_2=lm(Y2~t,data=subset(temp_case_data,ID==case_ID[j]))
    increasing_2[j]=as.numeric(coef(model_2)[2]>0)
    model_3=lm(Y3~t,data=subset(temp_case_data,ID==case_ID[j]))
    increasing_3[j]=as.numeric(coef(model_3)[2]>0)
  }
}
remove_ID=case_ID[increasing_1==0 & increasing_2==0 & increasing_3==0 & is.na(increasing_1)==F & is.na(increasing_2)==F & is.na(increasing_3)==F]
keep_ID=case_ID[!(case_ID %in% remove_ID)]

#Only keep visits within pre-clinical period of interest
case_data=subset(sdata,sdata$Cancer==1 & (sdata$d-sdata$t)<=T_int)
case_data=subset(case_data,ID %in% keep_ID)

model_cases <- lme(value ~ variable - 1, random = ~variable - 1 | ID, weights = varIdent(form = ~1 | variable), correlation = corSymm(form = ~1 | ID/t), control=lmeControl(msMaxIter = 200, msMaxEval = 500, opt="optim"), data = case_data)
#summary(model_cases)

#Mean vector of multivariate normal distribution
mu_theta_star_hat=as.numeric(model_cases$coefficients$fixed)
write.csv(mu_theta_star_hat,"mu_star.csv",row.names=FALSE)

#Between-subject covariance matrix
temp=suppressWarnings(as.numeric(VarCorr(model_cases)))
sigma_theta_star_hat=diag(temp[1:3])
sigma_theta_star_hat[1,2]=sigma_theta_star_hat[2,1]=temp[10]*temp[5]*temp[6]
sigma_theta_star_hat[1,3]=sigma_theta_star_hat[3,1]=temp[11]*temp[5]*temp[7]
sigma_theta_star_hat[2,3]=sigma_theta_star_hat[3,2]=temp[15]*temp[6]*temp[7]
write.csv(sigma_theta_star_hat,"sigma_theta_star.csv",row.names=FALSE)

#Within-subject covariance matrix
temp2=coef(model_cases$modelStruct$corStruct,unconstrained=FALSE)
temp3=coef(model_cases$modelStruct$varStruct,unconstrained=FALSE)

sigma_star_hat=diag(c(temp[8]*temp[8],temp[8]*temp[8]*temp3[1]*temp3[1],temp[8]*temp[8]*temp3[2]*temp3[2]))
sigma_star_hat[1,2]=sigma_star_hat[2,1]=temp[8]*temp[8]*temp3[1]*temp2[1]
sigma_star_hat[1,3]=sigma_star_hat[3,1]=temp[8]*temp[8]*temp3[2]*temp2[2]
sigma_star_hat[2,3]=sigma_star_hat[3,2]=temp[8]*temp3[1]*temp[8]*temp3[2]*temp2[3]
write.csv(sigma_star_hat,"sigma_star.csv",row.names=FALSE)


