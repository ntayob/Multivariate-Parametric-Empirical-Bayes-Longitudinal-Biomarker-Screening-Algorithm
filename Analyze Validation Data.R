########################################################################
#Load Validation Data after running mPEB in Matlab
########################################################################
#Output from Matlab
data_validation_mPEB=read.table("data_validation_output.txt")
names(data_validation_mPEB)=c("ID","t","f_075","f_080","f_085","f_090","f_095","f_100","f_105","f_110","f_115","f_120","f_125") 

#Original validation dataset
data_validation=read.csv("data_validation.csv",na.strings = F)

#Reduce number of significant digits to facilitate merging files
data_validation_mPEB$t=round(data_validation_mPEB$t,5)
data_validation$t=round(data_validation$t,5)

temp_data=subset(data_validation,select=c("ID","D","d","t"))
validation_analysis_data=merge(temp_data,data_validation_mPEB,by=c("ID","t"))

remove(data_validation,data_validation_mPEB,temp_data)

#########################################################################################
#Extract ROC(0.1) for multiple time-frames of interest prior to clinical diagnosis

#fpr: false positive rate (screening-level)
#tpr: true positive rate (patient-level)
#########################################################################################
specificity_fixed=0.9 #1-0.1=1-fpr=specificity

control_data=subset(validation_analysis_data,D==0)
#observed fpr
fpr=apply(control_data[,c("f_075","f_080","f_085","f_090","f_095","f_100","f_105","f_110","f_115","f_120","f_125")],2,mean)
index=seq(from=1,to=length(fpr)) 
selected_index=max(index[(abs(fpr-(1-specificity_fixed))==min(abs(fpr-(1-specificity_fixed))))]) #select index with observed specificity=1-fpr that matches target specificity

#Positive screen at any time counts towards patient-level sensitivity calculation
case_data=subset(validation_analysis_data,D==1)
subject_ID=unique(case_data$ID)
number_positive=array(NA,c(length(subject_ID),length(fpr)))
for(j in 1:length(subject_ID))
{
  subject_data=subset(case_data,ID==subject_ID[j])
  number_positive[j,]=apply(subject_data[,c("f_075","f_080","f_085","f_090","f_095","f_100","f_105","f_110","f_115","f_120","f_125")],2,sum)
}  
tpr=apply(number_positive>0,2,mean)
roc_t=tpr[selected_index]
est_t=fpr[selected_index]

#Print results
round(roc_t*100,2)
#81.63
round(est_t*100,2)
#10.18 

#Only positive screens within 1-year prior to clinical diagnosis count towards patient-level sensitivity calculation
case_data=subset(validation_analysis_data,D==1 & (d-t)<=1)
subject_ID=unique(case_data$ID) #only keep patients with at least one screening visit during this period
number_positive=array(NA,c(length(subject_ID),length(fpr)))
for(j in 1:length(subject_ID))
{
  subject_data=subset(case_data,ID==subject_ID[j])
  number_positive[j,]=apply(subject_data[,c("f_075","f_080","f_085","f_090","f_095","f_100","f_105","f_110","f_115","f_120","f_125")],2,sum)
}  
tpr=apply(number_positive>0,2,mean)
roc_t_1year=tpr[selected_index]
est_t_1year=fpr[selected_index]

#Print results
round(roc_t_1year*100,2)
#67.35
round(est_t_1year*100,2)
#10.18 

#Only positive screens between 1 to 2 years prior to clinical diagnosis count towards patient-level sensitivity calculation
case_data=subset(validation_analysis_data,D==1 & (d-t)>1 & (d-t)<=2)
subject_ID=unique(case_data$ID) #only keep patients with at least one screening visit during this period
number_positive=array(NA,c(length(subject_ID),length(fpr)))
for(j in 1:length(subject_ID))
{
  subject_data=subset(case_data,ID==subject_ID[j])
  number_positive[j,]=apply(subject_data[,c("f_075","f_080","f_085","f_090","f_095","f_100","f_105","f_110","f_115","f_120","f_125")],2,sum)
}  
tpr=apply(number_positive>0,2,mean)
roc_t_1_2year=tpr[selected_index]
est_t_1_2year=fpr[selected_index]

#Print results
round(roc_t_1_2year*100,2)
#28.95
round(est_t_1_2year*100,2)
#10.18 

#Only positive screens more than 2 years prior to clinical diagnosis count towards patient-level sensitivity calculation
case_data=subset(validation_analysis_data,D==1 & (d-t)>2)
subject_ID=unique(case_data$ID) #only keep patients with at least one screening visit during this period
number_positive=array(NA,c(length(subject_ID),length(fpr)))
for(j in 1:length(subject_ID))
{
  subject_data=subset(case_data,ID==subject_ID[j])
  number_positive[j,]=apply(subject_data[,c("f_075","f_080","f_085","f_090","f_095","f_100","f_105","f_110","f_115","f_120","f_125")],2,sum)
}  
tpr=apply(number_positive>0,2,mean)
roc_t_gt2year=tpr[selected_index]
est_t_gt2year=fpr[selected_index]

#Print results
round(roc_t_gt2year*100,2)
#41.94
round(est_t_gt2year*100,2)
#10.18 


