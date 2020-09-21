clear all; %clear memory
close all; %closes figures
clc; %clear console

options.Algorithm = 'sqp';
options.MaxFunEvals = 1000;


%Read in validation data
Y=csvread('data_validation.csv',1,6);
Y1=Y(:,1);
Y2=Y(:,2);
Y3=Y(:,3);
[N, b]=size(Y); %N=Number of screenings and b=number of biomarkers
%D=csvread('data_validation.csv',1,1,[1,1,N,1]);
%d=csvread('data_validation.csv',1,2,[1,2,N,2]);
obs_number=csvread('data_validation.csv',1,5,[1,5,N,5]);
ID=csvread('data_validation.csv',1,0,[1,0,N,0]);
Time=csvread('data_validation.csv',1,3,[1,3,N,3]);
subject_ID=ID(obs_number==1); %unique subject ids in dataset
n=length(subject_ID);

%disease-free model parameters
mu=csvread('mu.csv',1,0);
sigma=csvread('sigma.csv',1,0);
sigma_theta=csvread('sigma_theta.csv',1,0);

%diseased model parameters
mu_star=csvread('mu_star.csv',1,0);
sigma_star=csvread('sigma_star.csv',1,0);
sigma_theta_star=csvread('sigma_theta_star.csv',1,0);

%lower and upperbounds of biomarker values observed
x_lb=[quantile(Y1,0);quantile(Y2,0);quantile(Y3,0)];
x_ub=[quantile(Y1,1);quantile(Y2,1);quantile(Y3,1)];

%Explore multiple population-level false positive rates around target of 10
%percent
f_0=0.075:0.005:0.125;
number_positive_tests=zeros(N,length(f_0)); %counter of number of positive tests at each visit
for i=1:length(f_0)
    l=0;
    for j=1:1:n
        %Get biomarker levels for jth patient
        temp1=Y1.*(ID==j);
        temp1(temp1==0) = [];
        
        temp2=Y2.*(ID==j);
        temp2(temp2==0) = [];
        
        temp3=Y3.*(ID==j);
        temp3(temp3==0) = [];
        for k=1:1:sum(ID==j) %iterate through each screening visit
            if k==1
                Y_in_bar=[0; 0; 0]; %sample average of prior screenings is a 0 vector if first screening visit
            else
                temp1_1=temp1(linspace(1,k-1,k-1));
                temp2_1=temp2(linspace(1,k-1,k-1));
                temp3_1=temp3(linspace(1,k-1,k-1));
                Y_in_bar=[mean(temp1_1); mean(temp2_1); mean(temp3_1)]; %sample average of prior screenings
            end
            
            %Conditional distribution in disease-free patients
            theta_hat_n=(inv(inv(sigma_theta) + (k-1)*inv(sigma)))*(inv(sigma_theta)*mu+(k-1)*inv(sigma)*Y_in_bar);
            V_n=inv(inv(sigma_theta) + (k-1)*inv(sigma)) + sigma;
            
            %Conditional distribution in diseased patients
            theta_hat_n_star=(inv(inv(sigma_theta_star) + (k-1)*inv(sigma_star)))*(inv(sigma_theta_star)*(mu_star) + (k-1)*inv(sigma_star)*Y_in_bar);
            V_n_star=inv(inv(sigma_theta_star) + (k-1)*inv(sigma_star)) + sigma_star;
            
            %consider multiple starting values
            x_matrix=zeros(19,3);
            fval_matrix=zeros(19,1);
            exitflag_matrix=zeros(19,1);
            for w=1:1:19
                p_set=0.05*w; %starting points defined by 5th, 10th, 15th ..., 95th quantiles of the biomarker levels
                x_0=[quantile(Y1,p_set);quantile(Y2,p_set);quantile(Y3,p_set)];
                %implement numerical optimization
                probability_cases2=@(xin) probability_cases(xin,theta_hat_n_star,V_n_star);
                probability_controls2=@(xin) probability_controls(xin,theta_hat_n,V_n,f_0(i));
                [x,fval,exitflag]=fmincon(probability_cases2,x_0,[],[],[],[],x_lb,x_ub,probability_controls2,options);
                %x: thresholds
                % -(sensitivity) at thresholds
                %exitflag: errors in optimization (only keep result if
                %exitflag_matrix==1 indicating local minimum found)
                x_matrix(w,:)=x;
                fval_matrix(w,1)=fval;
                exitflag_matrix(w,1)=exitflag;
            end
            
            %find global optimal solution
            fval_completed=fval_matrix.*(exitflag_matrix==1);
            x_completed=x_matrix.*([exitflag_matrix exitflag_matrix exitflag_matrix]==1);
            fval_min=min(fval_completed);
            if fval_min<0
                x_min=x_completed.*([fval_completed fval_completed fval_completed]==fval_min);
                x_final=transpose(max(unique(x_min,'rows')));
                Y_i_n_1=[temp1(k); temp2(k); temp3(k)];
                number_positive_tests(l+k,i)=sum(Y_i_n_1>x_final); %count the number of biomarkers that exceed their threshold at each screening visit
            end
        end
        l=l+sum(ID==j);
    end
end

temp_input=[ID,Time,number_positive_tests];
%output results for analysis in R
dlmwrite('data_validation_output.txt',temp_input,'delimiter','\t','precision',15);










