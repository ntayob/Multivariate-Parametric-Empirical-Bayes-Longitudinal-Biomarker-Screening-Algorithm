function [c ceq] = probability_controls2(xin,theta_hat_n,V_n,f_0)
c=[];
ceq = 1-mvncdf(xin,theta_hat_n,V_n)-f_0;

