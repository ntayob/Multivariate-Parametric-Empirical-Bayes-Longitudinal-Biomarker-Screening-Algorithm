function f = probability_cases2(xin,theta_hat_n_star,V_n_star)
f = 1-mvncdf(xin,theta_hat_n_star,V_n_star);
f=-f;