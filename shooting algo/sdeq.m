function k_t2 = sdeq(k_t,k_t1,t_t,t_t1,alpha, beta,delta,gama,n)

k_t2 = (k_t1*(1-delta)+k_t1^alpha-t_t1-(k_t*(1-delta)+k_t^alpha-t_t-k_t1*(1+gama)*(1+n))*(alpha*k_t1^(alpha-1)+1-delta)*beta/(1+gama))/(1+gama)/(1+n);



% k_t2 =(k_t1*(1-delta)+k_t1^alpha-beta/(1+gama)*(alpha*k_t1^(alpha-1)+1-delta)*(k_t^alpha+k_t*(1-delta)-k_t1*(1+gama)*(1+n)))/((1+gama)*(1+n));
end
% (k2*(1-delta)+k2^alpha-beta/(1+gamma)*(alpha*k1^(alpha-1)+1-delta)*(k1*(1-delta)+k1^alpha-k2*(1+gamma)*(1+n)))/((1+gamma)*(1+n));