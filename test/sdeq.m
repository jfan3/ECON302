function k_t2 = simaa(k_t,k_t1,phi)
global alpha beta delta
r_t=alpha*k_t^(alpha-1);
r_t1=alpha*k_t1^(alpha-1);
y_t= k_t^alpha; g_t = 0.2*y_t;
y_t1 = k_t1^alpha; g_t1=0.2*y_t1;
pi_t= (y_t-alpha*y_t)-g_t;
pi_t1=(1-alpha)*y_t1-g_t1;

c_t=(1-delta)*k_t+r_t*k_t+pi_t-k_t1-phi*(k_t1-k_t)^2;
aa=beta*c_t/(1+2*phi*(k_t1-k_t));
bb=(1-delta)*k_t1+r_t1*k_t1+pi_t1;

if (phi==0)
    k_t2=(bb-aa*(r_t1+(1-delta)));
else
    a=phi;
    b=(-2)*phi*k_t1+1+aa*2*phi;
    c=aa*(r_t1-2*phi*k_t1+1-delta)-bb+phi*k_t1^2;
    k_t2=(-b+(b^2-4*a*c)^0.5)/(2*a);


% syms x 
% S=solve ((1-delta)*k_t1+r_t1*k_t1+pi_t1-x-phi*(x-k_t1)^2-...
% c_t*beta*((r_t1+2*phi*(x-k_t1)+1-delta)/(1+2*phi*(k_t1-k_t))));
% k_t2=double(S);
end