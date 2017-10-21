function wage = y(t,z,eta)
global lambda;
global t_rt;
global z_final;
if t<47
    ft = -2.17+0.16*(t+18)-0.032/10*(t+18)^2+0.002/100*(t+18)^3;
    wage=exp(ft+z+eta);
else
    ft_w=-2.17+0.16*(t_rt+18)-0.032/10*(t_rt+18)^2+0.002/100*(t_rt+18)^3;
    wage=lambda*exp(ft_w+z_final);
end

