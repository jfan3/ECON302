
function k_series = test(phi)
global T beta delta alpha;
T =100;
beta=0.96;
delta=0.06;
alpha=0.36;
% phi = 0;
% phis=[0 0.1 0.2];

k_star= (alpha/(beta^(-1)-(1-delta)))^(1/(1-alpha));
k0=0.5*k_star;

% for n=1:length(phis)
    k_series= zeros(1,T);
    c_series= zeros(1,T-1);
    k_series(1)=k0;
    lower = k0;
    upper = k_star;
    % ki = (lower+upper)/2;
    i=1;


    difference = k_star-k0;

    while abs(difference)>0.001

        k_series(3:T)= zeros(1,T-2);
        ki = (lower+upper)/2;
        k_series(2)=ki;

        for a=3:100
        k_series(a)=sdeq(k_series(a-2),k_series(a-1),phi);
            if k_series(a)<0 ||k_series(a)>k_star
                break
            end
        end


        if k_series(a)>k_star
            upper=ki;

        else 
            lower=ki;
        end

        i=i+1;
        difference = k_series(T)-k_star;
    end




    for b = 1:99
        c_series(b)= (1-delta)*k_series(b)+0.8*k_series(b)^alpha-k_series(b+1)-phi...
            *(k_series(b+1)-k_series(b))^2;


    end
