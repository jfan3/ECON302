global lambda;
global t_rt; 
global t_d;
t_rt = 47;
t_d = 60;
gama = 2; %risk aversion
beta= 0.96; 
q=beta;
sigma_zeta=0.0106^0.5;
sigma_eta=0.0738^0.5;
lambda= 0.68;
z=[-2*sigma_zeta -sigma_zeta 0 sigma_zeta sigma_zeta*2];
eta=[-sigma_eta sigma_eta];
% a_t = zeros(1,length(z)*length(eta));
a=-2:0.1:15;
m=length(a);
n=length(z);
l=length(eta);
N= 10000;
a_star=zeros(m,n,l,t_d);
v_star=zeros(m,n,l,t_d);
global z_final;
z_final = z(randi([1,5], 1,1));
values=zeros(m,n,l,t_d);
count=0;
count1=0;

%final period condition

for k3=1:m
    for i3=1:n
        for j3=1:l
            values(k3,i3,j3,t_d)=util(gating(a(k3)+y(t_d,z_final,eta(j3))));
        end
    end
end

%policy rule
for t = t_d-1:-1:1
%     consumption= zeros(m,n,l);
%     v2=zeros(m,n,l);
       
    for k2=1:m %a
        for i=1:n %z
            for j=1:l %eta
                values_temp=zeros(m,1);
                for k1=1:m %a'
                    if t<=t_rt
                        values_temp(k1)=util(gating(a(k2)+y(t,z(i),eta(j))-q*a(k1)))+...
                            0.25*beta*(values(k1,sub(i-1,n),sub(j-1,l),t+1)+...
                            values(k1,sub(i-1,n),sub(j+1,l),t+1)+...
                            values(k1,sub(i+1,n),sub(j-1,l),t+1)+...
                            values(k1,sub(i+1,n),sub(j+1,l),t+1));
    %                        +0.25*beta*(util(gating(a(k1)+y(t,z(sub(i-1,n)),eta(sub(j-1,l)))))+...
    %                        util(gating(a(k1)+y(t,z(sub(i-1,n)),eta(sub(j+1,l)))))+...
    %                        util(gating(a(k1)+y(t,z(sub(i+1,n)),eta(sub(j-1,l)))))+...
    %                        util(gating(a(k1)+y(t,z(sub(i+1,n)),eta(sub(j+1,l))))));
    %                         +0.25*beta*(v2(k1,sub(i-1,n),sub(j-1,l)))...
    %                         +v2(k1,sub(i-1,n),sub(j+1,l))...
    %                         +v2(k1,sub(i+1,n),sub(j-1,l))...
    %                         +v2(k1,sub(i+1,n),sub(j+1,l)));
                    else
                        values_temp(k1)=util(gating (a(k2)+y(t,z(i),eta(j)-q*a(k1))))+...
                            beta*values(k1,sub(i,n),sub(j,l),t+1);
                    end
                    
                    count1=count1+1;
                end
                [values(k2,i,j,t), a_star(k2,i,j,t)]=max(values_temp);
            end
        end 
    end
end
% data = zeros(1,t_d);
% data (1,:)=a_star(1,1,1,:);
% figure (1);
% plot(1:t_d,data)

%simulation
%matrices for all N individuals accross T time
zAll = zeros(N,t_d);  
eAll = zeros(N,t_d);
yAll = zeros(N,t_d);
aAll = zeros(N,t_d);
cAll = zeros(N,t_d);

f = @(t) -2.17 + 0.168*(t+18) - 0.032*(t+18)^2/10 + 0.002*(t+18)^3/100;

for n = 1:N

   
    zv = zeros(1,t_d);
    zv(1) = 3;
    for i = 2:t_d
        if zv(i-1) == 1
            zv(i) = 1+randi(0:1); % 50% stay at i = 1; 50% move to i = 2
        elseif zv(i-1) == 5
            zv(i) = 5-randi(0:1); % 50% stay at i = 5; 50% move to i = 4
        else
            if randi(0:1) == 0        
                zv(i) = zv(i-1)-1; % 50% i = i-1
            else
                zv(i) = zv(i-1)+1; % 50% i = i+1
            end
        end
    end
    zAll(n,:) = zv;

    % create a vector to hold all epsilon's
    ev = zeros(1,t_d);
    ev(1) = randi(1:2);
    for i = 1:t_d
        ev(i) = randi(1:2);
    end
    eAll(n,:) = ev;

    % create a vector to hold all a's
    av = zeros(1,t_d);
    [~,av(1)] = min(abs(a)); 
    for i = 1:t_d-1
        av(i+1) = a_star(av(i),zv(i),ev(i),i);
    end
    aAll(n,:) = av;
    
    % calculate y based on z's and epsilon's
    yv = zeros(1,t_d);
    for i = 1:t_d
        if i > t_rt
            yv(i) = lambda*exp(f(t_rt)+z(zv(t_rt)));
        else
    		yv(i) = exp(f(i)+z(zv(i))+eta(ev(i))); 
        end
    end
    yAll(n,:) = yv;
    
    % create a vector to hold all c's
    cv = zeros(1,t_d);
    for i = 1:t_d-1
        cv(i) = a(av(i)) + yv(i) - q*a(av(i+1));
    end
    cv(T) = a(av(T)) + yv(T);
    cAll(n,:) =cv;
    


end

% graph average income, consumption, and wealth
ave_y = mean(yAll);
ave_c = mean(cAll);
ave_a = mean(a(aAll));
figure;
plot([1:60],ave_y,'b-', [1:60],ave_c, 'r--', [1:60],ave_a, 'k-.');
legend('y', 'c', 'a');
title('Average Income, Consumption, and Wealth');
xlabel('time');
ylabel('Income/Consumption/Wealth');

% yAvg = zeros(1,t_d);
% cAvg = zeros(1,t_d);
% aAvg = zeros(1,t_d);
% 
% for t = 1:t_d
%     yAvg(t) = mean(yAll(:,t));
%     cAvg(t) = mean(cAll(:,t));
%     aAvg(t) = mean(aAll(:,t));
% end
% 
% figure(1)
% plot(1:t_d, yAvg, 1:t_d, cAvg,'--', 1:t_d, aAvg, 'o-');
% legend('Income', 'Consumption', 'Wealth')
% xlabel('Time')
% % 
% % yOrig = yAvg;
% % cOrig = cAvg;
% % aOrig = aAvg;
% % 
% % yLamba = yAvg;
% % cLambda = cAvg;
% % aLambda = aAvg;
% 
% % yZ = yAvg;
% % cZ = cAvg;
% % aZ = aAvg;
% 
% yE = yAvg;
% cE = cAvg;
% aE = aAvg;
% 
% figure(2)
% plot(1:t_d, yAvg, 1:t_d, yLamba,'--', 1:t_d, yZ, 'o-', 1:t_d, yE);
% legend('Original', 'Lambda', 'Z', 'E')
% xlabel('Time')
% ylabel('Income')
% 
% figure(3)
% plot(1:t_d, cAvg, 1:t_d, cLambda,'--', 1:t_d, cZ, 'o-', 1:t_d, cE);
% legend('Original', 'Lambda', 'Z', 'E')
% xlabel('Time')
% ylabel('Consumption')
% 
% figure(4)
% plot(1:t_d, aAvg, 1:t_d, aLambda,'--', 1:t_d, aZ, 'o-', 1:t_d, aE);
% legend('Original', 'Lambda', 'Z', 'E')
% xlabel('Time')
% ylabel('Wealth')
% 
% load handel;
% player = audioplayer(y, Fs);
% play(player);
