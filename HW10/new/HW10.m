% It will take around 5-10 min to run.

global gamma 
global a z e W lambda q T beta

gamma = 2;
beta = 0.96;
q = beta;
T = 60;
W = 47;

sigmaZeta = sqrt(0.0106);
sigmaEpsilon = sqrt(0.0738);
lambda = 0.68;


z = [-2*sigmaZeta -sigmaZeta 0 sigmaZeta 2*sigmaZeta];
e = [-2*sigmaEpsilon 2*sigmaEpsilon];

%as given in class
a = -2:0.1:20;

A = length(a);
Z = length(z);
E = length(e);
N = 10000;

a_star = zeros(A,Z,E,T);
vtemp = zeros(A,Z,E,T);
zAll = zeros(N,T);  
eAll = zeros(N,T);
yAll = zeros(N,T);
aAll = zeros(N,T);
cAll = zeros(N,T);

for i = 1:A         %looping through possible wealths
    for j = 1:Z     %looping through permanent shocks
        for k = 1:E %looping through transitory shocks
            vtemp(i,j,k,T) = util(y(j,k,T)+a(i));
        end
    end
end

for t = T-1:-1:1            %looping through the different age
    for i = 1:A             %looping through possible wealths
        for j = 1:Z
            for k = 1:E
                values = zeros(A,1);    
                for m = 1:A
                    if t <= W    
                        values(m)=util(a(i)+y(j,k,t)-q*a(m))...
                            +0.25*beta*(vtemp(m,gating(j-1,Z),1,t+1)...
                            +vtemp(m,gating(j-1,Z),2,t+1)...
                            +vtemp(m,gating(j+1,Z),1,t+1)...
                            +vtemp(m,gating(j+1,Z),2,t+1));
                    else
                        values(m)=util(a(i)+y(j,k,t)-q*a(m))...
                            +beta*(vtemp(m,j,k,t+1));
                    end  
                end
                [vtemp(i,j,k,t), a_star(i,j,k,t)] = max(values);
            end
        end
    end
end




for n = 1:N
    
    z_p = zeros(1,T);
    e_p = zeros(1,T);
    y_p = zeros(1,T);
    a_p = zeros(1,T);
    c_p = zeros(1,T);
    
    for t = 2:T
        if t <= W
            prevVal = find(z==z_p(t-1));
            if prevVal == 1
                z_p(t) = z(randi(2)); % 50/50 staying or going to z(2)
            elseif prevVal == 5
                z_p(t) = z(5-randi(2)+1); % 50/50 staying or going to z(4)
            else
                random = randi(2); %50/50 for going one up or one down index
                if random == 1
                    z_p(t) = z(prevVal+1);
                else
                    z_p(t) = z(prevVal-1);
                end
            end
        else
            z_p(t) = z_p(W);
        end
    end
    
    zAll(n,:) = z_p;
    
    for t = 1:T
        e_p(t) = e(randi(2));
        y_p(t) = y(find(z==z_p(t)),find(e==e_p(t)),t);
    end
    
    eAll(n,:) = e_p;
    yAll(n,:) = y_p;
    
    for t = 2:T
        aIndex = find(a== a_p(t-1));
        zIndex = find(z== z_p(t-1));
        eIndex = find(e== e_p(t-1));
        a_p(t) = a(a_star(aIndex, zIndex, eIndex, t-1));
    end
    aAll(n,:) = a_p;
    
    for t = 1:T-1
        c_p(t) = y_p(t)+a_p(t)-q*a_p(t+1);
    end
    c_p(T) = y_p(T)+a_p(T);

    cAll(n,:) = c_p;
end

yAvg = zeros(1,T);
cAvg = zeros(1,T);
aAvg = zeros(1,T);

for t = 1:T
    yAvg(t) = mean(yAll(:,t));
    cAvg(t) = mean(cAll(:,t));
    aAvg(t) = mean(aAll(:,t));
end

figure(1)
plot(1:T, yAvg, 1:T, cAvg,'r--', 1:T, aAvg, 'k-');
legend('Income', 'Consumption', 'Wealth')
xlabel('Time')
