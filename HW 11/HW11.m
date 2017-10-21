
clear all
gamma = 2;
beta = 0.96;
time_T = 60;
time_W = 47;

sigmaZeta = sqrt(0.0106);
sigmaEpsilon = sqrt(0.0738);
lambda = 0.68;
A = 1;
alpha = 1/3;
delta = 0.10;

z = [-2*sigmaZeta -sigmaZeta 0 sigmaZeta 2*sigmaZeta];
e = [-sigmaEpsilon sigmaEpsilon];


a = -10:0.1:25;

numA = length(a);
numZ = length(z);
numE = length(e);
numN = 10000;

w_up = 1.25;
w_low = 0.75; %1
k_diff = 1;
n_diff = 1;

while abs(n_diff)>0.001
    w = (w_up+w_low)/2;
    upperr = 0.1;
    lowerr = 0; %0.05
    while abs(k_diff)>0.05
        r = (upperr+lowerr)/2;
        q = 1/(1+r);
        ai = zeros(numA,numZ,numE,time_T);
        vtemp = zeros(numA,numZ,numE,time_T);
        
        for i = 1:numA         %looping through possible wealths
            for j = 1:numZ     %looping through permanent shocks
                for k = 1:numE %looping through transitory shocks
                    c = a(i);
                    if c<0
                        c = 10^-8;
                    end
                    vtemp(i,j,k,time_T) = c^(1-gamma)/(1-gamma);
                end
            end
        end
        
        for t = time_T-1:-1:1            %looping through the different age
            for i = 1:numA             %looping through possible wealths
                for j = 1:numZ
                    for k = 1:numE
                        values = zeros(numA,1);
                        for m = 1:numA
                            if t <= time_W
                                income = exp(-2.17 + 0.168*(t+18) - 0.032*(t+18)^2/10 + 0.002*(t+18)^3/100 + z(j) + e(k));
                                c = a(i)+w*income-q*a(m);
                                if c<0
                                    c = 10^-8;
                                end
                                util = c^(1-gamma)/(1-gamma);
                                if j==1
                                    values(m)=util +0.25*beta*(vtemp(m,gating(j,numZ),1,t+1)...
                                        +vtemp(m,gating(j,numZ),2,t+1)...
                                        +vtemp(m,gating(j+1,numZ),1,t+1)...
                                        +vtemp(m,gating(j+1,numZ),2,t+1));
                                elseif j==5
                                    values(m)=util...
                                        +0.25*beta*(vtemp(m,gating(j-1,numZ),1,t+1)...
                                        +vtemp(m,gating(j-1,numZ),2,t+1)...
                                        +vtemp(m,gating(j,numZ),1,t+1)...
                                        +vtemp(m,gating(j,numZ),2,t+1));
                                else
                                    values(m)=util...
                                        +0.25*beta*(vtemp(m,gating(j-1,numZ),1,t+1)...
                                        +vtemp(m,gating(j-1,numZ),2,t+1)...
                                        +vtemp(m,gating(j+1,numZ),1,t+1)...
                                        +vtemp(m,gating(j+1,numZ),2,t+1));
                                end
                            else
                                income = 0;
                                c = a(i)+w*income-q*a(m);
                                if c<0
                                    c = 10^-8;
                                end
                                util = c^(1-gamma)/(1-gamma);
                                values(m)=util+beta*(vtemp(m,j,k,t+1));
                            end
                        end
                        [vtemp(i,j,k,t), ai(i,j,k,t)] = max(values);
                    end
                end
            end
        end
        
        %matrices for all N individuals accros T time
        zAll = zeros(numN,time_T);
        eAll = zeros(numN,time_T);
        aAll = zeros(numN,time_T);
        cAll = zeros(numN,time_T);
        hAll = zeros(numN,time_T);
        
        for n = 1:numN
            %creating simulation for one person
            zPerson = zeros(1,time_T);
            ePerson = zeros(1,time_T);
            aPerson = zeros(1,time_T);
            cPerson = zeros(1,time_T);
            hPerson = zeros(1,time_T);
            
            for t = 2:time_T
                if t <= time_W
                    prevVal = find(z==zPerson(t-1));
                    if prevVal == 1
                        zPerson(t) = z(randi(2)); % 50/50 staying or going to z(2)
                    elseif prevVal == 5
                        zPerson(t) = z(5-randi(2)+1); % 50/50 staying or going to z(4)
                    else
                        random = randi(2); %50/50 for going one up or one down index
                        if random == 1
                            zPerson(t) = z(prevVal+1);
                        else
                            zPerson(t) = z(prevVal-1);
                        end
                    end
                else
                    zPerson(t) = zPerson(time_W);
                end
            end
            
            zAll(n,:) = zPerson;
            
            for t = 1:time_T
                ePerson(t) = e(randi(2));
                tempZ = find(z==zPerson(t));
                tempE = find(e==ePerson(t));
                if t <= time_W
                    income = exp(-2.17 + 0.168*(t+18) - 0.032*(t+18)^2/10 + 0.002*(t+18)^3/100 + z(tempZ) + e(tempE));
                else
                    income = 0;
                end
                hPerson(t) = income;
            end
            
            eAll(n,:) = ePerson;
            hAll(n,:) = hPerson;
            
            for t = 2:time_T
                aIndex = find(a== aPerson(t-1));
                zIndex = find(z== zPerson(t-1));
                eIndex = find(e== ePerson(t-1));
                aPerson(t) = a(ai(aIndex, zIndex, eIndex, t-1));
            end
            aAll(n,:) = aPerson;
            
            for t = 1:time_T-1
                cPerson(t) = w*hPerson(t)+aPerson(t)-q*aPerson(t+1);
            end
            cPerson(time_T) = w*hPerson(time_T)+aPerson(time_T);
            cAll(n,:) = cPerson;
            
        end %end of one person simulation looped for N people
        bigK = sum(sum(aAll));
        if bigK<0
            bigK = 0.00000001;
        end
        bigN = sum(sum(hAll));
        if bigN<0
            bigN = 0.00000001;
        end
        mpk = alpha*A*bigK^(alpha-1)*bigN^(1-alpha);
        k_diff = mpk-(r+delta);
        if k_diff>0
            lowerr = r;
        else
            upperr = r;
        end
    end
    bigK = sum(sum(aAll));
    if bigK<0
        bigK = 0.00000001;
    end
    bigN = sum(sum(hAll));
    if bigN<0
        bigN = 0.00000001;
    end
    mpn = (1-alpha)*A*bigK^alpha*bigN^(-alpha);
    n_diff = mpn-w;
    if n_diff>0
        w_low = w;
    else
        w_up = w;
    end
end




wealthDist = zeros(numN*time_T,1);
for i=1:numN
    for j = 1:time_T
        wealthDist((i-1)*time_T+j) = aAll(i,j);
    end
end

wealthDist = sort(wealthDist);
cumulativeWealth = cumsum(wealthDist);


%sanity check
plot(1:numN*time_T, cumulativeWealth)

p50 = round(numN*time_T*0.5);
p90 = round(numN*time_T*0.9);
p95 = round(numN*time_T*0.95);
p99 = round(numN*time_T*0.99);

tableResults = zeros(5,1);
tableResults(1) = cumulativeWealth(p50-1);
tableResults(2) = cumulativeWealth(p90-1)-tableResults(1);
tableResults(3) = cumulativeWealth(p95-1)-tableResults(2);
tableResults(4) = cumulativeWealth(p99-1)-tableResults(3);
tableResults(5) = cumulativeWealth(numN*time_T)-tableResults(4);

tableResults = tableResults./max(cumulativeWealth);


for i = 1:numN*time_T
    graphing(i) = i;
end

lorenz = trapz(graphing, cumulativeWealth);
height = cumulativeWealth(numN*time_T) - cumulativeWealth(1);
area = height*numN*time_T/2;
gini_wealth = (area-lorenz)/area


totalIncome = sum(sum(hAll))*w;
incDist = zeros(numN*time_T,1);
for i=1:numN
    for j = 1:time_T
        incDist((i-1)*time_T+j) = hAll(i,j)*w;
    end
end

incDist = sort(incDist);
cumulativeIncome = cumsum(incDist);


lorenzInc = trapz(graphing, cumulativeIncome);
heightInc = max(cumulativeIncome) - cumulativeIncome(1);
areaInc = heightInc*numN*time_T/2;
gini_income = (areaInc-lorenzInc)/areaInc



