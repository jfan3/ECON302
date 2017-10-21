t_end = 2;
gama = 2; %risk aversion
global beta 
beta= 0.96; 
q=beta^(-1);
y_s1 = 1;
y_s2 = 0;
m=21;
n=2;
a= -1:0.1:1;
v1=zeros(m,n);
a2=zeros(m,n);

for i=1:m
    for j=1:n
        val=zeros (m,1);
        
        for k=1:m
            
            val(k)=(gating(a(i)+(n==1)-q*a(k)))^(1-gama)/(1-gama)+ 0.5*beta*...
                ((gating(a(k)+1))^(1-gama)/(1-gama)+(gating(a(k)))^(1-gama)/(1-gama));
        end
        [v1(i,j),a2(i,j)]=max(val);
    end 
end
figure (1)
plot(a,v1(:,1),a,a2(:,1));%,a,v1(:,2),a,a2(:,1),a,a2(:,2)
