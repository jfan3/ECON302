
k_series_0=test(0);
k_series_1=test(0.1);
k_series_2=test(0.2);

figure (1)
plot(1:100,k_series_0,'b',1:100,k_series_1,'r',1:100,k_series_2,'k')
legend('\Phi=0','\Phi=0.1','\Phi=0.2')
xlabel('time')
ylabel('k_t')
title('k_t evolution')


dk_1=zeros(1,99);
dk_2=zeros(1,99);
dk_3=zeros(1,99);
for i=1:99
    dk_1(i)=k_series_0(i+1)-k_series_0(i);
    dk_2(i)=k_series_1(i+1)-k_series_1(i);
    dk_3(i)=k_series_2(i+1)-k_series_2(i);
end

figure (2)
plot(1:99,dk_1,'b',1:99,dk_2,'r',1:99,dk_3,'k')
legend('\Phi=0','\Phi=0.1','\Phi=0.2')
xlabel('time')
ylabel('dk')
title('\Deltak_t/dt')

%The intuition behind the evolution of k_t can be explained by the fact
%that when its is more costly to make adjustments to k_t (phi increases),
%They tend to make less investment changes (explains why when phi=0.2, dk_t 
% and k_t are the smallest at first). The agent will reach a point where k_t has somewhat 
% spbilized, then the agent is somewhat more indifferent to making changes
% to k_t, and all three lines converge to k*.

% ddk_1=zeros(1,98);
% ddk_2=zeros(1,98);
% ddk_3=zeros(1,98);
% for i=1:98
%     ddk_1(i)=dk_1(i+1)-dk_1(i);
%     ddk_2(i)=dk_2(i+1)-dk_2(i);
%     ddk_3(i)=dk_3(i+1)-dk_3(i);
% end
% figure (3)
% plot(1:98,ddk_1,'b',1:98,ddk_2,'r',1:98,ddk_3,'k')
% legend('\Phi=0','\Phi=0.1','\Phi=0.2')
% xlabel('time')
% ylabel('ddk')