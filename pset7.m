ro= [0.95 0.98];
sigma_e = [0.007 0.014];
z0 = 0;
delta = 0.0224;
error_1 = sigma_e(1)*randn(1,100);
error_3 = sigma_e(2)*randn(1,100);
z_series_1=zeros(1,100);
z_series_2=zeros(1,100);
z_series_3=zeros(1,100);
stoch_val=[-delta;delta];
transition= [0.975 0.025 ;0.975 0.025];
states=zeros(2,100);
states(:,1)=stoch_val;
for a=2:length(z_series_1)
    z_series_1(a)= ro(1)*z_series_1(a-1)+error_1(a);
    z_series_2(a)= ro(2)*z_series_2(a-1)+error_1(a);
    z_series_3(a)= ro(1)*z_series_3(a-1)+error_3(a);
    states(:,a)=transition*states(:,a-1);
end
figure (1)
plot(1:100,z_series_1,1:100,z_series_2,1:100,z_series_3,1:100,states(1,:),1:100,states(2,:));
legend('ro=0.95,\sigma=0.007','ro=0.98,\sigma=0.007','ro=0.95,\sigma=0.014',...
    '2-state markov (low)','2-state markov (high)')
ylabel('total prductivity factor')
xlabel('time')

