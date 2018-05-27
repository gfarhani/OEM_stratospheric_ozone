OHP_2017 = load('SignalOHP2017.mat');
alt = (OHP_2017.alt)'*1000;
zmsis = alt';  
data = load('OHP_26072017.mat');
real_332 = data.Sig_Matrix(:,6)*728143 ;
D = find(zmsis >=50000);
D1 = min(D);
F = find(zmsis <= 160000);
F1 = max(F);
bg_332 = mean(real_332(D1:F1));
real_332_corr = real_332 - bg_332;
% semilogx(real_332_corr, zmsis)
% figure;
counts = real_332_corr (85:250);
H = zmsis(85:250);
tt = counts'.*(H .^2);
% semilogx(tt, H); hold all;
nair = X.x(m+1:2*m-20);
zret = Q.zret(1:end-20);
semilogx(nair, zret)
scale = interp1(zret, nair, H);
semilogx(scale, H, 'r', tt, H, 'b')
constant = scale(50)/tt(50);
comp = ((constant)*tt);
% comp = tt*10^(12.85397);
figure;
semilogx(comp, H, 'b',scale, H, 'r')
avg = (comp+scale)/2;
plot(((comp-scale)./avg)*100, H/1000)

%%%%%%%%%%
 %%%% scale the OEM to the MSIS model at high altitudes
 figure; 
%  semilogx(X.x(m+1:2*m), Q.zret/1000, 'r');
 
 semilogx(Q.n_air_p, Q.zret/1000, 'b', 'LineWidth', 1.5);hold on
 MSIS = sum(Q.n_air_p(80:100));
 ret = sum(X.x(m+80:m+100));
 scaled_constant = MSIS./ret;
 air_OEM = X.x(m+1:2*m)*scaled_constant;
 semilogx(air_OEM, Q.zret/1000, 'r', 'LineWidth', 1.5);
 legend('MSIS', 'retrieved relative air density')
 ylim([12.6, 70])
 air_OEM_new = interp1(Q.zret/1000,  air_OEM, H/1000);
 constant_scaled = air_OEM_new(50)/tt(50);
 comp_scaled = ((constant_scaled)*tt);
 figure;
semilogx(comp_scaled, H, 'b',air_OEM_new, H, 'r')
avg = (comp_scaled+air_OEM_new)/2;
plot(((comp_scaled-air_OEM_new)./avg)*100, H/1000)

 %%%%%%%%%%%%
 