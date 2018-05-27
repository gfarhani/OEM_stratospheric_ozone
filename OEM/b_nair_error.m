%%%%%%%%%%%%%%%%%%%%%% errors %%%%%%%%%%
figure;
 X.vmr = X.x(m+1:2*m);
 obsErrO3 = X.eo(m+1:2*m) ./ X.vmr;
 plot(obsErrO3(2:end-1).*100,Q.zret(2:end-1)./1000);hold on;
 

R  = makeParameterJacobians_OHP_test1(Q,x);
dfacSigmaR = 0.003; % ISSI recommend
SsigmaR_308 = (dfacSigmaR.*Q.sigma_Rayleigh_308).^2;
Sxsigma_308 = X.G(:,1:m1+m2)*R.KsigmaRay(:,1:m1+m2)'*SsigmaR_308*R.KsigmaRay(:,1:m1+m2)*X.G(:,1:m1+m2)';
sigmaRayErrq_308 = sqrt(diag(Sxsigma_308(m+1:2*m,m+1:2*m))) ./ X.vmr;
plot(sigmaRayErrq_308(2:end).*100,Q.zret(2:end)/1000,'--', 'LineWidth', 1.5)

SsigmaR_353 = (dfacSigmaR.*Q.sigma_Rayleigh_353).^2;
Sxsigma_353 = X.G(:,m1+m2+1:end)*R.KsigmaRay(:,m1+m2+1:end)'*SsigmaR_353*R.KsigmaRay(:,m1+m2+1:end)*X.G(:,m1+m2+1:end)';
sigmaRayErrq_353 = sqrt(diag(Sxsigma_353(m+1:2*m,m+1:2*m))) ./ X.vmr;
hold all
plot(sigmaRayErrq_353(2:end).*100,Q.zret(2:end)/1000,'--', 'LineWidth', 1.5)

xlabel('Ozone Density Uncertainty (%)')
ylabel('Altitude (km)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Lidar Ration only for on-line channel %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % % dfacLR = 0.03; % ISSI recommend
% % % LR = (dfacLR.*Q.LR).^2;
% % % SxLR_308 = X.G(:,1:m1+m2)*R.KLR(:,1:m1+m2)'*LR*R.KLR(:,1:m1+m2)*X.G(:,1:m1+m2)';
% % % sigmaLR_308 = sqrt(diag(SxLR_308(m+1:2*m,m+1:2*m))) ./ X.vmr;
% % % plot(sigmaLR_308(2:end).*100,Q.zret(2:end)/1000,'--', 'LineWidth', 1.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% angestrm only for 308 nm %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % dfacanges= 0.1; % ISSI recommend
% % % Q.angestrom = (dfacanges.*Q.angestrom).^2;
% % % Sxanges_308 = X.G(:,1:m1+m2)*R.Kanges(:,1:m1+m2)'*Q.angestrom*R.Kanges (:,1:m1+m2)*X.G(:,1:m1+m2)';
% % % sigma_ang_308 = sqrt(diag(Sxanges_308(m+1:2*m,m+1:2*m))) ./ X.vmr;
% % % plot(sigma_ang_308(2:end).*100,Q.zret(2:end)/1000,'--', 'LineWidth', 1.5)

%%%%%%%%%%%%%%%%%%%%%%%%%% sigma O3 %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dfacSimaO3_308 = 0.02; % ISSI recommend
sigo3_308 = [Q.sigo3_308_med, Q.sigo3_b_parameter];



SsigmaO3_308 = (dfacSimaO3_308.*sigo3_308).^2;
SsigmaO3_308 = diag(SsigmaO3_308);
Sxsigma_O3_308 = X.G(:,1:m1+m2)*R.Ksigma_o3(1:m1+m2,1:m1+m2)*SsigmaO3_308*R.Ksigma_o3(1:m1+m2,1:m1+m2)'*X.G(:,1:m1+m2)';
sigmaErrq_O3_308 = sqrt(diag(Sxsigma_O3_308(1:m,1:m))) ./ X.vmr;
plot(sigmaErrq_O3_308(2:end).*100,Q.zret(2:end)/1000,'--', 'LineWidth', 1.5)

dfacSimaO3_353 = 0.2; % ISSI recommend
sigo3_353 = 1.01611270000000e-26;
SsigmaO3_353 = (dfacSimaO3_353.*sigo3_353).^2;

Sxsigma_o3_353 = X.G(:,m1+m2+1:end)*R.Ksigma_o3(:,m1+m2+1:end)'*SsigmaO3_353*R.Ksigma_o3(:,m1+m2+1:end)*X.G(:,m1+m2+1:end)';
sigma_o3_Errq_353 = sqrt(diag(Sxsigma_o3_353(m+1:2*m,m+1:2*m))) ./ X.vmr;
hold all
total = sigma_o3_Errq_353+sigmaRayErrq_353+sigmaRayErrq_308+obsErrO3;
plot(total(2:end).*100,Q.zret(2:end)/1000, 'LineWidth', 3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% aerosol back scattering coefecient%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % dfacSima_beta_308 = 0.2; % ISSI recommend
% % sig_beta_308 = [Q.Beta_aerosol_308_med, Q.Beta_aerosol_308_low];
% % Ssigma_beta_308 = (dfacSima_beta_308.*sig_beta_308).^2;
% % Ssigma_beta_308 = diag(Ssigma_beta_308);
% % Sxsigma_beta_308 = X.G(:,1:m1+m2)*R.Ksigma_o3(1:m1+m2,1:m1+m2)*Ssigma_beta_308*R.Ksigma_o3(1:m1+m2,1:m1+m2)'*X.G(:,1:m1+m2)';
% % sigmaErrq_beta_308 = sqrt(diag(Sxsigma_beta_308(m+1:2*m,m+1:2*m))) ./ X.vmr;
% % plot(sigmaErrq_beta_308(2:end).*100,Q.zret(2:end)/1000,'--', 'LineWidth', 1.5)

hleg = legend('Statistical', '\sigma_{Rayleigh_308}',  '\sigma_{Rayleigh_353}', ....
     '\sigma_{o3_{308}}', 'Total');
set(hleg,'FontSize',8,'Box','off');
pltx = get(gca,'XLim');
plot(pltx,[Q.zret(fini) Q.zret(fini)]./1000,'k--')
ylim([1 50])
xlim([0 .7])