function [S_med_308, S_low_308, S_med_353, S_low_353]=forwardmodel_test7_308_353(Q,x)
%length of retrieval 
m = length(Q.zret);
%x(1:m) is the ozone retrieval
%x(m+1:2*m) is the air density retrieval 
%x((2*m)+1:3*m) is the backscattering 
xj = (interp1(Q.zret,x(1:m),Q.zmsis,'linear'));
nj = (interp1(Q.zret,x(m+1:2*m),Q.zmsis,'linear'));
% to build the optical depth befor the zret 
%%%%%%%%%%% constant of transmission before the retrieval %%%%%%%%%%%%%%% 
%%%%% connstant of transmision for O3 at 308  nm channel %%%%%%%%%
% % z0 = 640:Q.delta_data:Q.zmsis(1);
% % no30 = xj(1) .*ones(size(z0)); 
gamma_norm_308 = exp(-2*x(end));
gamma_no3_308=exp(-2.*cumtrapz(xj.*Q.sigo3_308_low).*(Q.zmsis(2)-Q.zmsis(1)));
gamma_no3_308 = gamma_no3_308 .* gamma_norm_308;
%%%%% connstant of transmision for O3 at 353 nm cbannel %%%%%%%%%
gamma_norm_353 = exp(-2*x(end-1));
gamma_no3_353=exp(-2.*cumtrapz(xj.*Q.sigo3_353_low).*(Q.zmsis(2)-Q.zmsis(1)));
gamma_no3_353 = gamma_no3_353 .* gamma_norm_353;
%%%%%%%%%%% constant of transition for air density at 308 nm channel %%%%%
gamma_norm_air_308 = exp(-2*x(end-2));
gamma_nair_308=exp(-2.*cumtrapz((Q.sigma_Rayleigh_308) * nj).*(Q.zmsis(2)-Q.zmsis(1)));
gamma_nair_308 = gamma_nair_308 .* gamma_norm_air_308;
%%%% constant of transition for air density at 353 nm channel %%%%%
gamma_norm_air_353 = exp(-2*x(end-3));
gamma_nair_353=exp(-2.*cumtrapz((Q.sigma_Rayleigh_353) * nj).*(Q.zmsis(2)-Q.zmsis(1)));
gamma_nair_353 = gamma_nair_353 .* gamma_norm_air_353;

%%%%% beta for 353 nm channel %%%%%
beta_mol_353 = (Q.sigma_Rayleigh_353 *(Q.Pray/4/pi)) * nj;
%%%%% beta for 308 nm channel %%%%% Q.Pray/4/pi
beta_mol_308 = (Q.sigma_Rayleigh_308 * (Q.Pray/4/pi)) * nj;
%Here is the forwrad model for the Hich channel 308
SM_308= ((x(end-4))*(beta_mol_308).*gamma_nair_308.*gamma_no3_308)./((Q.zmsis).^2) ;
% S_med_308= SM +Q.bg_med_308;
S_med_308_true = SM_308+(x(end-11));
%S_med_308 = SM_308 + exp(x(end-13)+(x(end-12)*Q.zmsis/1000))+ x(end-11);
% S_med_308 = SM_308 + x(end-3)+exp((-Q.zmsis/1000)*x(end-2));
SL_308= ((x(end-5)).*(beta_mol_308).*gamma_nair_308.*gamma_no3_308)./((Q.zmsis).^2) ;
% S_med_308= SM +Q.bg_med_308;
% S_low_308 = SL_308 +x(end);
%S_low_308 = SL_308+(x(end-1)+((Q.zmsis/1000)*x(end)));
S_low_308_true= SL_308 + exp(x(end-10)+(x(end-9)*Q.zmsis/1000))+ x(end-8);
SM_353= ((x(end-6)).*(beta_mol_353).*gamma_nair_353.*gamma_no3_353)./((Q.zmsis).^2) ;
S_med_353_true = SM_353 + x(end-12);
SL_353= ((x(end-7)).*(beta_mol_353).*gamma_nair_353.*gamma_no3_353)./((Q.zmsis).^2) ;
S_low_353_true = SL_353 + x(end-13);
% % % % 
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% adding dead time %%%%%%%%%%%%
% % x = (Q.tau)*1e-9;
% % S_med_308 = 1+ ((1-x)*S_med_308_true-1).*(exp(-x*S_med_308_true));
% % S_med_353 = 1+ ((1-x)*S_med_353_true-1).*(exp(-x*S_med_308_true));
% % S_low_308 = 1+ ((1-x)*S_low_308_true-1).*(exp(-x*S_low_308_true));
% % S_low_353 = 1+ ((1-x)*S_low_353_true-1).*(exp(-x*S_low_308_true));
% % % % % % 
% % % % % deltaTime = 80;
% % % % % Rate_308 = 100;
% % % % % Rate_353 = 50;
clight = 3e8;
y2HzRaw_308 = clight ./ (2.*( 1608326 ).*Q.delta_data);
y2HzRaw_353 = clight ./ (2.*(  728143 ).*Q.delta_data);
%%%%%%%%
% % % % % % % % % % % % % % 
S_med_308_HZ = y2HzRaw_308 * S_med_308_true;
S_med_353_HZ = y2HzRaw_353 * S_med_353_true;
S_low_308_HZ = y2HzRaw_308 * S_low_308_true;
S_low_353_HZ = y2HzRaw_353 * S_low_353_true;

% % % % % % % 
S_med_308 = (S_med_308_true).*exp(-S_med_308_HZ .* x(end-15)*1e-9);
S_med_353 = (S_med_353_true).*exp(-S_med_353_HZ .* x(end-16)*1e-9);
S_low_353 = (S_low_353_true).*exp(-S_low_353_HZ .* x(end-14)*1e-9);
S_low_308 = (S_low_308_true).*exp(-S_low_308_HZ .* x(end-17)*1e-9);
% % % % 
% % % S_med_308 = (S_med_308_HZ./(1+(x(end-14)*1e-9)*S_med_308_HZ))/y2HzRaw_308;
% % % S_low_308 = S_low_308_HZ./(1+(x(end-14)*1e-9)*S_low_308_HZ)/y2HzRaw_308;
% % % S_med_353 = S_med_353_HZ./(1+(x(end-14)*1e-9)*S_med_353_HZ)/y2HzRaw_353;
% % % S_low_353 = S_low_353_HZ./(1+(x(end-14)*1e-9)*S_low_353_HZ)/y2HzRaw_353;
% % % % 
% % % % S_med_308 = S_med_308_Obs*1200000;
% % % % S_low_308 = S_low_308_Obs*1200000;
% % % % S_med_353 = S_med_353_Obs*600000;
% % % % S_low_353 = S_low_353_Obs*600000;
% % % % 
% % % % S_med_308 = S_med_308_HZ_Obs./y2HzRaw ;
% % % % S_med_353 = S_med_353_HZ_Obs./y2HzRaw ;
% % % % S_low_308 = S_low_308_HZ_Obs./y2HzRaw ;
% % % % S_low_353 = S_low_353_HZ_Obs./y2HzRaw ;

for i = 1:length(S_low_353)
    if ~isempty(find(isnan(S_low_353(:,i))) == 1);
        S_low_353(:,i) = S_low_353(:,i-2);
    end
end
for i = 1:length(S_low_308)
    if ~isempty(find(isnan(S_low_308(:,i))) == 1);
        S_low_308(:,i) = S_low_308(:,i-1);
    end
end
for i = 1:length(S_med_353)
    if ~isempty(find(isnan(S_med_353(:,i))) == 1);
        S_med_353(:,i) = S_med_353(:,i-1);
    end
end

for i = 1:length(S_med_308)
    if ~isempty(find(isnan(S_med_308(:,i))) == 1);
        S_med_308(:,i) = S_med_308(:,i-1);
    end
end

