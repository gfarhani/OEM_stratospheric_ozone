function R  = makeParameterJacobians_OHP_test1(Q,x)
% makeParameterJacobians
%   pulled this out of makeR.m to speed up code, only do b parameter jacobians
%   at end
%	current b parameters:
%       
%       back scattering coefficent
%       Rayleigh scatter cross section
%       Ozone sctattrr cross section
%       Lidar Ratio
%       Angestrom coefficent
% -Usage-
%	R = makeRarameterJacobians(Q,x)
%
% -Inputs-
%	Q retrieval a priori information
%	x retrieved parameters
%
% -Outputs-
%	R Jacobians for b parameters
Zi = Q.zmsis;
Zj = Q.zret;
m = length(Q.zret);
xj = (interp1(Q.zret,x(1:m),Q.zmsis,'linear'));
nj = (interp1(Q.zret,x(m+1:2*m),Q.zmsis,'linear'));
%%%%%%%%%% Tansmission for ozone 308 %%%%%
gamma_norm_308 = exp(-2*x(end));
gamma_no3_308=exp(-2.*cumtrapz(xj.*Q.sigo3_308_low).*(Q.zmsis(2)-Q.zmsis(1)));
gamma_no3_308 = gamma_no3_308 .* gamma_norm_308;
%%%%%%%%%%% constant of transition for air density at 308 nm channel %%%%%
gamma_norm_air_308 = exp(-2*x(end-2));
gamma_nair_308=exp(-2.*cumtrapz(nj.*Q.sigma_Rayleigh_308).*(Q.zmsis(2)-Q.zmsis(1)));
gamma_nair_308 = gamma_nair_308 .* gamma_norm_air_308;
%%%%% backscattering coefficient for ir molecules
beta_mol_308 = (Q.sigma_Rayleigh_308 * Q.Pray/4/pi) * nj;
%%%%%%%%%%%%%%%%%%
term1_308 = (beta_mol_308./Q.sigma_Rayleigh_308);
term2_308 = -2*Q.Beta_aerosol_308.*cumtrapz(nj.*(Q.zmsis(2)-Q.zmsis(1)));
term3_308 = -2*beta_mol_308.*cumtrapz(nj.*(Q.zmsis(2)-Q.zmsis(1)));
comon1_308_med = (x(end-4).*gamma_nair_308.*gamma_no3_308.*Q.gamma_aerosol_308_low)./((Q.zmsis).^2);
dSMsigmaRay_308_med = comon1_308_med.*(term1_308+term2_308+term3_308);
dSMsigmaRay_308_med = dSMsigmaRay_308_med(Zi>=Q.bot_data_high);
for i = 1:length(dSMsigmaRay_308_med)
    if ~isempty(find(isnan(dSMsigmaRay_308_med(:,i))) == 1);
        dSMsigmaRay_308_med(:,i) = dSMsigmaRay_308_med(:,i-1);
    end
end
%%%%%%%%%%%%%%%%%%%%%
comon1_308_low = (x(end-5).*gamma_nair_308.*gamma_no3_308.*Q.gamma_aerosol_308_low)./((Q.zmsis).^2);
dSMsigmaRay_308_low = comon1_308_low.*(term1_308+term2_308+term3_308);
dSMsigmaRay_308_low  = dSMsigmaRay_308_low (Zi<=Q.top_data_low_308);
for i = 1:length(dSMsigmaRay_308_low)
    if ~isempty(find(isnan(dSMsigmaRay_308_low(:,i))) == 1);
      dSMsigmaRay_308_low(:,i) = dSMsigmaRay_308_low(:,i-1);
    end
end
%%%%%%%%%%%%%%%%%%%
%%%%% connstant of transmision for O3 at 353 nm cbannel %%%%%%%%%
gamma_norm_353 = exp(-2*x(end-1));
gamma_no3_353=exp(-2.*cumtrapz(xj.*Q.sigo3_353_low).*(Q.zmsis(2)-Q.zmsis(1)));
gamma_no3_353 = gamma_no3_353 .* gamma_norm_353;
%%%%%%%%%%%%%%%%%
%%%% constant of transition for air density at 353 nm channel %%%%%
gamma_norm_air_353 = exp(-2*x(end-3));
gamma_nair_353=exp(-2.*cumtrapz(nj.*Q.sigma_Rayleigh_353).*(Q.zmsis(2)-Q.zmsis(1)));
gamma_nair_353 = gamma_nair_353 .* gamma_norm_air_353;
%%%%%%%%%%%%%%%%
beta_mol_353 = (Q.sigma_Rayleigh_353 * Q.Pray/4/pi) * nj;
%%%%%%%%%%%%%%%%%%%%
term1_353 = (beta_mol_353./Q.sigma_Rayleigh_353);
term2_353 = -2*Q.Beta_aerosol_353.*cumtrapz(nj.*(Q.zmsis(2)-Q.zmsis(1)));
term3_353 = -2*beta_mol_353.*cumtrapz(nj.*(Q.zmsis(2)-Q.zmsis(1)));
comon1_353_med = (x(end-6).*gamma_nair_353.*gamma_no3_353.*Q.gamma_aerosol_353_low)./((Q.zmsis).^2);
dSMsigmaRay_353_med = comon1_353_med.*(term1_353+term2_353+term3_353);
dSMsigmaRay_353_med = dSMsigmaRay_353_med(Zi>=Q.bot_data_high);
for i = 1:length(dSMsigmaRay_353_med)
    if ~isempty(find(isnan(dSMsigmaRay_353_med(:,i))) == 1);
        dSMsigmaRay_353_med(:,i) = dSMsigmaRay_353_med(:,i-1);
    end
end
%%%%%%
term1_353 = (beta_mol_353./Q.sigma_Rayleigh_353);
term2_353 = -2*Q.Beta_aerosol_353.*cumtrapz(nj.*(Q.zmsis(2)-Q.zmsis(1)));
term3_353 = -2*beta_mol_353.*cumtrapz(nj.*(Q.zmsis(2)-Q.zmsis(1)));
comon1_353_low = (x(end-7).*gamma_nair_353.*gamma_no3_353.*Q.gamma_aerosol_353_low)./((Q.zmsis).^2);
dSMsigmaRay_353_low= comon1_353_low.*(term1_353+term2_353+term3_353);
dSMsigmaRay_353_low  = dSMsigmaRay_353_low (Zi>=Q.bot_data_low);
for i = 1:length(dSMsigmaRay_353_low)
    if ~isempty(find(isnan(dSMsigmaRay_353_low(:,i))) == 1);
      dSMsigmaRay_353_low(:,i) = dSMsigmaRay_353_low(:,i-1);
    end
end

KsigmaRay = [dSMsigmaRay_308_med'; dSMsigmaRay_308_low'; dSMsigmaRay_353_med'; dSMsigmaRay_353_low'];
R.KsigmaRay = diag(KsigmaRay);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Lidar Ratio only for 308 nm %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
term1_LR =  -2*cumtrapz(Q.Beta_aerosol_353.*(Q.zmsis(2)-Q.zmsis(1)));
term2_LR = Q.Beta_aerosol_353 + beta_mol_308;
dSM_LR_308 = comon1_308_med.*term1_LR.*term2_LR;

for i = 1:length(dSM_LR_308)
    if ~isempty(find(isnan(dSM_LR_308(:,i))) == 1);
        dSM_LR_308(:,i) = dSM_LR_308(:,i-1);
    end
end


dSL_LR_308 = comon1_308_low.*term1_LR.*term2_LR;

for i = 1:length(dSL_LR_308)
    if ~isempty(find(isnan(dSL_LR_308(:,i))) == 1);
        dSL_LR_308(:,i) = dSL_LR_308(:,i-1);
    end
end
KLR = [dSM_LR_308'; dSL_LR_308'];
R.KLR = diag(KLR);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% angestrome coeficient only for 353  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
term_alpha = -2*Q.LR.*cumtrapz(Q.Beta_aerosol_353.*(Q.zmsis(2)-Q.zmsis(1))).*log(353/308).*((353/308)^(Q.angestrom));
term1_alpha_med = comon1_308_med.*Q.Beta_aerosol_308.*log(353/308);
term2_alpha_med = comon1_308_med.*term_alpha.*beta_mol_308;
term3_alpha_med = comon1_308_med.*term_alpha.*Q.Beta_aerosol_308;

dSM_anges_308 =term1_alpha_med+term2_alpha_med+term3_alpha_med;

for i = 1:length(dSM_anges_308)
    if ~isempty(find(isnan(dSM_anges_308(:,i))) == 1);
        dSM_anges_308(:,i) = dSM_anges_308(:,i-1);
    end
end

term1_alpha_low = comon1_308_low.*Q.Beta_aerosol_308.*log(353/308);
term2_alpha_low = comon1_308_low.*term_alpha.*beta_mol_308;
term3_alpha_low = comon1_308_low.*term_alpha.*Q.Beta_aerosol_308;

dSL_anges_308 =term1_alpha_low+term2_alpha_low+term3_alpha_low;

for i = 1:length(dSL_anges_308)
    if ~isempty(find(isnan(dSL_anges_308(:,i))) == 1);
        dSL_anges_308(:,i) = dSL_anges_308(:,i-1);
    end
end

Kanges = [dSM_anges_308'; dSL_anges_308'];
R.Kanges = diag(Kanges);

%%%%%%%%%%%%%%%%%%%%%%% sigma O3 which is numerical %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is the real forward model
[S_med_308, S_low_308, S_med_353, S_low_353]=forwardmodel_test7_308_353(Q,x);
gamma_norm_308 = exp(-2*x(end));
dn = 1.e-23; % 1e-5
sigo3j = Q.sigo3_308_low + Q.sigo3_308_low *dn;
gamma_no3_308j=exp(-2.*cumtrapz(xj.*sigo3j).*(Q.zmsis(2)-Q.zmsis(1)));
gamma_no3_308j = gamma_no3_308j .* gamma_norm_308;

gamma_norm_353 = exp(-2*x(end-1));
dn1 = 1.e-19; % 1e-5
sigo3_353j = Q.sigo3_353_low + dn1;
gamma_no3_353j=exp(-2.*cumtrapz(xj.*sigo3_353j).*(Q.zmsis(2)-Q.zmsis(1)));
gamma_no3_353j = gamma_no3_353j .* gamma_norm_353;


SM_308j= ((x(end-4))*(beta_mol_308+Q.Beta_aerosol_308).*gamma_nair_308.*gamma_no3_308j.*Q.gamma_aerosol_308_low)./((Q.zmsis).^2) ;
S_med_308_truej = SM_308j+(x(end-11));

SL_308j= ((x(end-5)).*(beta_mol_308+Q.Beta_aerosol_308).*gamma_nair_308.*gamma_no3_308j.*Q.gamma_aerosol_308_low)./((Q.zmsis).^2) ;
S_low_308_truej= SL_308j + exp(x(end-10)+(x(end-9)*Q.zmsis/1000))+ x(end-8);

SM_353j= ((x(end-6)).*(beta_mol_353+Q.Beta_aerosol_353).*gamma_nair_353.*gamma_no3_353j.*Q.gamma_aerosol_353_low)./((Q.zmsis).^2) ;
S_med_353_truej = SM_353j + x(end-12);
SL_353j= ((x(end-7)).*(beta_mol_353+Q.Beta_aerosol_353).*gamma_nair_353.*gamma_no3_353j.*Q.gamma_aerosol_353_low)./((Q.zmsis).^2) ;
S_low_353_truej = SL_353j + x(end-13);

dS_M_308_o3 = (S_med_308_truej - S_med_308) ./ dn;
dS_L_308_o3 = (S_low_308_truej - S_low_308) ./ dn;
dS_M_353_o3 = (S_med_353_truej - S_med_353) ./ dn;
dS_L_353_o3 = (S_low_353_truej - S_low_353) ./ dn;

dS_M_308_o3=  dS_M_308_o3(Zi>=Q.bot_data_high);
dS_L_308_o3  = dS_L_308_o3(Zi<=Q.top_data_low_308);
dS_M_353_o3 = dS_M_353_o3(Zi>=Q.bot_data_high);
dS_L_353_o3  = dS_L_353_o3 (Zi>=Q.bot_data_low);
for i = 1:length(dS_M_308_o3)
    if ~isempty(find(isnan(dS_M_308_o3(:,i))) == 1);
        dS_M_308_o3(:,i) = dS_M_308_o3(:,i-1);
    end
end

for i = 1:length(dS_L_308_o3)
    if ~isempty(find(isnan(dS_L_308_o3(:,i))) == 1);
        dS_L_308_o3(:,i) = dS_L_308_o3(:,i-1);
    end
end

for i = 1:length(dS_M_353_o3)
    if ~isempty(find(isnan(dS_M_353_o3(:,i))) == 1);
        dS_M_353_o3(:,i) = dS_M_353_o3(:,i-1);
    end
end

for i = 1:length(dS_L_353_o3 )
    if ~isempty(find(isnan(dS_L_353_o3 (:,i))) == 1);
        dS_L_353_o3 (:,i) = dS_L_353_o3 (:,i-1);
    end
end


Ksigma_o3 = [dS_M_308_o3'; dS_L_308_o3'; dS_M_353_o3'; dS_L_353_o3'];
R.Ksigma_o3 = diag(Ksigma_o3);

sigo3_T_high=  Q.delta_sigma(Zi>=Q.bot_data_high);
sigo3_T_low  = Q.delta_sigma(Zi<=Q.top_data_low_308);

for i = 1:length(sigo3_T_high)
    if ~isempty(find(isnan(sigo3_T_high(1,i))) == 1);
        sigo3_T_high(1,i) = 1e-26;
    end
end

for i = 1:length(sigo3_T_low)
    if ~isempty(find(isnan(sigo3_T_low(1,i))) == 1);
        sigo3_T_low(1,i) = 1e-26;
    end
end
Ksigma_o3_T = [dS_M_308_o3'.*sigo3_T_high'; dS_L_308_o3'.*sigo3_T_low'];
R.Ksigma_o3_T = diag(Ksigma_o3_T);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% aresol backscattring coefficent %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beta_aerosol_353 = 1e-8*exp(-(Q.alt_data - 10000)/15000);
gamma_aerosol_353 = exp(-2.*cumtrapz((beta_aerosol_353)*Q.LR).*(Q.zmsis(2)-Q.zmsis(1)));
beta_aerosol_308 = beta_aerosol_353 * (308/353)^(-Q.angestrom);
gamma_aerosol_308 = exp(-2.*cumtrapz((beta_aerosol_308)*Q.LR).*(Q.zmsis(2)-Q.zmsis(1)));

dn2_353 = 1.e-5*beta_aerosol_353;
beta_aerosol_353j = beta_aerosol_353 +dn2_353;
gamma_aerosol_353j = exp(-2.*cumtrapz((beta_aerosol_353j)*Q.LR).*(Q.zmsis(2)-Q.zmsis(1)));
beta_aerosol_308j = beta_aerosol_353j * (308/353)^(-Q.angestrom);
gamma_aerosol_308j = exp(-2.*cumtrapz((beta_aerosol_308j)*Q.LR).*(Q.zmsis(2)-Q.zmsis(1)));

gamma_aerosol_353j_low = gamma_aerosol_353j(Q.A:Q.B);
gamma_aerosol_308j_low = gamma_aerosol_308j(Q.A:Q.B);
% % gamma_aerosol_353j_med = gamma_aerosol_353j(Q.C:Q.B);
% % gamma_aerosol_308j_med = gamma_aerosol_308j(Q.C:Q.B);

z0 = 755:300:Q.bot_ret; l =length(z0);
beta_aerosol_3530 = beta_aerosol_353j(1:l);
odnorm_aerosol_353 = trapz(z0, beta_aerosol_3530*Q.LR);
beta_aerosol_3080 = beta_aerosol_308j(1:l);
odnorm_aerosol_308 = trapz(z0, beta_aerosol_3080*Q.LR);

gamma_norm_aerosol_353 = exp(-odnorm_aerosol_353 );
gamma_aerosol_353_j = gamma_aerosol_353j_low.*gamma_norm_aerosol_353;
gamma_norm_aerosol_308 = exp(-odnorm_aerosol_308 );
gamma_aerosol_308_j = gamma_aerosol_308j_low.*gamma_norm_aerosol_308;

beta_aerosol_308j_low = beta_aerosol_308j(Q.A:Q.B);
beta_aerosol_353j_low = beta_aerosol_353j(Q.A:Q.B);
SM_308j= ((x(end-4))*(beta_mol_308+ beta_aerosol_308j_low).*gamma_nair_308.*gamma_no3_308.*gamma_aerosol_308_j)./((Q.zmsis).^2) ;
S_med_308_truej = SM_308j+(x(end-11));

SL_308j= ((x(end-5)).*(beta_mol_308+beta_aerosol_308j_low).*gamma_nair_308.*gamma_no3_308.*gamma_aerosol_308_j)./((Q.zmsis).^2) ;
S_low_308_truej= SL_308j + exp(x(end-10)+(x(end-9)*Q.zmsis/1000))+ x(end-8);

SM_353j= ((x(end-6)).*(beta_mol_353+beta_aerosol_353j_low).*gamma_nair_353.*gamma_no3_353.*gamma_aerosol_353_j)./((Q.zmsis).^2) ;
S_med_353_truej = SM_353j + x(end-12);
SL_353j= ((x(end-7)).*(beta_mol_353+beta_aerosol_353j_low).*gamma_nair_353.*gamma_no3_353.*gamma_aerosol_353_j)./((Q.zmsis).^2) ;
S_low_353_truej = SL_353j + x(end-13);

dS_M_308_aerosol = (S_med_308_truej - S_med_308) ./ dn;
dS_L_308_aerosol = (S_low_308_truej - S_low_308) ./ dn;
dS_M_353_aerosol = (S_med_353_truej - S_med_353) ./ dn;
dS_L_353_aerosol = (S_low_353_truej - S_low_353) ./ dn;

dS_M_308_aerosol=  dS_M_308_aerosol(Zi>=Q.bot_data_high);
dS_L_308_aerosol = dS_L_308_aerosol(Zi<=Q.top_data_low_308);
dS_M_353_aerosol= dS_M_353_aerosol(Zi>=Q.bot_data_high);
dS_L_353_aerosol = dS_L_353_aerosol(Zi>=Q.bot_data_low);


for i = 1:length(dS_M_308_aerosol)
    if ~isempty(find(isnan(dS_M_308_aerosol(:,i))) == 1);
        dS_M_308_aerosol(:,i) = dS_M_308_aerosol(:,i-1);
    end
end

for i = 1:length(dS_L_308_aerosol)
    if ~isempty(find(isnan(dS_L_308_aerosol(:,i))) == 1);
        dS_L_308_aerosol(:,i) = dS_L_308_aerosol(:,i-1);
    end
end

for i = 1:length(dS_M_353_aerosol)
    if ~isempty(find(isnan(dS_M_353_aerosol(:,i))) == 1);
        dS_M_353_aerosol(:,i) = dS_M_353_aerosol(:,i-1);
    end
end

for i = 1:length(dS_L_353_aerosol)
    if ~isempty(find(isnan(dS_L_353_aerosol(:,i))) == 1);
        dS_L_353_aerosol (:,i) = dS_L_353_aerosol (:,i-1);
    end
end

Ksigma_aerosol = [dS_M_308_aerosol'; dS_L_308_aerosol'; dS_M_353_aerosol'; dS_L_353_aerosol'];
R.Ksigma_aerosol = diag(Ksigma_aerosol);

%%%%%%%%%%%%%%%%%%%%%%% temperaure %%%%%%%%%%%%%%%
%This is the real forward model
[S_med_308, S_low_308, S_med_353, S_low_353]=forwardmodel_test7_308_353(Q,x);
gamma_norm_308 = exp(-2*x(end));
% dn = 1.e-23; % 1e-5
sigo3j = Q.sigo3_308_lowj;
gamma_no3_308j=exp(-2.*cumtrapz(xj.*sigo3j).*(Q.zmsis(2)-Q.zmsis(1)));
gamma_no3_308j = gamma_no3_308j .* gamma_norm_308;



SM_308j= ((x(end-4))*(beta_mol_308+Q.Beta_aerosol_308).*gamma_nair_308.*gamma_no3_308j.*Q.gamma_aerosol_308_low)./((Q.zmsis).^2) ;
S_med_308_truej = SM_308j+(x(end-11));

SL_308j= ((x(end-5)).*(beta_mol_308+Q.Beta_aerosol_308).*gamma_nair_308.*gamma_no3_308j.*Q.gamma_aerosol_308_low)./((Q.zmsis).^2) ;
S_low_308_truej= SL_308j + exp(x(end-10)+(x(end-9)*Q.zmsis/1000))+ x(end-8);


dS_M_308_o3 = (S_med_308_truej - S_med_308) ./ dn;
dS_L_308_o3 = (S_low_308_truej - S_low_308) ./ dn;


dS_M_308_o3=  dS_M_308_o3(Zi>=Q.bot_data_high);
dS_L_308_o3  = dS_L_308_o3(Zi<=Q.top_data_low_308);

for i = 1:length(dS_M_308_o3)
    if ~isempty(find(isnan(dS_M_308_o3(:,i))) == 1);
        dS_M_308_o3(:,i) = dS_M_308_o3(:,i-1);
    end
end

for i = 1:length(dS_L_308_o3)
    if ~isempty(find(isnan(dS_L_308_o3(:,i))) == 1);
        dS_L_308_o3(:,i) = dS_L_308_o3(:,i-1);
    end
end




Ksigma_sigmao3_T = [dS_M_308_o3'; dS_L_308_o3'];
R.Ksigma_sigmao3_T = diag(Ksigma_sigmao3_T);

return;



