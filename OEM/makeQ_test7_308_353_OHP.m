% %%%%%%%%%%% Testtemp Cov is used for O3 %%%%%
% lengthcT = 300;
% % %  lc(Zj<=410000) = 900;
% % %  lc(Zj>40000) = 150;X
% % % 
% Tfac = ones(1, length(Zj));
% Tfac(1,1:11) = .5*x(1:11);
% Tfac(1,12:25) = .5*x(12:25);
% Tfac(1,25:end) = .2*x(25:end);
% %%%%%%%%%
% %%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Air density
% lengthcT_low = 300;
% lengthcT_med = 900;
% lengthcT_high = 900;
% % %  lc(Zj<=40000) = 900;
% % %  lc(Zj>40000) = 150;
% % % 
% Tfac = ones(1, length(Zj));
% Tfac(1,1:11) = .15*x(1:11);
% Tfac(1,11:17) = .15*x(11:17);
% % Tfac(1,17:26) = .3*x(17:26);
% Tfac(17:end) = .15*x(17:end);
 
%%%%%%%%%%%%%%%
function [Q, y_real, Sy_real] = makeQ_test7_308_353_OHP( )
 
[nmsis, Tmsis, zmsis]=ATMOSPHERIC( ); 
% MSIS data as a priori for the air density1e 
zmsis1 = zmsis; % SI units
nmsis1 = nmsis;  %SI units
Tmsis1 = Tmsis;  % SI units
%%%%%%%%%%%
% % % % % % % test = load('sigtestOHP_Ghazal.mat');
% % % % % % % data = test.sigmatrix(90:100, :, :);
% % % % % % % real_med_308 = data(5,:,1)*886170;
% % % % % % % real_med_353 = data(5,:,2)*886170;
% % % % % % % real_low_308 = data(5,:,3)*886107;
% % % % % % % real_low_353 = data(5,:,4)*886107;
% % % % % % % alt = (test.alt)'*1000;
%%%%%%%%%%%
OHP_2017 = load('SignalOHP2017.mat');
 
alt = (OHP_2017.alt)'*1000;
data = load('OHP_14072017.mat');
real_med_308 = data.Sig_Matrix(:,1)* 1608326 ;
real_med_353 = data.Sig_Matrix(:,2)*  804551  ;
real_low_308 = data.Sig_Matrix(:,3)*  1608326 ;
real_low_353 = data.Sig_Matrix(:,4)*  804551 ;
% % % % % % 
% % % [real_med_308, alt] = coadd(real_med_308, alt_OHP, 4);
% % % [real_low_308, alt] = coadd(real_low_308, alt_OHP, 4);
% % % [real_med_353, alt] = coadd(real_med_353, alt_OHP, 4);
% % % [real_low_353, alt] = coadd(real_low_353, alt_OHP, 4);
% % % % % % 
real_med_308= real_med_308(1:end-25)';
real_med_353 = real_med_353(1:end-25)';
real_low_308 = real_low_308(1:end-25)';
real_low_353 = real_low_353(1:end-25)';
alt = alt(1:end-25);
% % % data = OHP_2017.sigmatrix(15, :, :);
% % % real_med_308 = data(:,:,1)*OHP_2017.ntirtab(12,1);
% % % real_med_353 = data(:,:,2)*OHP_2017.ntirtab(12,2);
% % % real_low_308 = data(:,:,3)*OHP_2017.ntirtab(12,1);
% % % real_low_353 = data(:,:,4)*OHP_2017.ntirtab(12,2);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmsis = interp1(zmsis1, nmsis, alt');
Tmsis = interp1(zmsis1, Tmsis, alt');
zmsis = alt';  
alt_data = alt';
Q.alt_data = alt_data;
%%%%%%
% Rayleigh cross section for 308 and 353 nm wavelength the units for sigma
% is m^2
lambda_308 = 308e-9; lambda_353 = 353e-9; % convert in m 
ro_n=0.0301;%depol factor
N_s=2.547e25;% Molecular density(m-3)
n_s_308 =(5791817/(238.0185-(1/(lambda_308*1e6))^2)+167909/(57.362-(1/(lambda_308*1e6))^2))*1e-8+1;
n_s_353 =(5791817/(238.0185-(1/(lambda_353*1e6))^2)+167909/(57.362-(1/(lambda_353*1e6))^2))*1e-8+1;
sigma_Rayleigh_308 = ((24*(pi^3)*((n_s_308^2-1)^2))/((lambda_308^4)*(N_s^2)*(n_s_308^2+2)^2))*((6+3*ro_n)/(6-7*ro_n));
sigma_Rayleigh_353 = ((24*(pi^3)*((n_s_353^2-1)^2))/((lambda_353^4)*(N_s^2)*(n_s_353^2+2)^2))*((6+3*ro_n)/(6-7*ro_n));
% Total voume scattering coeffcient 
sigma_Rayleigh_308 = sigma_Rayleigh_308;
beta_mol_tot_353 = sigma_Rayleigh_353 * nmsis;
beta_mol_tot_308 = sigma_Rayleigh_308 * nmsis;
% Backscattering coefficient (angle of theta = pi)
gamma=ro_n/(2-ro_n);
Pray=3/4/(1+2*gamma)*(1+3*gamma+(1-gamma)*cos(pi())^2);
beta_mol_308 =beta_mol_tot_308/4/pi*Pray;
beta_mol_353 =beta_mol_tot_353/4/pi*Pray;
%%%%% transmission for air %%%%%
gamma_air_308 = exp(cumtrapz(zmsis,-2*sigma_Rayleigh_308.*nmsis));
gamma_air_353 = exp(cumtrapz(zmsis,-2*sigma_Rayleigh_353.*nmsis));
%Ozone Cross Section for 308 nm and 353 nm for the latter wavelength it is
%a constant value 
[sigo3_308 status, message]=Eureka_calc_sec_eff(Tmsis,308);
sigo3_308 = sigo3_308*1e-4; % unit conversion to [m^2]
Tmsisj = Tmsis+Tmsis*.001;
[sigo3_308j status, message]=Eureka_calc_sec_eff(Tmsisj,308);
sigo3_308j = sigo3_308j*1e-4; % unit conversion to [m^2]
delta_sigma = (sigo3_308j-sigo3_308)/.001;
Q.delta_sigma = delta_sigma;
[sigo3_353 status, message]=Eureka_calc_sec_eff(Tmsis,353);
sigo3_353 = sigo3_353*1e-4; % unit conversion to [m^2]
%sigo3_308 = 1.26E-23;
% a priori value for the ozone concentration
% % % O3_dens = 5.81e12*exp(-3.9e-3*(zmsis/1000-15.8).^2)*1e6;
OHP_sonde_compare;
% % 
A = load('O3data2.mat');
Oz_dens = A.O3den;
oz_alt = A.zO3den;
oz_alt = [1000:1030:5000+A.zO3den(end)];
 
% % % % % % % % O3 = O3;
Oz_dens = A.O3den(1:end);
% % % % % %  O3_dens = interp1(oz_alt ,O3, zmsis);
% % % % % % A = load('O3data2.mat');
% % % % Ozone = load('O3_new.mat');
% % % % Oz_dens = Ozone.unnamed1(:,1);
% % % % oz_alt = A.zO3den;
% % % O3_dens = interp1(oz_alt ,Oz_dens, zmsis);
A = load('OHP_climatology.mat');
Oz_dens = A.unnamed1(:,2);
oz_alt = A.unnamed1(:,1);
o3_sonde = interp1(alts*1000, conco3s, oz_alt);
O3 = [o3_sonde(1:15)*1e6; 1*Oz_dens(16:end)];
 O3_dens = interp1(oz_alt ,Oz_dens  , zmsis);
for i = 1:length(O3_dens)
    if ~isempty(find(isnan(O3_dens(1,i))) == 1);
        O3_dens(1,i) = 0;
    end
end
%% calculating the ozone transmission 
gamma_no3_308=exp(-2.*cumtrapz(O3_dens.*sigo3_308).*(zmsis(2)-zmsis(1)));
gamma_no3_353=exp(-2.*cumtrapz(O3_dens.*sigo3_353).*(zmsis(2)-zmsis(1)));
%%%%%%%%%%%%%%%%%%%%
% % %  beta_aero_532 = ncread('ahsrl_20060201T0000_20060301T0000_180s_45m.nc', 'profile_beta_a_backscat');
% % %  altitude = ncread('ahsrl_20060201T0000_20060301T0000_180s_45m.nc', 'altitude');
% % % % %  m = 28; % Length of the month
% % % % %  d = 26;% day of the month
% % % % %  t = length(beta_aero)./m;
% % % % %  start_point = t*d;
% % % % %  end_point = start_point + t;
% % % % %  % here is the portion of the whole beta profiles that is related to that
% % % % %  % specific night of interest
% % % % %  beta_day = beta_aero(:,start_point+1:end_point);
% % % % % % % for i = 1:length(beta_day(:,1));
% % % % % % % for j = 1:length(beta_day(1,:));
% % % % % % % if ~isempty(find(isnan(beta_day(i,j))) == 1);
% % % % % % % beta_day(i,j) = 0;
% % % % % % % end
% % % % % % % end
% % % % % % % end
% % % % % % for i = 1:length(beta_day)
% % % % % %     if beta_day(:,i) <0
% % % % % %         beta_day(:,i) = 0;
% % % % % %     end
% % % % % % end
% % % % % for i=1:length(beta_day(:,1))
% % % % % beta_532(i) = nanmean(beta_day(:,i));
% % % % % end
% % % angestrom = 1;
% % % lambda_532 = 532e-9;
% % % beta_aero_308 = beta_aero_532 * (lambda_308/lambda_532)^(-angestrom);
% % % beta_aero_353 = beta_aero_532 * (lambda_353/lambda_532)^(-angestrom);
% % % LR = 50;
% % % Beta_aerosol_308 = interp1(altitude, beta_aero_308, zmsis);
% % % Beta_aerosol_353 = interp1(altitude, beta_aero_353, zmsis);
% % % for i = 1:length(Beta_aerosol_308)
% % %     if Beta_aerosol_308(:,i) <0
% % %        Beta_aerosol_308(:,i) = 0;
% % %     end
% % % end
% % % 
% % % for i = 1:length(Beta_aerosol_308)
% % %     if ~isempty(find(isnan(Beta_aerosol_308(:,i)) == 1));
% % %         Beta_aerosol_308(:,i) = 0;
% % %     end
% % % end
% % % 
% % % for i = 1:length(Beta_aerosol_353)
% % %     if Beta_aerosol_353(:,i) <0
% % %        Beta_aerosol_353(:,i) = 0;
% % %     end
% % % end
% % % 
% % % for i = 1:length(Beta_aerosol_353)
% % %     if ~isempty(find(isnan(Beta_aerosol_353(:,i)) == 1));
% % %         Beta_aerosol_353(:,i) = 0;
% % %     end
% % % end
% % % gamma_aerosol_308 = exp(-2.*cumtrapz((Beta_aerosol_308(2:end))*LR).*(zmsis(2)-zmsis(1)));
% % % gamma_aerosol_353 = exp(-2.*cumtrapz((Beta_aerosol_353(2:end))*LR).*(zmsis(2)-zmsis(1)));
% % % 
% % % altitude = altitude';
% % % % z = zmsis(zmsis<=20000);
% % % 
% % % for i= 1:length(gamma_aerosol_308)
% % % if ~isempty(find(isnan(gamma_aerosol_308(i))) == 1);
% % %     gamma_aerosol_308(i) = gamma_aerosol_308(i-1);
% % % end
% % % end
% % % 
% % % for i= 1:length(gamma_aerosol_353)
% % % if ~isempty(find(isnan(gamma_aerosol_353(i))) == 1);
% % %     gamma_aerosol_353(i) = gamma_aerosol_353(i-1);
% % % end
% % % end
% % % gamma_aerosol_353 = [gamma_aerosol_353, gamma_aerosol_353(end)];
% % % gamma_aerosol_308 = [gamma_aerosol_308, gamma_aerosol_308(end)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta_aerosol_353 = 1e-8*exp(-(zmsis - 10000)/15000); % This comes from Frans paper. *3 * (10e-8) for aerosol free night
%and 3e-7 for dusty nights like 2006 data
% here is the vector of it
% beta_aerosol_353 = 1e-8*ones(1, length(zmsis));
LR = 50; % just a rough guess 
%% calculate the transmision for the aerosols at 308 nm.
gamma_aerosol_353 = exp(-2.*cumtrapz((beta_aerosol_353)*LR).*(zmsis(2)-zmsis(1)));
angestrom = 1;
beta_aerosol_308 = beta_aerosol_353 * (lambda_308/lambda_353)^(-angestrom);
gamma_aerosol_308 = exp(-2.*cumtrapz((beta_aerosol_308)*LR).*(zmsis(2)-zmsis(1)));
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = find(zmsis >=80000);
D1 = min(D);
F = find(zmsis <= 150000);
F1 = max(F);
K = find(zmsis >=80000);
K1 = min(K);
R = find(zmsis <= 150000);
R1 = max(R);
%Background consideration for both high and low gains at the two channels
bg_med_308 = mean(real_med_308(D1:F1));
var_bg_med_308 = std((real_med_308(D1:F1))).^2;
bg_med_353 = mean(real_med_353(K1:R1));
% bg_med_353 = real_med_353(end);
var_bg_med_353 = std((real_med_353(K1:R1))).^2;
bg_low_308 = mean(real_low_308(D1:F1));
var_bg_low_308 = std((real_low_308(D1:F1))).^2;
bg_low_353 = mean(real_low_353(D1:F1));
var_bg_low_353 = std((real_low_353(D1:F1))).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Calculating the conctant of the lidar for all channels and their gains
M1_308 = ((beta_mol_308+beta_aerosol_308).*gamma_air_308.*gamma_no3_308.*gamma_aerosol_308)./((zmsis).^2) ;
M1_353 = ((beta_mol_353+beta_aerosol_353).*gamma_air_353.*gamma_no3_353.*gamma_aerosol_353)./((zmsis).^2) ;
A = find(zmsis>=22000);
B = min(A);
A1 = find(zmsis>=15000);
B1 = min(A1);
C_low_308 = (real_low_308(B)./M1_308(B));
C_med_308 = (real_med_308(B)./M1_308(B));
C_low_353 = (real_low_353(B)./M1_353(B));
C_med_353 = (real_med_353(B)./M1_353(B));
Q.C_low_308 = C_low_308;
Q.C_med_308 = C_med_308;
Q.C_low_353 = C_low_353;
Q.C_med_353 = C_med_353;
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MM1_med_308 = C_med_308 *M1_308;
MM1_med_308 = MM1_med_308 + bg_med_308;
MM1_low_308 = C_low_308 *M1_308;
MM1_low_308 = MM1_low_308 + bg_low_308;
MM1_med_353 = C_med_353 *M1_353;
MM1_med_353 = MM1_med_353 + bg_med_353;
MM1_low_353 = C_low_353 *M1_353;
MM1_low_353 = MM1_low_353 + bg_low_353;
Q.MM1_med_308 = MM1_med_308;
Q.MM1_low_308 = MM1_low_308;
Q.MM1_med_353 = MM1_med_353;
Q.MM1_low_353 = MM1_low_353;
zmsis1 = zmsis;
Q.zmsis1 = zmsis1;
% semilogx(MM1_low_308, zmsis1); hold on;
%%Adding Noise to the System
% % for i = 1: length(MM1_med_308);
% %     [mdeviate_med_308(i)] = SLS_PoissonDeviate(MM1_med_308(i),1);
% %     
% %     
% % end
% % 
% % for i = 1: length(MM1_low_308);
% %     [mdeviate_low_308(i)] = SLS_PoissonDeviate(MM1_low_308(i),1);
% %     
% %     
% % end
% % 
% % 
% % for i = 1: length(MM1_med_353);
% %     [mdeviate_med_353(i)] = SLS_PoissonDeviate(MM1_med_353(i),1);
% %     
% %     
% % end
% % 
% % for i = 1: length(MM1_low_353);
% %     [mdeviate_low_353(i)] = SLS_PoissonDeviate(MM1_low_353(i),1);
% %     
% %     
% % end
 
counts_med_308 = MM1_med_308;
counts_low_308 = MM1_low_308;
counts_med_353 = MM1_med_353;
counts_low_353 = MM1_low_353;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
bot_ret = 9600;%11400;
Q.bot_ret = bot_ret;
top_ret = 98500;%89000; 
Q.top_ret = top_ret;
top_data = 98500; %89000;         
Q.top_data = top_data;
top_data_low_308 = 67000; %67400; 
Q.top_data_low_308 = top_data_low_308;
bot_data_low =9600; %11400; 
Q.bot_data_low = bot_data_low;
bot_data_high = 17100; %21300
Q.bot_data_high = bot_data_high;
% zret = bot_ret:350:18500;

zret = bot_ret:500:14000;
zret1 = 14000:700:21200;
% zret5 = 18000:700:21000;
% zret2 = 19500:900:25200;
zret3 = 21200:700:35300;
% zret1 = 25200:900:top_ret;
zret4 = 35300:900:top_ret;
%zret2 = 18200:900:top_ret;

% % zret4 = 23200:15000:25200;
% % zret5 = 25200:25000:27200;
% % zret6 = 27200:25000:29200;
% % zret7 = 29200: 2500: top_ret;
% % zret2 = 22200:10000:23200;
% % zret4 = 23200:15000:25200;
% % zret5 = 25200:25000:27200;
% % zret6 = 27200:25000:29200;
% % zret7 = 29200: 2500: top_ret;
% zret6 = 27200:2200:29200;
% zret7 = 29200:2300:32200;
% zret8 = 32200:2400:35200;
% zret9 = 35200:2500:top_ret;
% zret1 = 25200:900:top_ret;
% zret4 = 35300:800:top_ret;
% zret4 = 42000:4000:top_ret;
zret = [zret,zret1, zret3, zret4];
%%%%%%%%%%%%%%%%%%%%%%%%%
 
% O3_dens = interp1(zmsis, O3_dens, zret);
 
 
Q.O3_dens = O3_dens;
Q.zmsis = zmsis;
Q.Tmsis = Tmsis;                                                             
Q.sigma_Rayleigh_308 = sigma_Rayleigh_308;
Q.sigma_Rayleigh_353 = sigma_Rayleigh_353;
 
Q.zret = zret;
% zmsis = [7000:300:49000];
 
A = find(zmsis>=bot_data_low, 1 );
B = max(find(zmsis<=top_data));
C = min(find(zmsis>=bot_data_high));
D = max(find(zmsis<=top_data_low_308));
Q.A = A;
Q.B = B;
Q.C = C;
Q.D = D;
Q.Pray = Pray;
Q.lambda_308 = lambda_308;
Q.lambda_353 = lambda_353;
Q.ro_n = ro_n;
Q.N_s = N_s;
Q.n_s_308 = n_s_308; 
Q.n_s_353 = n_s_353; 
Q.sigo3_308 = sigo3_308;
sigo3_308_low = sigo3_308(A:B);
Q.sigo3_308_low = sigo3_308_low;
sigo3_308_lowj = sigo3_308j(A:B);
Q.sigo3_308_lowj = sigo3_308_lowj;
sigo3_b_parameter = sigo3_308(A:D);
Q.sigo3_b_parameter = sigo3_b_parameter;
T2_b_parameter = Tmsis(A:D);
Q.T2_b_parameter = T2_b_parameter;
sigo3_b_parameterj = sigo3_308j(A:D);
Q.sigo3_b_parameterj = sigo3_b_parameterj;
sigo3_308_med = sigo3_308(C:B);
Q.sigo3_308_med = sigo3_308_med;
sigo3_308_medj = sigo3_308(C:B);
Q.sigo3_308_medj = sigo3_308_medj;
T1_b_parameter = Tmsis(C:B);
Q.T1_b_parameter = T1_b_parameter;
Q.sigo3_353 = sigo3_353;
sigo3_353_low = sigo3_353(A:B);
Q.sigo3_353_low = sigo3_353_low;
sigo3_353_med = sigo3_353(C:B);
Q.sigo3_353_med = sigo3_353_med;
gamma_aerosol_308_low = gamma_aerosol_308(A:B);
Q.gamma_aerosol_308_low = gamma_aerosol_308_low;
gamma_aerosol_353_low = gamma_aerosol_353(A:B);
Q.gamma_aerosol_353_low = gamma_aerosol_353_low;
Beta_aerosol_353= beta_aerosol_353(A:B);
Q.Beta_aerosol_353 = Beta_aerosol_353;
Beta_aerosol_353_med= beta_aerosol_353(C:B);
Q.Beta_aerosol_353_med = Beta_aerosol_353_med;
Beta_aerosol_308_med= beta_aerosol_308(C:B);
Q.Beta_aerosol_308_med = Beta_aerosol_308_med;
Beta_aerosol_308 = beta_aerosol_308(A:B);
Q.Beta_aerosol_308 = Beta_aerosol_308;
Beta_aerosol_308_low = beta_aerosol_308(A:D);
Q.Beta_aerosol_308_low = Beta_aerosol_308_low;
delta_data = zmsis(2)-zmsis(1);
Q.delta_data = delta_data;
zmsis = bot_data_low:delta_data: top_data;
Q.zmsis = zmsis;
 
% % % 
% % % k = find(Q.zmsis<=16000);
% % % h = k(end);
% % % delta_ret = 300;
% % % Z1 =bot_ret:delta_ret:Q.zmsis(h);
% % % Z2= Q.zmsis(h):delta_ret*3:top_ret;  
% % % zret = [Z1 Z2];
% % % Q.zret = zret;
zmsis2 = bot_data_high:delta_data:top_data;
Q.zmsis2 = zmsis2;
zmsis3 = bot_data_low:delta_data:top_data_low_308;
Q.zmsis3 = zmsis3;
z0 = 755:150:bot_ret; l =length(z0);
no30  = Q.O3_dens(1:l);
sigo30_308 = Q.sigo3_308(1:l);
odnorm_308 = trapz(z0, no30.*sigo30_308);
Q.odnorm_308 = odnorm_308;
sigo30_353 = Q.sigo3_353(1:l);
odnorm_353 = trapz(z0, no30.*sigo30_353);
Q.odnorm_353 = odnorm_353;
nair0  = nmsis(1:l);
odnorm_nair_308 = trapz(z0, nair0.*sigma_Rayleigh_308);
Q.odnorm_nair_308= odnorm_nair_308;
odnorm_nair_353 = trapz(z0, nair0.*sigma_Rayleigh_353);
Q.odnorm_nair_353= odnorm_nair_353;
 
counts_med_308 = counts_med_308(C:B);
counts_low_308 = counts_low_308(A:D);
counts_med_353 = counts_med_353(C:B);
counts_low_353 = counts_low_353(A:B);
 
counts = [counts_med_308, counts_low_308, counts_med_353, counts_low_353];
for i = 1:length(counts);
    if counts(i)==0;
        counts(i) = 1;
    end
end
y = counts';
Q.y = y;
Sy = diag(y);
Q.Sy = Sy;
Q.counts_med_308 = counts_med_308;
Q.counts_low_308 = counts_low_308;
Q.counts_med_353 = counts_med_353;
Q.counts_low_353 = counts_low_353;
nM = length(counts_med_308);
nL = length(counts_low_308);
Q.nM = nM;
Q.nL = nL;
Q.bg_med_308 = bg_med_308;
Q.var_bg_med_308 = var_bg_med_308;
Q.bg_low_308 = bg_low_308;
Q.var_bg_low_308 = var_bg_low_308;
Q.bg_med_353 = bg_med_353;
Q.var_bg_med_353 = var_bg_med_353;
Q.bg_low_353 = bg_low_353;
Q.var_bg_low_353 = var_bg_low_353;
 
 
real_med_308 =  real_med_308(C:B);
real_low_308 =  real_low_308(A:D);
real_med_353 =  real_med_353(C:B);
real_low_353 =  real_low_353(A:B);
Q.real_med_308 = real_med_308;
Q.real_low_308 = real_low_308;
Q.real_med_353 = real_med_353;
Q.real_low_353 = real_low_353;
real_data = [real_med_308, real_low_308, real_med_353, real_low_353];
for i = 1:length(real_data);
    if real_data(i)==0;
        real_data(i) = 1;
    end
end
y_real = real_data';
Q.y_real = y_real;
%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
 
% % % % Var_med_308 = zeros(1, length(real_med_308));
% % % % Var_low_308 = zeros(1, length(real_low_308));
% % % % Var_med_353 = zeros(1, length(real_med_353));
% % % % Var_low_353 = zeros(1, length(real_low_353));
% % % % 
% % % %  
% % % % [Var_med_308,goH308] =bobpoissontest2(real_med_308,zmsis2);
% % % % [Var_med_353,goH353] =bobpoissontest2(real_med_353,zmsis2);
% % % % [Var_low_308,goL308] =bobpoissontest2(real_low_308,zmsis3);
% % % % [Var_low_353,goL353] =bobpoissontest2(real_low_353,zmsis);
% % % % 
% % % % r1 = ones(1,goH308-1).* Var_med_308(1);
% % % % r2 = ones(1,goH308-1).* Var_med_308(end);
% % % % 
% % % % R1 = ones(1,goH353-1).* Var_med_353(1);
% % % % R2 = ones(1,goH353-1).* Var_med_353(end);
% % % % 
% % % % r3 = ones(1,goL308-1).* Var_low_308(1);
% % % % r4 = ones(1,goL308-1).* Var_low_308(end);
% % % % 
% % % % 
% % % % R3 = ones(1,goL353-1).* Var_low_353(1);
% % % % R4 = ones(1,goL353-1).* Var_low_353(end);
% % % % 
% % % % Var_med_308 = [r1, Var_med_308 ,r2];
% % % % Var_med_353 = [R1, Var_med_353 ,R2];
% % % % Var_low_308 = [r3, Var_low_308, r4];
% % % % Var_low_353 = [R3, Var_low_353, R4];
% % % % 
% % % % for i = 1: length(zmsis3)
% % % % 
% % % %     if zmsis3(i) <= 13000
% % % % 
% % % %         YY_L_308(i) = Var_low_308(i);
% % % % 
% % % %     else
% % % % 
% % % %         YY_L_308(i) = real_low_308(i);
% % % % 
% % % %     end
% % % % 
% % % % end
% % % % 
% % % % 
% % % % for i = 1: length(zmsis)
% % % % 
% % % %     if zmsis(i) <= 13000
% % % % 
% % % %         YY_L_353(i) = Var_low_353(i);
% % % % 
% % % %     else
% % % % 
% % % %         YY_L_353(i) = real_low_353(i);
% % % % 
% % % %     end
% % % % 
% % % % end
% % % % 
% % % % for i = 1: length(zmsis2)
% % % % 
% % % %     if zmsis2(i) <= 21000
% % % % 
% % % %         YY_M_308(i) = Var_med_308(i);
% % % %         
% % % %     else
% % % % 
% % % %         YY_M_308(i) = real_med_308(i);
% % % % 
% % % %     end
% % % % 
% % % % end
% % % % 
% % % % 
% % % % for i = 1: length(zmsis2)
% % % % 
% % % %     if zmsis2(i) <= 21000
% % % % 
% % % %         YY_M_353(i) = Var_med_353(i);
% % % %         
% % % %     else
% % % % 
% % % %         YY_M_353(i) = real_med_353(i);
% % % % 
% % % %     end
% % % % 
% % % % end
 
 
% % % % 
% % % % Yvar =[YY_M_308 YY_L_308 YY_M_353 YY_L_353];
% % % % 
% % % % Sy_real = diag(Yvar);
% % % % Q.Sy_real = Sy_real;
 
lzA = length(Q.zmsis3);
lzB = length(Q.zmsis);
lzC = length(Q.zmsis2);
go = 7; % 6; 12; %24
    stop = go-1;
    j = 0;
    for i = go:lzA-stop
        j = j + 1;
        [pp,spp,ppregress] = fitlinenp(Q.zmsis3(i-stop:i+stop),Q.real_low_308(i-stop:i+stop));
        tmp = pp(1).*Q.zmsis3(i-stop:i+stop) + pp(2);
        varl308(i) = (std(Q.real_low_308(i-stop:i+stop) - tmp)).^2;
        
    end
    for i = go:lzB-stop
        j = j + 1;
        [pp,spp,ppregress] = fitlinenp(Q.zmsis(i-stop:i+stop),Q.real_low_353(i-stop:i+stop));
        tmp = pp(1).*Q.zmsis(i-stop:i+stop) + pp(2);
        varl353(i) = (std(Q.real_low_353(i-stop:i+stop) - tmp)).^2;
        
    end
    
    for i = go:lzC-stop
        j = j + 1;
        [pp,spp,ppregress] = fitlinenp(Q.zmsis2(i-stop:i+stop),Q.real_med_353(i-stop:i+stop));
        tmp = pp(1).*Q.zmsis2(i-stop:i+stop) + pp(2);
        varm353(i) = (std(Q.real_med_353(i-stop:i+stop) - tmp)).^2;
        
    end
  go1 = 7;  
      for i = go1:lzC-stop
        j = j + 1;
        [pp,spp,ppregress] = fitlinenp(Q.zmsis2(i-stop:i+stop),Q.real_med_308(i-stop:i+stop));
        tmp = pp(1).*Q.zmsis2(i-stop:i+stop) + pp(2);
        varm308(i) = (std(Q.real_med_308(i-stop:i+stop) - tmp)).^2;
        
    end
%%%%%%%%%%%%%%%%%%%%%%%
varlow308= zeros(size(Q.real_low_308));
varlow308(go:lzA-stop) = varl308(go:lzA-stop);
    varlow308(1:go-1) = varl308(go);
    varlow308(lzA-stop+1:end) = varlow308(lzA-stop);
 
    varlow353 =  zeros(size(Q.real_low_353));
varlow353(go:lzA-stop) = varl353(go:lzA-stop);
    varlow353(1:go-1) = varl353(go);
    varlow353(lzA-stop+1:end) = varlow353(lzA-stop);
 
 
    varmed353 =  zeros(size(Q.real_med_353));
varmed353(go1:lzC-stop) = varm353(go1:lzC-stop);
    varmed353(1:go1-1) = varm353(go1);
    varmed353(lzC-stop+1:end) = varmed353(lzC-stop);  
    
    
   varmed308 =  zeros(size(Q.real_med_308));
varmed308(go1:lzC-stop) = varm308(go1:lzC-stop);
    varmed308(1:go1-1) = varm308(go1);
    varmed308(lzC-stop+1:end) = varmed308(lzC-stop);     
    
var_low_308 = [varlow308(1:10), real_low_308(11:end)] ;   
var_low_353 = [varlow353(1:10), real_low_353(11:end)] ; 
var_med_353 = [varmed353(1:30), real_med_353(31:end)] ; 
var_med_308 = [varmed308(1:30), real_med_308(31:end)] ; 
var = [var_med_308, real_low_308, var_med_353 , real_low_353];
Sy_real = diag(var);
for i =1:length(real_data)
    if real_data(i) <15;
        Sy_real(i) = 15;
        real_data(i) = real_data(i) + rand(1)*sqrt(Sy_real(i));
    end
end
 
Sy_real = diag(real_data);      
% Q.Sy_real = Sy_real;
 
O3_dens = interp1(alt_data, O3_dens, zret);
 
% no3_p = 4.5e12*exp(-1.2e-3*(Q.zret/1000-24).^2)*1e6;
% Q.no3_p = 1*O3_dens;
% % O3_dens(1:9) = O3_dens(1:9)*1.6;
% % O3_dens(9:21) = O3_dens(9:21)*.8;
% % o3_trad = interp1(objectdata.O3*1000, objectdata.O3*1e6, Q.zret);
% % no3 = (O3_dens./(O3_dens(1)).*o3_trad); 
Q.no3_p = 1.0*O3_dens;
 
nmsis = nmsis(A:B);
nair_sonde = interp1(alts*1000, air_density, Q.zret);
n_air_p = interp1(Q.zmsis, nmsis, Q.zret);
% air = [nair_sonde(1:40), n_air_p(41:end)];
Q.n_air_p = n_air_p;
 
Q.tau1= 4.6;%7 
Q.tau2 = 4.6;
Q.tau3 = 4.6;
Q.tau4 = 4.6;  
 
 
% % % % Sx(n-14, n-14) = (.4*(Q.tau1)).^2; low
% % % % Sx(n-15, n-15) = (.4*(Q.tau2)).^2; med
% % % % Sx(n-16, n-16) = (.4*(Q.tau3)).^2; med
% % % % Sx(n-17, n-17) = (.4*(Q.tau4)).^2; low
Q.LR = LR;
Q.angestrom = angestrom;
Q.a1 = -.25;
Q.a0 = 1.6;
Q.a2 = -3.5;
Q.a3 = 1;
Q.a4 = 12;
Q.a5 = -0.002;
Q.date = date;
 
% Q.zmsis1 = zmsis1;
 
 


