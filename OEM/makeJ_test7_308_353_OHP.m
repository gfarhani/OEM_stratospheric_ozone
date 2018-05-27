

% Jacobian for T
function  [R, yf, J] = makeJ_test7_308_353_OHP(Q, R, x, iter)
% Zi is data vector
% Zj is ret vector
Zi = Q.zmsis;
Zj = Q.zret;
zmsis1 = Zi(Zi>=Q.bot_data_high);
zmsis2 = Zi(Zi>=Q.bot_data_low);

[SM_308, SL_308, SM_353, SL_353]=forwardmodel_test7_308_353(Q,x);
SM_308 = SM_308(Zi>=Q.bot_data_high);
SM_353 = SM_353(Zi>=Q.bot_data_high);
SL_308 = SL_308(Zi<=Q.top_data_low_308);
SL_353 = SL_353(Zi>=Q.bot_data_low);
m1 = length(SM_308); m2 = length(SL_308); 
SM_308 = SM_308'; SL_308 = SL_308'; SM_353 = SM_353'; SL_353 = SL_353'; 
yf = [SM_308; SL_308; SM_353; SL_353];
mdata = length(yf);

m = length(Q.zret);n = 2*m+18;
Kernel = zeros(mdata, n);

for j = 1:2*m
    [dS_M_308, dS_L_308, dS_M_353, dS_L_353] = derivSHSN5(Q,x,j,@forwardmodel_test7_308_353);
    dSM_308 = dS_M_308(Zi>=Q.bot_data_high);
    dSL_308 = dS_L_308(Zi<=Q.top_data_low_308);
    dSM_353 = dS_M_353(Zi>=Q.bot_data_high);
    dSL_353 = dS_L_353(Zi>=Q.bot_data_low);
    Kernel(1:m1,j) = dSM_308;
    Kernel(m1+1:m1+m2,j) = dSL_308;
    Kernel(m1+m2+1: m2+(2*m1), j) = dSM_353;
    Kernel(m2+(2*m1)+1 : mdata, j) = dSL_353;
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Ozone %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % figure;
% % % % % % % % % % set(gca,'fontsize',9)
% % % subplot(1,2,1)
% % % plot(Kernel(1:m1,2:m), zmsis1/1000)
% % % title ('O3 jacobian 308 nm channel')
% % % ylim([10,50])
% % % subplot(1,2,2)
% % % plot(Kernel(m1+1:m1+m2,2:m), zmsis2/1000)
% % % ylim([10, 50])
% % % title ('O3 jacobian 308 nm channel')
% % % % % % % % % % % % % % % 
% % % % % % % % % % figure
% % % % % % % % % % set(gca,'fontsize',9)
% % % % % % % % % % subplot(1,2,1)
% % % % % % % % % % plot(Kernel(m1+m2+1:m2+(2*m1),2:m), zmsis1/1000)
% % % % % % % % % % title ('O3 jacobian 353 nm channel')
% % % % % % % % % % ylim([10,50])
% % % % % % % % % % subplot(1,2,2)
% % % % % % % % % % plot(Kernel(m2+(2*m1)+1:mdata,2:m), zmsis2/1000)
% % % % % % % % % % ylim([10, 50])
% % % % % % % % % % title ('O3 jacobian 353 nm channel')
% % % % % % % % % % % % % % % % % 
% % % % % % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % %%%%%%%%%%%%%%%%%%%% air density %%%%%%%%%%%%%%%%%%
% % % % % % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % figure;
% % % % % % % % % % set(gca,'fontsize',9)
% % % % % % % % % % subplot(1,2,1)
% % % % % % % % % % plot(Kernel(1:m1,m+1:2*m), zmsis1/1000)
% % % % % % % % % % title ('nair jacobian 308 nm channel')
% % % % % % % % % % ylim([10,50])
% % % % % % % % % % subplot(1,2,2)
% % % % % % % % % % plot(Kernel(m1+1:m1+m2,m+1:2*m), zmsis2/1000)
% % % % % % % % % % ylim([10, 50])
% % % % % % % % % % title ('nair jacobian 308 nm channel')
% % % % % % % % % % % % % % % 
% % % % % % % % % % figure
% % % % % % % % % % set(gca,'fontsize',9)
% % % % % % % % % % subplot(1,2,1)
% % % % % % % % % % plot(Kernel(m1+m2+1:m2+(2*m1),m+1:2*m), zmsis1/1000)
% % % % % % % % % % title ('nair jacobian 353 nm channel')
% % % % % % % % % % ylim([10,50])
% % % % % % % % % % subplot(1,2,2)
% % % % % % % % % % plot(Kernel(m2+(2*m1)+1:mdata,m+1:2*m), zmsis2/1000)
% % % % % % % % % % ylim([10, 50])
% % % % % % % % % % title ('nair jacobian 353 nm channel')
% % % % % % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % %%%%%%%%%%%%%%%%%%%%%  beta_aerosol 353 and 308 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % figure;
% % % % % % % % % % set(gca,'fontsize',9)
% % % % % % % % % % subplot(1,2,1)
% % % % % % % % % % plot(Kernel(1:m1,(2*m)+1:3*m), zmsis1/1000)
% % % % % % % % % % title ('beta_aerosol_308 jacobian 308 nm channel')
% % % % % % % % % % ylim([10,50])
% % % % % % % % % % subplot(1,2,2)
% % % % % % % % % % plot(Kernel(m1+1:m1+m2,(2*m)+1:3*m), zmsis2/1000)
% % % % % % % % % % ylim([10, 50])
% % % % % % % % % % title ('beta_aerosol_308 jacobian 308 nm channel')
% % % % % % % % % % % % % % % 
% % % % % % % % % % figure
% % % % % % % % % % set(gca,'fontsize',9)
% % % % % % % % % % subplot(1,2,1)
% % % % % % % % % % plot(Kernel(m1+m2+1:m2+(2*m1),(2*m)+1:3*m), zmsis1/1000)
% % % % % % % % % % title ('beta_aerosol_353 jacobian 353 nm channel')
% % % % % % % % % % ylim([10,50])
% % % % % % % % % % subplot(1,2,2)
% % % % % % % % % % plot(Kernel(m2+(2*m1)+1:mdata,(2*m)+1:3*m), zmsis2/1000)
% % % % % % % % % % ylim([10, 50])
% % % % % % % % % % title ('beta_aerosol_353 jacobian 353 nm channel')
% % % % % % % 
% % % % % % % 
% % % % % % % 
[dSM_308, dSL_308, dSM_353, dSL_353] = derivSHSN5(Q,x,n,@forwardmodel_test7_308_353);
Kernel(1:m1, n) = dSM_308(Zi>=Q.bot_data_high);
Kernel(m1+1:m1+m2, n) = dSL_308(Zi<=Q.top_data_low_308);
[dSM_308, dSL_308, dSM_353, dSL_353] = derivSHSN5(Q,x,n-1,@forwardmodel_test7_308_353);
Kernel(m1+m2+1:(2*m1)+m2, n-1) = dSM_353(Zi>=Q.bot_data_high);
Kernel((2*m1)+m2+1:mdata, n-1) = dSL_353;
[dSM_308, dSL_308, dSM_353, dSL_353] = derivSHSN5(Q,x,n-2,@forwardmodel_test7_308_353);
Kernel(1:m1, n-2) = dSM_308(Zi>=Q.bot_data_high);
Kernel(m1+1:m1+m2, n-2) = dSL_308(Zi<=Q.top_data_low_308);
[dSM_308, dSL_308, dSM_353, dSL_353] = derivSHSN5(Q,x,n-3,@forwardmodel_test7_308_353);
Kernel(m1+m2+1:(2*m1)+m2, n-3) = dSM_353(Zi>=Q.bot_data_high);
Kernel((2*m1)+m2+1:mdata, n-3) = dSL_353;
[dSM_308, dSL_308, dSM_353, dSL_353] = derivSHSN5(Q,x,n-4,@forwardmodel_test7_308_353);
Kernel(1:m1, n-4) = dSM_308(Zi>=Q.bot_data_high);
[dSM_308, dSL_308, dSM_353, dSL_353] = derivSHSN5(Q,x,n-5,@forwardmodel_test7_308_353);
Kernel(m1+1:m1+m2, n-5) = dSL_308(Zi<=Q.top_data_low_308);
[dSM_308, dSL_308, dSM_353, dSL_353] = derivSHSN5(Q,x,n-6,@forwardmodel_test7_308_353);
Kernel(m1+m2+1:(2*m1)+m2, n-6) = dSM_353(Zi>=Q.bot_data_high);
[dSM_308, dSL_308, dSM_353, dSL_353] = derivSHSN5(Q,x,n-7,@forwardmodel_test7_308_353);
Kernel((2*m1)+m2+1:mdata, n-3) = dSL_353;
[dSM_308, dSL_308, dSM_353, dSL_353] = derivSHSN5(Q,x,n-8,@forwardmodel_test7_308_353);
Kernel(m1+1:m1+m2, n-8) = dSL_308(Zi<=Q.top_data_low_308);
[dSM_308, dSL_308, dSM_353, dSL_353] = derivSHSN5(Q,x,n-9,@forwardmodel_test7_308_353);
Kernel(m1+1:m1+m2, n-9) = dSL_308(Zi<=Q.top_data_low_308);
[dSM_308, dSL_308, dSM_353, dSL_353] = derivSHSN5(Q,x,n-10,@forwardmodel_test7_308_353);
Kernel(m1+1:m1+m2, n-10) = dSL_308(Zi<=Q.top_data_low_308);
[dSM_308, dSL_308, dSM_353, dSL_353] = derivSHSN5(Q,x,n-11,@forwardmodel_test7_308_353);
Kernel(1:m1, n-11) = dSM_308(Zi>=Q.bot_data_high);
[dSM_308, dSL_308, dSM_353, dSL_353] = derivSHSN5(Q,x,n-12,@forwardmodel_test7_308_353);
Kernel(m1+m2+1:(2*m1)+m2, n-12) = dSM_353(Zi>=Q.bot_data_high);
[dSM_308, dSL_308, dSM_353, dSL_353] = derivSHSN5(Q,x,n-13,@forwardmodel_test7_308_353);
Kernel(m1+1:m1+m2, n-13) = dSL_308(Zi<=Q.top_data_low_308);
[dSM_308, dSL_308, dSM_353, dSL_353] = derivSHSN5(Q,x,n-14,@forwardmodel_test7_308_353);
Kernel((2*m1)+m2+1:mdata, n-14) = dSL_353;
%%%%%%
[dSM_308, dSL_308, dSM_353, dSL_353] = derivSHSN5(Q,x,n-15,@forwardmodel_test7_308_353);
Kernel(1:m1, n-15) = dSM_308(Zi>=Q.bot_data_high);
%%%%%%%%%
[dSM_308, dSL_308, dSM_353, dSL_353] = derivSHSN5(Q,x,n-16,@forwardmodel_test7_308_353);
Kernel(m1+m2+1:(2*m1)+m2, n-16) = dSM_353(Zi>=Q.bot_data_high);
%%%%%%%%%%%%%
[dSM_308, dSL_308, dSM_353, dSL_353] = derivSHSN5(Q,x,n-17,@forwardmodel_test7_308_353);
Kernel(m1+1:m1+m2, n-17) = dSL_308(Zi<=Q.top_data_low_308);


J = Kernel;
end