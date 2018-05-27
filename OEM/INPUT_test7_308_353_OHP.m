% This is collects all the inputs requires for the oem.m code

function [O,Q,R,Sx,Se,y,x] = INPUT_test7_308_353_OHP


%set Q

%  [Q, y, Se] = makeQ_test7_308_353( );
%
[Q, y, Se] = makeQ_test7_308_353_OHP( );
% [Q, y, Se] = makeQ_test5_308_353_alpha_Sophie( );

%set O 
%%%%O = defOreal;
[O] = makeO_DIAL;

%set R, retrieval structure
R = [];
R.jq = {};
R.ji = {};
iter = 1;

x = [(Q.no3_p), (Q.n_air_p),Q.tau4,Q.tau3,Q.tau2,(Q.tau1),Q.bg_low_353,Q.bg_med_353,Q.bg_med_308,... 
    Q.a0,Q.a1,Q.a4,Q.C_low_353,(Q.C_med_353),(Q.C_low_308),...
    (Q.C_med_308), Q.odnorm_nair_353 ,...
    Q.odnorm_nair_308, Q.odnorm_353,Q.odnorm_308];

x = x';
m = length(Q.zret); %% length for the retrieval grid
n = (2*m)+18;
% Sx = zeros(n,n);
S_empty = zeros(m,m);
[So3]=TestTempCov_new(Q.zret,(Q.no3_p));
[Sair]=TestTempCov_n_air(Q.zret,(Q.n_air_p));

S1 = [So3; S_empty];
S2 = [S_empty; Sair];

Sa1 = [ S1, S2]; 
T1 = zeros(18, 2*m); 
Sa2 = [Sa1; T1];
T2 = zeros(n, 18);
Sa3 = [Sa2, T2];
Sx = Sa3;
Sx(n,n) = (.6*(Q.odnorm_308)).^2;
Sx(n-1,n-1) = (.6*(Q.odnorm_353)).^2;
Sx(n-2, n-2) = (.6*(Q.odnorm_nair_308)).^2;
Sx(n-3, n-3) = (.6*(Q.odnorm_nair_353)).^2;
Sx(n-4, n-4) = (0.6*(Q.C_med_308)).^2;
Sx(n-5, n-5) = (0.6*(Q.C_low_308)).^2;
Sx(n-6, n-6) = (0.6*(Q.C_med_353)).^2;
Sx(n-7, n-7) = (0.6*(Q.C_low_353)).^2;
Sx(n-8, n-8) = 0.2*(Q.a4).^2;
Sx(n-9, n-9) = 0.2*(Q.a1).^2;
Sx(n-10, n-10) = 0.2*(Q.a0).^2;
Sx(n-11, n-11) =  Q.var_bg_med_308;
Sx(n-12, n-12) = Q.var_bg_med_353;
Sx(n-13, n-13) = Q.var_bg_low_353;
Sx(n-14, n-14) = (.5*(Q.tau1)).^2;
Sx(n-15, n-15) = (.5*(Q.tau2)).^2;
Sx(n-16, n-16) = (.5*(Q.tau3)).^2;
Sx(n-17, n-17) = (.5*(Q.tau4)).^2;
% % Sx(n-17, n-17) = 0.2*(Q.var_bg_med_353)^2;
% % Sx(n-18, n-18) = 0.2*(Q.var_bg_low_353)^2;

% [R, yf, J] = makeJ_test7_308_353(Q, R, x, iter);
% [R, yf, J] = makeJ_test7_308_353_OHP(Q, R, x, iter);
end 