function [dS_M_308, dS_L_308, dS_M_353, dS_L_353] = derivSHSN5(Q,x,j,forwardmodel_test7_308_353)
% derivative of forward model with respect to x

[SM308, SL308, SM353, SL353] = forwardmodel_test7_308_353(Q,x);

if ~isempty(find(isnan(x)) == 1)
    'after FM: Nans in retrieval vector derivSHSN2'
    stop
end

dn = 1.e-1 .* x(j); % 1e-5
xpert = x;
if x(j) == 0 % trap for tau's where tau(1) = 0
    dn = 1.e-1 .* x(j+1);
end
xpert(j) = x(j) + dn;

[SM308j, SL308j, SM353j, SL353j] = forwardmodel_test7_308_353(Q,xpert);
%dFdRat = (ratj - rat) ./ dn;
dS_M_308 = (SM308j - SM308) ./ dn;
dS_L_308 = (SL308j - SL308) ./ dn;
dS_M_353 = (SM353j - SM353) ./ dn;
dS_L_353 = (SL353j - SL353) ./ dn;
%dFdRat = dSHdx./SH - dSNdx./SN;

return