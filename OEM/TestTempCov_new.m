
function [Sa_T]=TestTempCov_new(Zj,x)
 
 m = length(x);
 n = m;
%% for Sophie's data, Length of 900 and uncertainty of .1*x is the best fit with h
%% het smoothed anf filtered data
%  lengthcT = 900; % =3000; % only need m of these
% Tfac =  0.40*x; % The diffrence btw two ozone for OHP .15 is the best
 % for Eureka should go higher to .3
% lc = ones(1,m);
lengthcT_low = 500;            %2000;   
lengthcT_med = 1400;            %5000;
lengthcT_high = 1400;             %7500;
% %  lc(Zj<=40000) = 900;
% %  lc(Zj>40000) = 150;
% % 
Tfac = ones(1, length(Zj));
Tfac(1,1:9) = .5*x(1:9);
Tfac(1,9:26) = .5*x(9:26);
% Tfac(1,17:26) = .3*x(17:26);
Tfac(26:end) = .2*x(26:end);
%%%%%%%

% Tfac = .5*x;
Tmodvar= Tfac;
% Tmodvar = (Tfac.*ones(size(x)));
%Tmodvar = (Tfac.*Sa).^2;
 vars2 = Tmodvar;
lc_low = lengthcT_low.*ones(1,9);
%h = m-25;
lc_med = lengthcT_med.*ones(1,17);
%v = m-h;
lc_high = lengthcT_high.*ones(1,m-(17+9));
lc = [lc_low , lc_med, lc_high];
 Sa_T =zeros(n,n);

for i = 1:m
    for j = 1:m
        
        sigprod = (vars2(i).*vars2(j));
        diffz = Zj(i) - Zj(j);
        sumlc = lc(i) + lc(j);
        shape(3) = (1-(1-exp(-1)).*2.*abs(diffz)./sumlc);
        
        if shape(3) < 0
            shape(3) = 0;
        end
        
        Sa_T(i,j) = sigprod.*shape(3);

    end
end

 Sa_T(n,n) = vars2(n);


        