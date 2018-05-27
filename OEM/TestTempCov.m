
function [Sa_T]=TestTempCov(Zj,x)
 
 m = length(x);
 n = m;
%% for Sophie's data, Length of 900 and uncertainty of .1*x is the best fit with h
%% het smoothed anf filtered data
%  lengthcT = 900; % =3000; % only need m of these
% Tfac =  0.40*x; % The diffrence btw two ozone for OHP .15 is the best
 % for Eureka should go higher to .3
% lc = ones(1,m);
%% 
% lengthcT = 10;
% %  lc(Zj<=410000) = 900;
% %  lc(Zj>40000) = 150;X
% % 
lengthcT = 300;
% %  lc(Zj<=410000) = 900;
% %  lc(Zj>40000) = 150;X
% % 
Tfac = ones(1, length(Zj));
Tfac(1,1:11) = .5*x(1:11);
Tfac(1,12:25) = .5*x(12:25);
Tfac(1,25:end) = .2*x(25:end);
% Tfac = .5*x;
Tmodvar= Tfac;
% Tfac = .5*x;
% Tmodvar= Tfac;
% Tfac = .5*x;

% Tfac = .5*x;
% Tfac = .5*x;
% Tfac = .5*x;
Tmodvar= Tfac;
%Tmodvar = (Tfac.*Sa).^2;
 vars2 = Tmodvar;
lc = lengthcT.*ones(1,m);
 Sa_T =zeros(n,n);

for i = 1:m
    for j = 1:m
        
        sigprod = (vars2(i).*vars2(j));
        diffz = Zj(i) - Zj(j);
        sumlc = lc(i) + lc(j);
        shape(3) = (1-(1-exp(-1)).*2.*abs(diffz)./sumlc);
        shape(3) = exp(-(2.*diffz./sumlc).^2);
        if shape(3) < 0
            shape(3) = 0;
        end
        
        Sa_T(i,j) = sigprod.*shape(3);

    end
end

 Sa_T(n,n) = vars2(n);


        