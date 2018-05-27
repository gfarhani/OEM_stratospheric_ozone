
function [Sa_T]=TestTempCov_n_air(Zj,x)
 
 m = length(x);
 n = m;
lengthcT_low = 300;
lengthcT_med = 900;
lengthcT_high = 1500;
% %  lc(Zj<=40000) = 900;
% %  lc(Zj>40000) = 150;
% % 
Tfac = ones(1, length(Zj));
Tfac(1,1:11) = .3*x(1:11);
Tfac(1,11:17) = .3*x(11:17);
% Tfac(1,17:26) = .3*x(17:26);
Tfac(17:end) = .3*x(17:end);
% % Tfac = ones(1, length(Zj));
% % Tfac(1,1:20) = 0.2*x(1:20);
% % Tfac(1,21:40) = .3*x(21:40);
% % Tfac(41:end) = .55*x(41:end);
% % Tmodvar= Tfac;
Tmodvar= Tfac;
% Tmodvar = (Tfac.*ones(size(x)));
%Tmodvar = (Tfac.*Sa).^2;
 vars2 = Tmodvar;
lc_low = lengthcT_low.*ones(1,9);
h = m-9;
lc_med = lengthcT_med.*ones(1,h);
v = m-h;
lc_high = lengthcT_high.*ones(1,v);
lc = [lc_low , lc_med, lc_high];
 Sa_T =zeros(n,n);
% lc = lengthcT_low.*ones(1,m);
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