OHP_sonde_compare;
air = X.x(m+1:m+1+30);
zret = Q.zret(1:31);
air_p = Q.n_air_p(1:31);
figure;
plot(((air_p-air')./air_p)*100, zret)
title('difference btw the a priori and the OEM retrieval')
sonde = interp1(alts(2500:end)*1000, air_density(2500:end), zret);
figure;
plot(((air'-sonde)./sonde)*100, zret)
title('differene btw the sonde and the OEM retrieval for air density')
figure; 
plot(((air_p-sonde)./sonde)*100, zret)
title('difference btw a priori and teh sonde')
O3 = X.x(1:31)/1e6;
zret = Q.zret(1:31);
sondeO3 = interp1(alts(2000:end)*1000, conco3s(2000:end), zret);
figure;
plot(((O3'-sondeO3)./O3')*100, zret)
title('difference btw the sonde and OEM retrieval for ozone')
%%%%%%%%%% 
figure;
scale = sonde(25)./air(25);
air_scale = air*scale;
semilogx(air_scale, zret); hold on; semilogx(sonde, zret)
%%%%%
figure; 
avg = (sonde(1:end)'+air(1:end))/2;
plot(((sonde(1:end)'-air(1:end))./avg)*100, zret(1:end)/1000)

title('OEM and sonde difference')