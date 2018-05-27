clear all;
close all;
warning off;
[O,Q,R,Sx,Se,y,x] = INPUT_test7_308_353_OHP;
m = length(Q.zret);
I = ones(size(Sx));

[X,R] = oem_try2(O,Q,R,@makeJ_test7_308_353_OHP,Sx,Se,[],[],x,y);

subplot(1,2,1)
uFakvec = 2;
degF1 = trace(X.A(1:m,1:m));
degF2 = trace(X.A(m+1:2*m,m+1:2*m));
finiDoF = round(min(degF1,degF2));

X.A = X.A';
hk = plot(X.A(3:3:m, 3:m), Q.zret(3:end)/1000, 'LineWidth', 1.5); hold on; 
unit = ones(size(Q.zret(3:end)));
%unit = unit';
response = X.A(3:m,3:m)*unit';
plot(response,Q.zret(3:m)./1000,'r', 'LineWidth', 1.5);

fak = find(diag(X.A(2:m,2:m)) >= 0.8);
fak2 = find(response >= 0.8);
if isempty(fak2)
    fak2 = 1;
end
fakvec = [fak(end); fak2(end); finiDoF]';
if uFakvec == 0
    fini = finiDoF;
else
    fini = fakvec(uFakvec);
end
set(hk,'LineWidth',1.5)
xlabel('Ozone AK')
ylabel('Altitude (km)')
set(gca,'fontsize',14)
%axis([-0.1 1.1 Q.zRET(1)./1000 Q.zRET(end)./1000])
xlim([-0.01 1.01])
ylim([12.6 70])
pltx = get(gca,'XLim');
plot(pltx,[Q.zret(fini) Q.zret(fini)]./1000,'k--')
% % xlim([-1 1])
subplot(1,2,2)

hk = plot(X.A(m+1:3:2*m-1-10, m+1:2*m-10), Q.zret(1:end-10)/1000, 'LineWidth', 1.5); hold on;
unit = ones(size(Q.zret(1:end-10)));
% unit = unit;
response = X.A(m+1:2*m-10, m+1:2*m-10)*unit';
plot(response,Q.zret(1:m-10)./1000,'r');
set(hk,'LineWidth',1.5)
xlabel('Air Density AK')
ylabel('Altitude (km)')
%axis([-0.1 1.1 Q.zRET(1)./1000 Q.zRET(end)./1000])
xlim([-0.01 1.02])
ylim([12.6 90])
pltx = get(gca,'XLim');
set(gca,'fontsize',14)
plot(pltx,[Q.zret(fini) Q.zret(fini)]./1000,'k--')
% % ylim([1 100])
% % xlim([-1 1])
% title('Averaging Kernel for Ozone in the 308 nm channel');
% figure; plot((X.x(1:end-8)), Q.zret/1000, ':o'); hold on; plot(Q.no3_p, Q.zret/1000, ':*');
% title('Ozone retrieval is compared with the a_priori profile date 2006,02,24');
% [objectdata, numberdensity, air_density, sondealt] = readingsonde('sonde20050301.txt');
% numberdensity =  sgolayfilt(numberdensity, 1, 11);
% [objectdata, numberdensity, sondealt] = readingozonesonde2016('scan1.txt');
%  numberdensity = smooth(numberdensity, 10);
figure; plot((X.x(2:m))./1e6, Q.zret(2:end)/1000, 'LineWidth', 1.5);hold all; plot(Q.no3_p./1e6, Q.zret/1000, 'LineWidth', 1.5);  
% hold on; plot(numberdensity, objectdata.alt/1000, 'black', 'LineWidth', 1.5);
% legend('retrieval', 'a priori', 'Ozonesonde');
legend('retrieval', 'a priori')
ylim([7 60])                       
xlim([0, 10e12]);
xlabel ('Ozone (mol/cm^3)'); ylabel('Altitude (km)');
pltx = get(gca,'XLim');
plot(pltx,[Q.zret(fini) Q.zret(fini)]./1000,'k--')
% title('Ozone retrieval is campoared with the a_priori profile and the ozone sonde profile')
% % figure; plot((X.x(2:m))./1e6, Q.zret(2:end)/1000, 'LineWidth', 1.5); 
% % % hold on; plot(Q.no3_p./1e6, Q.zret/1000, 'LineWidth', 1.5);
% % hold on;
% % pltx = get(gca,'XLim');
% % plot(pltx,[Q.zret(fini) Q.zret(fini)]./1000,'k--')
% % filename = 'sro1702071803.c1l';
% % [ objectdata ] = OHP_O3(filename);
% % plot(objectdata.O3, objectdata.alt)
% % legend('OEM', 'Traditinal') 
% % ylim([1 70])
% % xlim([0, 10e12]);
% % xlabel ('Ozone (mol/cm^3)'); ylabel('Altitude (km)');
figure; semilogx((X.x(m+1:2*m)), Q.zret(1:end)/1000, 'LineWidth', 1.5); hold all; semilogx(Q.n_air_p, Q.zret/1000, 'LineWidth', 1.5);
pltx = get(gca,'XLim');
% plot(pltx,[Q.zret(fini) Q.zret(fini)]./1000,'k--')
ylim([1 70])
legend('OEM', 'a priori')
xlabel ('Air density (#/m^3)'); ylabel('Altitude (km)');
%%%%%%%%%%%%%%%%%%% residulas%%%%%%%%%%%%%%%%%%%%

figure;
set(gca,'fontsize',18)
subplot(2,2,1)
hold on;
yf = X.yf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(((y(1:Q.nM)-yf(1:Q.nM))./yf(1:Q.nM)).*100,Q.zmsis(Q.zmsis>=Q.bot_data_high)./1000,'b');
hold on;
plot(-sqrt(y(1:Q.nM))./yf(1:Q.nM).*100,Q.zmsis((Q.zmsis>=Q.bot_data_high))./1000,'r',sqrt(y(1:Q.nM))./yf(1:Q.nM).*100,Q.zmsis((Q.zmsis>=Q.bot_data_high))./1000,'r');
hold off;
xlim([-30 30])
ylim([20 75])
xlabel('high gain 308 nm channel(%)')
ylabel('Altitude(km)')
%%%%%%%%%%%%%%%%%%%
m1 = Q.nM ;
m2 = Q.nL;
m3 = length(Q.real_low_308);
subplot(2,2,2)
plot(((y(m1+1:m1+m3)-yf(m1+1:m1+m3))./yf(m1+1:m1+m3)).*100,Q.zmsis(Q.zmsis<=Q.top_data_low_308)./1000,'b');
%  plot((y(1:Q.nL) - yf(1:Q.nL)./yf(1:Q.nL)).*100, Q.zmsis(Q.zmsis>=Q.bot_data_low)/1000, 'b');
hold on;
plot(-sqrt(y(m1+1:m1+m3))./yf(m1+1:m1+m3).*100,Q.zmsis(Q.zmsis<=Q.top_data_low_308)./1000,'r',sqrt(y(m1+1:m1+m3))./yf(m1+1:m1+m3).*100,Q.zmsis(Q.zmsis<=Q.top_data_low_308)./1000,'r');
hold off;
xlabel('low gain 308 nm channel(%)')
xlim([-30 30])
ylim([10 75])
ylabel('Altitude(km)')
%%%%%

% set(gca,'fontsize',18)
subplot(2,2,3)
hold on;
plot(((y(m1+m3+1:(2*m1)+m3)-yf(m1+m3+1:(2*m1)+m3))./yf(m1+m3+1:(2*m1)+m3)).*100,Q.zmsis(Q.zmsis>=Q.bot_data_high)./1000,'b');
hold on;
plot(-sqrt(y(m1+m3+1:(2*m1)+m3))./yf(m1+m3+1:(2*m1)+m3).*100,Q.zmsis((Q.zmsis>=Q.bot_data_high))./1000,'r',sqrt(y(m1+m3+1:(2*m1)+m3))./yf(m1+m3+1:(2*m1)+m3).*100,Q.zmsis((Q.zmsis>=Q.bot_data_high))./1000,'r');
hold off;
ylim([20 75])
xlim([-30 30])
xlabel('high gain 353 nm channel(%)')
ylabel('Altitude(km)')
%%%%%%%%%%%%%%%%%%%
mdata = length(yf);
subplot(2,2,4)
plot(((y(2*m1+m3+1:mdata)-yf(2*m1+m3+1:mdata))./yf(2*m1+m3+1:mdata)).*100,Q.zmsis(Q.zmsis>=Q.bot_data_low)./1000,'b');
hold on;
plot(-sqrt(y(2*m1+m3+1:mdata))./yf(2*m1+m3+1:mdata).*100,Q.zmsis(Q.zmsis>=Q.bot_data_low)./1000,'r',sqrt(y(2*m1+m3+1:mdata))./yf(2*m1+m3+1:mdata).*100,Q.zmsis(Q.zmsis>=Q.bot_data_low)./1000,'r');
hold off;
xlabel('low gain 353 nm channel(%)')
xlim([-30 30])
ylim([10 75])
ylabel('Altitude(km)')
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Vertical Resolution  %%%%%%%%%%%%%%%%%%%%


wwidth = zeros(size(Q.zret));
k = 0;
for j = 1:m % 3:m-1
    k = k + 1;
    wwidth(k) = fwhmquiet(Q.zret,X.A(1:m,j));
    if isnan(wwidth(k))
        wwidth(k) = 0;

    elseif wwidth(k) == -999
        wwidth(k) = 0;
    end
end
width = wwidth(1:end-1);
zw = Q.zret(1:end-1);
fpltw = find(width ~= 0);
handfig(4) = figure;
plot(width(fpltw)/1000,zw(fpltw)./1000)
hold on
xlabel('Vertical Resolution (km)')
ylabel('Alitude (km)')
%axis([0 300 Q.zRET(1)./1000 Q.zRET(end)./1000])
pltx = get(gca,'XLim');
plot([0 pltx(2)],[Q.zret(fini) Q.zret(fini)]./1000,'k--')
ylim([12 50])

%%%%%%%%%%%%%%%%%%%%% error for nair %%%%%%%%%%%%%%
figure;
xAlpha = X.x(m+1:2*m);
axAlpha = (xAlpha);
obsErrair = X.eo(m+1:2*m) ./ axAlpha;
plot(obsErrair(2:end-2).*100,Q.zret(2:end-2)./1000)
hold on
hleg = legend('Statistical');
set(hleg,'FontSize',8,'Box','off');
%xlabel('Optical Depth Uncertainty (%)')
xlabel('Air Density Uncertainty (%)')
ylabel('Altitude (km)')
pltx = get(gca,'XLim');
ylim([1 70])
plot(pltx,[Q.zret(fini) Q.zret(fini)]./1000,'k--')
% % % % % % % err = (X.e(2:m));
% % % % % % % upper = err + exp(X.x(2:m)); lower = exp(X.x(2:m)) - err;
% % % % % % % figure;
% % % % % % % [fillhandle,msg]=jbfilly(Q.zret(2:m)/1000,upper'./1e6,lower'./1e6,rand(1,3),rand(1,3),0, 1);
% % % % % % % hold on;
% % % % % % % plot(exp(X.x(2:m))./1e6, Q.zret(2:m)/1000, 'LineWidth', 1.5);
% % % % % % % ylim([1 70])
% % % % % % % xlim([0, 10e12]);
% % % % % % % 
% % % % % % % xlabel ('Ozone (mol/cm^3)'); ylabel('Altitude (km)');
% % % % % % % figure; plot((X.eo(2:m)./X.x(2:m))*100, Q.zret(2:end)/1000);
% % % % % % % ylim([1 70])
% % % % % % % xlabel ('error %'); ylabel('Altitude (km)');
% % % % % % % %%%%%%%%%%%%%%%% error for air %%%%%%%%%%%%%%%
% % % % % % % err = X.e(m+1:2*m);
% % % % % % % upper = err + X.x(m+1:2*m); lower = X.x(m+1:2*m) - err;
% % % % % % % title('O3')m 
% % % % % % % figure;
% % % % % % % [fillhandle,msg]=jbfilly(Q.zret(1:m)/1000,upper',lower',rand(1,3),rand(1,3),0, 1);
% % % % % % % hold on;
% % % % % % % semilogx(X.x(m+1:2*m), Q.zret(1:m)/1000, 'LineWidth', 1.5);
% % % % % % % hold on;
% % % % % % % semilogx(Q.n_air_p, Q.zret(1:m)/1000)
% % % % % % % ylim([1 70])
% % % % % % % xlabel ('air density (#/m^3)'); ylabel('Altitude (km)');
% % % % % % % figure; plot((X.eo(m+1:2*m)./X.x(m+1:2*m))*100, Q.zret(1:end)/1000);
% % % % % % % ylim([1 70])
% % % % % % % xlim([0, 70]);
% % % % % % % xlabel ('error %'); ylabel('Altitude (km)');
% % % % % % % title('air density')
% % % % % % % %%%%%%%%%%%%
% % % % % % % figure; 
% % % % % % % plot(objectdata.resolution, objectdata.alt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% beta%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % I = ones(m, m);
% % H = -(1/2)*log(abs(I-X.A(1:m, 1:m)));
% % figure; plot(H, Q.zret); hold on;  plot(sum(H), Q.zret);