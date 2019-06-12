%% Dados
%   BP - (fa = 4000 Hz, f1 = 1000 Hz; f2 = 1200 Hz, f3 = 1250 Hz; f4 = 1300 Hz, Ap = 1 dB, As = 20 dB, GdB = 0 dB)
%   IIR - Eliptico, FIR - PM
%%

clear all;
close all;
clc;

% Projeto Filtro IIR - Eliptico

Ap = 1; % Ganho na banda de passagem em dB
As = 20; % Atenua��o no stopband em dB
fa = 4000; % Hz
f1 = 1000; % Hz
f2 = 1200; % Hz
f3 = 1250; % Hz
f4 = 1300; % Hz
GdB = 0; % dB

% Frequ�ncias
fs1 = f1;
fp1 = f2;
fp2 = f3;
fs2 = f4;
% Frequ�ncia -> omega
ws1 = 2*pi*fs1;
wp1 = 2*pi*fp1;
wp2 = 2*pi*fp2;
ws2 = 2*pi*fs2;
wa = fa*2*pi;
% C�lculo do tetha
tetha_s1 = ws1/(wa/2);
tetha_p1 = wp1/(wa/2);
tetha_s2 = ws2/(wa/2);
tetha_p2 = wp2/(wa/2);
% C�lculo do lambda
lambda_s1 = 2*tan(tetha_s1 * pi/2);
lambda_s2 = 2*tan(tetha_s2 * pi/2);
lambda_p1 = 2*tan(tetha_p1 * pi/2);
lambda_p2 = 2*tan(tetha_p2 * pi/2);
lambda_0 = sqrt(lambda_p2*lambda_p1);
lambda_s = min(lambda_s1,lambda_s2);
% Lowpass -> Bandpass
B = lambda_p2 - lambda_p1;
Os = abs((-lambda_s^2+lambda_0^2)/(B*lambda_s));
Op = 1;

% Filtro el�ptico
[n,Wn] = ellipord(Op,Os,Ap,As,'s')
[b,a] = ellip(n,Ap,As,Wn,'s');

% Plot prot�tipo filtro PB
figure(1)
[h1,w1] = freqs(b,a,logspace(-2,2,10000));
semilogx(w1,20*log10(abs(h1)));grid on; ylim([-30 5]);
title('H(p)');hold on;grid on;
plot([10^-2,Os,Os,10^2],[0,0,-As,-As], '--r')
plot([10^-2,1,1],[-Ap,-Ap,-80], '--r')

% Transforma��o de frequ�ncia Lowpass para Bandpass
ap = a; bp = b; 
syms p;
Np(p) = poly2sym(bp, p);
Dp(p) = poly2sym(ap, p);
Hp(p) = Np(p) / Dp(p);
pretty(vpa(collect(Hp(p)), 5))

% transforma��o de frequ�ncia
syms s;
Hs(s) = collect(subs(Hp(p), (s^2 + lambda_0^2)/(B*s))); %transforma��o lowpass/bandpass
[N, D] = numden(Hs(s));
pretty(vpa(Hs(s), 5))

% Normalizando de acordo com p^n
bs = sym2poly(N);
as = sym2poly(D);
an = as(1);
bsn = bs/an;
asn = as/an;
Hsn(s) = poly2sym(bsn, s)/poly2sym(asn, s);
pretty(vpa(Hsn(s), 5))

% Plot filtro BP
figure(2)
[h, w] = freqs(bsn,asn,linspace(0, 100, 10000));
plot(w/pi, 20*log10(abs(h))); grid on;hold on;ylim([-60 5]);xlim([0 2]);
title('H(s)')
% Fazer a mascara em cima do LAMBDA
plot([0,lambda_s1/pi,lambda_s1/pi,lambda_s2/pi,lambda_s2/pi,2],-[As,As,0,0,As,As], '--r')
plot([lambda_p1/pi,lambda_p1/pi,lambda_p2/pi,lambda_p2/pi],-[80,Ap,Ap,80], '--r')
hold off;

syms z;
aux = 2*((z-1)/(z+1));
Hz(z) = collect(subs(Hs(s), aux));
pretty(vpa(Hz(z),3))

[Nz,Dz] = numden(Hz(z));
bz = sym2poly(Nz);
az = sym2poly(Dz);

an = az(1);
bzn = bz/an;
azn = az/an;

Hzn(z) = poly2sym(bzn,z) / poly2sym(azn,z);
pretty(vpa(Hzn(z),5))


figure(3)
subplot(211)
[hz, wz] = freqz(bzn,azn, linspace(0, pi, 1000));
plot(wz/pi*fa/2, 20*log10(abs(hz))); grid on;hold on;ylim([-60 5])
title_txt = ['BP - Filtro IIR - EL�PTICO - N = ' num2str(n)];
title(title_txt);
% M�scara do filtro projetado
plot([0,f1,f1,f4,f4,2000],-[As,As,0,0,As,As], '--r')
plot([f2,f2,f3,f3],-[40,Ap,Ap,40], '--r')
hold off;

subplot(212)
plot(wz/pi*fa/2, 20*log10(abs(hz))); grid on;hold on;ylim([-5 2]);xlim([998 1302]);
title_txt = ['BP - Filtro IIR - EL�PTICO - N = ' num2str(n)];
title(title_txt);
% M�scara do filtro projetado
Amin = 40;
plot([0,f1,f1,f4,f4,1],-[As,As,0,0,As,As], '--r'); 
plot([f2,f2,f3,f3],-[Amin,Ap,Ap,Amin], '--m'); 
hold off;

figure(4)
subplot(1,2,1)
grpdelay(bzn,azn);
subplot(1,2,2)
zplane(bzn,azn);grid on;

%% Projeto BP - Filtro FIR - PM
%close all;
clear all;
clc;

Ap = 1; % Ganho na banda de passagem em dB
As = 20; % Atenua��o no stopband em dB
fa = 4000; % Hz
f1 = 1000; % Hz
f2 = 1200; % Hz
f3 = 1250; % Hz
f4 = 1300; % Hz
GdB = 0; % dB

f = [1000 1200 1250 1300]; % frequ�ncias em Hz
w = f/fa*(2*pi);
ws1 = w(1)/pi;
wp1 = w(2)/pi;
wp2 = w(3)/pi;
ws2 = w(4)/pi;
mags = [0 1 0];

% [n,fo,ao,w] = firpmord(f,a,dev,fs)
devAs = 10^(-(As+0.5)/20);
devAp = 1-10^(-(Ap/2-0.05)/20);
devs = [devAs devAp devAs];

% calculo da ordem com firpmord
%f = f + [0 -170 0 0];
f = f + [0 -150 0 0];
[n,f0,a0,w0] = firpmord(f,mags,devs,fa);

G0 = -Ap/2;
% calculo algoritmo PM
h_pm = firpm(n,f0,a0,w0);
h_pm = h_pm*10^(G0/20);

%clear Hw W
Amin = 60;
figure(4)
subplot(211);
[Hw,w] = freqz(h_pm,1,10000);
plot(w*fa/2/pi,20*log10(abs(Hw)));ylim([-Amin 5]);
title_txt = ['BP - Filtro FIR - PM - N = ' num2str(n)];
title(title_txt);
hold on
% M�scara
plot([0,ws1,ws1,ws2,ws2,1]*fa/2,-[As,As,0,0,As,As], '--r'); ylim([-60 5]);xlim([800 1500]);
plot([wp1,wp1,wp2,wp2]*fa/2,-[Amin,Ap,Ap,Amin], '--m'); grid on;
hold off;
subplot(212);
[Hw,w] = freqz(h_pm,1,10000);
plot(w*fa/2/pi,20*log10(abs(Hw)));ylim([-Amin 5]);
title_txt = ['BP - Filtro FIR - PM - N = ' num2str(n)];
title(title_txt);
hold on
% M�scara
plot([0,ws1,ws1,ws2,ws2,1]*fa/2,-[As,As,0,0,As,As], '--r'); ylim([-5 2]);xlim([998 1302]);
plot([wp1,wp1,wp2,wp2]*fa/2,-[Amin,Ap,Ap,Amin], '--m'); grid on;
hold off;

figure(5);
subplot(1,2,1)
grpdelay(h_pm);
subplot(1,2,2)
zplane(h_pm,1);



