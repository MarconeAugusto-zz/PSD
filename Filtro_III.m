%% Dados
%   BP - (fa = 4000 Hz, f1 = 1000 Hz; f2 = 1200 Hz, f3 = 1250 Hz; f4 = 1300 Hz, Ap = 1 dB, As = 20 dB, GdB = 0 dB)
%   IIR - Eliptico, FIR - PM
%%

clear all;
close all;
clc;

% Projeto Filtro IIR - Eliptico

Ap = 1; % Ganho na banda de passagem em dB
As = 20; % Atenuação no stopband em dB
fa = 4000; % Hz
f1 = 1000; % Hz
f2 = 1200; % Hz
f3 = 1250; % Hz
f4 = 1300; % Hz
GdB = 0; % dB

f0 = sqrt(f2*f3); 
B = f3-f2; %banda de passagem em Hz

f = [1000 1200 1250 1300]; % frequências em Hz
% substituindo de Hz para ômega
% w = f/fa*(2*pi); %w = 2*pi*f
% ws1 = w(1)/pi;
% wp1 = w(2)/pi;
% wp2 = w(3)/pi;
% ws2 = w(4)/pi;

% teste
fN = fa/2; %frequencia de Niquist
w = 2*pi*f; 
ws1 = w(1)/fN;
wp1 = w(2)/fN;
wp2 = w(3)/fN;
ws2 = w(4)/fN;

w0 = 2*pi*f0;
Bw = 2*pi*B;

%Ws2 =  abs(w0^2 - ws1^2)/(Bw*ws1); %utilidade disso ???
%Ws1 = abs(w0^2 - ws2^2)/(Bw*ws2);  %utilidade disso ??? 
Ws = min(ws2,ws1);
Wp = 1;

% Filtro elíptico
[n,Wn] = ellipord(Wp,Ws,Ap,As,'s')
[b,a] = ellip(n,Ap,As,Wn,'s');

% Plot filtro PB
figure(1)
[h1,w1] = freqs(b,a,logspace(-2,1,1000));
semilogx(w1,20*log10(abs(h1)));grid on; ylim([-60 5]);

% Transformação de frequência Lowpass para Bandpass
ap = a; bp = b; 
syms p;
Np(p) = poly2sym(bp, p);
Dp(p) = poly2sym(ap, p);
Hp(p) = Np(p) / Dp(p);
pretty(vpa(collect(Hp(p)), 5))

% transformação de frequência
syms s;
Hs(s) = collect(subs(Hp(p), (s^2 + w0^2)/(Bw*s))); %transformação lowpass/bandpass
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
subplot(211)
[h, w] = freqs(bsn,asn, linspace(0, fa*pi, 10000));
plot(w/2/pi, 20*log10(abs(h))); grid on;hold on;ylim([-60 5]);xlim([800 1500]);
title_txt = ['BP - Filtro IIR - ELÍPTICO - N = ' num2str(n)];
title(title_txt);
% Máscara do filtro projetado
Amin = 40;
plot([0,f1,f1,f4,f4,fa/2],-[As,As,0,0,As,As], '--r'); 
plot([f2,f2,f3,f3],-[Amin,Ap,Ap,Amin], '--m'); 
hold off;

subplot(212)
plot(w/2/pi, 20*log10(abs(h))); grid on;hold on;ylim([-5 2]);xlim([998 1302]);
title_txt = ['BP - Filtro IIR - ELÍPTICO - N = ' num2str(n)];
title(title_txt);
% Máscara do filtro projetado
Amin = 40;
plot([0,f1,f1,f4,f4,1],-[As,As,0,0,As,As], '--r'); 
plot([f2,f2,f3,f3],-[Amin,Ap,Ap,Amin], '--m'); 
hold off;

figure(3)
subplot(1,2,1)
grpdelay(b,a);
subplot(1,2,2)
zplane(b,a);grid on;



%% Projeto BP - Filtro FIR - PM
%close all;
clear all;
clc;

Ap = 1; % Ganho na banda de passagem em dB
As = 20; % Atenuação no stopband em dB
fa = 4000; % Hz
f1 = 1000; % Hz
f2 = 1200; % Hz
f3 = 1250; % Hz
f4 = 1300; % Hz
GdB = 0; % dB

f = [1000 1200 1250 1300]; % frequências em Hz
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
% Máscara
plot([0,ws1,ws1,ws2,ws2,1]*fa/2,-[As,As,0,0,As,As], '--r'); ylim([-60 5]);xlim([800 1500]);
plot([wp1,wp1,wp2,wp2]*fa/2,-[Amin,Ap,Ap,Amin], '--m'); grid on;
hold off;
subplot(212);
[Hw,w] = freqz(h_pm,1,10000);
plot(w*fa/2/pi,20*log10(abs(Hw)));ylim([-Amin 5]);
title_txt = ['BP - Filtro FIR - PM - N = ' num2str(n)];
title(title_txt);
hold on
% Máscara
plot([0,ws1,ws1,ws2,ws2,1]*fa/2,-[As,As,0,0,As,As], '--r'); ylim([-5 2]);xlim([998 1302]);
plot([wp1,wp1,wp2,wp2]*fa/2,-[Amin,Ap,Ap,Amin], '--m'); grid on;
hold off;

figure(5);
subplot(1,2,1)
grpdelay(h_pm);
subplot(1,2,2)
zplane(h_pm,1);



