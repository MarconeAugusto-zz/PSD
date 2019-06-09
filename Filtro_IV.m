%% Dados
%   BS - (fa = 6000 Hz, f1 = 1200 Hz; f2 = 1250 Hz, f3 = 1300 Hz; f4 = 1400 Hz, Ap = 0.5 dB, As = 60 dB, GdB = 0 dB)
%   IIR - Chebyshev I, FIR - PM

% Projeto Filtro IIR - Chebyshev I

clear all;
close all;
clc;

Ap = 0.5; % Ganho na banda de passagem em dB
As = 60; % Atenuação no stopband em dB
fa = 6000; % Hz
f1 = 1200; % Hz
f2 = 1250; % Hz
f3 = 1300; % Hz
f4 = 1400; % Hz
GdB = 0; % dB

%f0 = sqrt(((f3+f4)/2)*f1); 
%B = ((f3+f4)/2)-f1; %banda de rejeiçaõ em Hz
f0 = sqrt(f4*f1); 
B = f4-f1; %banda de rejeiçaõ em Hz

f = [1200 1250 1300 1400]; % frequências em Hz

% substituindo de Hz para ômega
% w = f/fa*(2*pi); %w = 2*pi*f
% ws1 = w(1)/pi;
% wp1 = w(2)/pi;
% wp2 = w(3)/pi;
% ws2 = w(4)/pi;

% teste
fN = fa/2; %frequencia de Niquist
w = 2*pi*f; 
wp1 = w(1)/fN;
ws1 = w(2)/fN;
ws2 = w(3)/fN;
wp2 = w(4)/fN;

w0 = 2*pi*f0;
Bw = 2*pi*B;


%Ws2 =  abs(w0^2 - ws1^2)/(Bw*ws1); %utilidade disso ???
%Ws1 = abs(w0^2 - ws2^2)/(Bw*ws2);  %utilidade disso ???
% Ws = min(ws2,ws1);
% Wp = 1;
Ws = min(ws2,ws1);
Wp = 1;

% Filtro Chebyshev 1
Rp = Ap; Rs = As;
[n,Wn] = cheb1ord(Wp,Ws,Rp,Rs,'s')
[b,a] = cheby1(n,Rp,Wn,'s');

% Plot filtro PB
figure(1)
[h1,w1] = freqs(b,a,logspace(-2,1,1000));
semilogx(w1,20*log10(abs(h1)));grid on; ylim([-60 5]);

% Transformação de frequência Lowpass para Bandstop
ap = a; bp = b; 
syms p;
Np(p) = poly2sym(bp, p);
Dp(p) = poly2sym(ap, p);
Hp(p) = Np(p) / Dp(p);
pretty(vpa(collect(Hp(p)), 5))

% transformação de frequência
syms s;
Hs(s) = collect(subs(Hp(p),((Bw*s)/(s^2 + w0^2))));%transformação lowpass/bandstop
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
plot(w/2/pi, 20*log10(abs(h))); grid on;hold on;ylim([-80 5]);xlim([1000 1600]);
title_txt = ['BP - Filtro IIR - CHEBYSHEV I - N = ' num2str(n)];
title(title_txt);
% Máscara do filtro projetado
Amin = 5;
plot([0,f1,f1,f4,f4,fa/2],-[Ap,Ap,80,80,Ap,Ap], '--r');
plot([f2,f2,f3,f3],[Amin,-As,-As,Amin], '--m');
hold off;

subplot(212)
[h, w] = freqs(bsn,asn, linspace(0, fa*pi, 10000));
plot(w/2/pi, 20*log10(abs(h)));grid on;hold on; ylim([-65 -55]);xlim([1180 1420]);
title_txt = ['BP - Filtro IIR - CHEBYSHEV I - N = ' num2str(n)];
title(title_txt);
% Máscara do filtro projetado
Amin = 5;
plot([0,f1,f1,f4,f4,fa/2],-[Ap,Ap,80,80,Ap,Ap], '--r');
plot([f2,f2,f3,f3],[Amin,-As,-As,Amin], '--m');
hold off;

figure(3)
subplot(1,2,1)
grpdelay(b,a);
subplot(1,2,2)
zplane(b,a);grid on;


%%
% Projeto Filtro FIR - PM

clear all;
clc;

Ap = 0.5; % Ganho na banda de passagem em dB
As = 60; % Atenuação no stopband em dB
fa = 6000; % Hz
f1 = 1200; % Hz
f2 = 1250; % Hz
f3 = 1300; % Hz
f4 = 1400; % Hz
GdB = 0; % dB

f = [1200 1250 1300 1400]; % frequências em Hz
w = f/fa*(2*pi);
ws1 = w(1)/pi;
wp1 = w(2)/pi;
wp2 = w(3)/pi;
ws2 = w(4)/pi;
mags = [1 0 1];

% [n,fo,ao,w] = firpmord(f,a,dev,fs)
devAs = 10^(-(As-3.5)/20);
devAp = 1-10^(-(Ap/2-0.05)/20);
devs = [devAp devAs devAp];

% calculo da ordem com firpmord
f = f + [0 0 30 0];
[n,f0,a0,w0] = firpmord(f,mags,devs,fa);

%G0 = -Ap/2;
% calculo algoritmo PM
h_pm = firpm(n,f0,a0,w0);
%h_pm = h_pm*10^(G0/20);

%clear Hw W
figure(4)
subplot(211)
[Hw,w] = freqz(h_pm,1,10000);
plot(w*fa/2/pi,20*log10(abs(Hw)))
title_txt = ['BS - Filtro FIR - PM - N = ' num2str(n)];
title(title_txt);
hold on;
% Máscara
Amin = 5;
plot([0,ws1,ws1,ws2,ws2,1]*fa/2,-[0,0,80,80,0,0], '--r'); ylim([-80 5]);xlim([1000 1600]);
plot([wp1,wp1,wp2,wp2]*fa/2,[Amin,-As,-As,Amin], '--m'); grid minor;
hold off;
subplot(212)
[Hw,w] = freqz(h_pm,1,10000);
plot(w*fa/2/pi,20*log10(abs(Hw)))
title_txt = ['BS - Filtro FIR - PM - N = ' num2str(n)];
title(title_txt);
hold on;
% Máscara
Amin = 5;
plot([0,ws1,ws1,ws2,ws2,1]*fa/2,-[0,0,80,80,0,0], '--r'); ylim([-65 -55]);xlim([1180 1420]);
plot([wp1,wp1,wp2,wp2]*fa/2,[Amin,-As,-As,Amin], '--m'); grid minor;
hold off;


figure(5);
subplot(1,2,1)
grpdelay(h_pm);
subplot(1,2,2)
zplane(h_pm,1);


