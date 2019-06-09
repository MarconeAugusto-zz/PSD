%% Dados
%   HP - (fa = 4000 Hz, f1 = 1000 Hz; f2 = 1300 Hz, Ap = 0.5 dB, As = 40 dB, GdB = 0 dB)
%   IIR - Chebyshev II, FIR - Janela Fixa

% Projeto Filtro IIR - Chebyshev II

clear all;
close all;
clc;

Ap = 0.5; % Ganho na banda de passagem em dB
As = 40; % Atenuação no stopband em dB
fa = 4000; % Hz
f1 = 1000; % Hz
f2 = 1300; % Hz
GdB = 0; % dB


f = [1000 1300]; % frequências em Hz

% substituindo de Hz para ômega
% w = f/fa*(2*pi); %w = 2*pi*f
% ws1 = w(1)/pi;
% wp1 = w(2)/pi;

% teste
fN = fa/2; %frequencia de Niquist
w = 2*pi*f; 
wp1 = w(1)/fN;
ws1 = w(2)/fN;

Ws = ws1/wp1; %Lowpass
Wp = 1;

% Filtro Chebyshev 2
Rp = Ap; Rs = As;
[n,Wn] = cheb2ord(Wp,Ws,Rp,Rs,'s')
[b,a] = cheby2(n,Rp,Wn,'s');
[b1,a1] = cheby2(n,Rp,1,'s');

% Plot filtro PB
figure(1)
[h1,w1] = freqs(b,a,logspace(-2,1,1000));
semilogx(w1,20*log10(abs(h1)));grid on; %ylim([-60 5]);
hold on;
[h1,w1] = freqs(b1,a1,logspace(-2,1,1000));
semilogx(w1,20*log10(abs(h1)));grid on; %ylim([-60 5]);
hold off;

%%
% Transformação de frequência Lowpass para Highpass
ap = a1; bp = b1; 
syms p;
Np(p) = poly2sym(bp, p);
Dp(p) = poly2sym(ap, p);
Hp(p) = Np(p) / Dp(p);
pretty(vpa(collect(Hp(p)), 5))

% transformação de frequência
syms s;
Hs(s) = collect(subs(Hp(p),wp1/s));%transformação lowpass/bandstop
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
%subplot(211)
[h, w] = freqs(bsn,asn, linspace(0, fa*pi, 10000));
plot(w/2/pi, 20*log10(abs(h))); grid on;hold on;
title_txt = ['BP - Filtro IIR - CHEBYSHEV I - N = ' num2str(n)];
title(title_txt);
%%

% Projeto Filtro FIR - Janela Fixa

Ap = 0.5; % Ganho na banda de passagem em dB
As = 40; % Atenuação no stopband em dB
fa = 4000; % Hz
f1 = 1000; % Hz
f2 = 1300; % Hz
GdB = 0; % dB