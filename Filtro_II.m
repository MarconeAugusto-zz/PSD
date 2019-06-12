% %% Dados
% %   HP - (fa = 4000 Hz, f1 = 1000 Hz; f2 = 1300 Hz, Ap = 0.5 dB, As = 40 dB, GdB = 0 dB)
% %   IIR - Chebyshev II, FIR - Janela Fixa
% 
% % Projeto Filtro IIR - Chebyshev II

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
w = f/fa*(2*pi); %w = 2*pi*f
wp = w(2)/pi;
ws = w(1)/pi;

Ws = wp/ws; %Lowpass -> Highpass
Wp = 1;
% Ws = 1;
% Wp = ws/wp;

% Filtro Chebyshev 2
Rp = Ap; Rs = As;
[n,Wn] = cheb2ord(Wp,Ws,Rp,Rs,'s')
[b,a] = cheby2(n,Rp,Wn,'s');
[b1,a1] = cheby2(n,Rp,1,'s');

% Plot filtro PB
figure(1)
[h1,w1] = freqs(b,a,logspace(-2,1,1000));
semilogx(w1,20*log10(abs(h1)));grid on; ylim([-40 5]);
hold on;
[h1,w1] = freqs(b1,a1,logspace(-2,1,1000));
semilogx(w1,20*log10(abs(h1)));grid on; %ylim([-60 5]);
hold off;


% Transformação de frequência Lowpass para Highpass
ap = a; bp = b; 
syms p;
Np(p) = poly2sym(bp, p);
Dp(p) = poly2sym(ap, p);
Hp(p) = Np(p) / Dp(p);
pretty(vpa(collect(Hp(p)), 5))

% transformação de frequência
syms s;
Hs(s) = collect(subs(Hp(p),wp/s));%transformação lowpass/bandstop
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
plot([0 wp/pi wp/pi 1],[-As -As 0 0], '--red')
plot([ws/pi,ws/pi,1],[-60 -Ap,-Ap], '--red')
%%

clear all;
close all;
clc;

% Projeto Filtro FIR - Janela Fixa

Ap = 0.5; % Ganho na banda de passagem em dB
As = 40; % Atenuação no stopband em dB
fa = 4000; % Hz
f1 = 1000; % Hz
f2 = 1300; % Hz
GdB = 0; % dB

g0 = GdB;

f = [1000 1300]; % frequências em Hz

% substituindo de Hz para ômega
w = f/fa*(2*pi); %w = 2*pi*f
ws = w(2);
wp = w(1);

%projeto original
Dw = ws - wp;
%M = ceil(3.32*pi/Dw); % ordem (3.32 tabela Hamming)
M = ceil(3.11*pi/Dw); % ordem (3.11 tabela Hann)

%Ajuste do ganho
    %levar o pico para abaixo de 0
    ganho = 0.05501; %ganho dB mediddo no plot do filtro
    g0 = ganho;
% 
% primeiro ajuste de M (N/2)
    wp1 = 0.4998*pi; ws1 = ws; % valores medidos no gráfico
    Dw1 = ws1 - wp1;
    M2 = ceil(M*Dw1/Dw);
    M = M2; 
%    
% % segundo ajuste de M (N/2)
%     wp2 = 0.4463*pi; ws2 = 0.5975*pi;
%     Dw2 = ws2 - wp2;
%     M3 = ceil(M*Dw2/Dw);
%     M = M3;
%   
% % Terceiro ajuste da frequência de corte
%     Dwp = 0.4*pi - wp2;
%     Dws = 0.6*pi - ws2;
%     wc3 = wc + (Dwp + Dws)/2;
%     wc = wc3;


wc = sqrt(wp*ws);   % frequência de corte, média das frequências

k = 1:M;

% Highpass
bi = -sin(k*wc)./(k*pi); 
b0 = 1 - ( wc/pi);
b = [flip(bi) b0 bi];

m = -M : M;
%dw = 0.04; % filtro hamming
dw = 0.00; % filtro hann
wk = (0.5+dw)+(0.5-dw)*cos(2*pi*m/(2*M+1)); %serve para filtro hamming ou hann
b = b.*wk*10^(-g0/20);

% Acertar a janela no plot
ws = w(2);
wp = w(1);

subplot(321)
[h, w] = freqz(b,1,linspace(0,pi,10000)); 
%plot(w/pi,abs(h)); grid on; xlim([0 1])    %dominio do tempo
hold on;
plot(w/pi,20*log10(abs(h))); grid on; xlim([0 1]);title('Resposta de magnitude de H(z)');ylim([-80 5]);
plot([0 wp/pi wp/pi 1],[-As -As 0 0], '--red')
plot([ws/pi,ws/pi,1],[-60 -Ap,-Ap], '--red')
hold off;

subplot(322)
stem([flip(bi) b0 bi]); grid on;title('Resposta ao impulso');

subplot(3,2,[4 6])
zplane(b, 1)
axis([-2 2 -2 2]);

subplot(323)
plot(w/pi, unwrap(angle(h))/pi); grid on;title('Resposta de fase de H(z)');
subplot(325)
grpdelay(b, 1);title('Atraso de grupo');

figure(2)
hold on;
plot(w/pi,20*log10(abs(h))); grid on; xlim([0 1]);title('Resposta de magnitude de H(z)');ylim([-80 5]);
plot([0 wp/pi wp/pi 1],[-As -As 0 0], '--red')
plot([ws/pi,ws/pi,1],[-60 -Ap,-Ap], '--red')
hold off;