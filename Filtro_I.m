% %% Dados
% %   LP - (fa = 4000 Hz, f1 = 1000 Hz; f2 = 1300 Hz, Ap = 2 dB, As = 30 dB, GdB = 5 dB)
% %   IIR - Butterworth, FIR - Janela Ajustável

%% Projeto Filtro IIR - Butterworth

clear aal;
close all;
clc;

Ap = 2; % Ganho na banda de passagem em dB
As = 30; % Atenuação no stopband em dB
fa = 4000; % Hz
f1 = 1000; % Hz
f2 = 1300; % Hz
GdB = 5; % dB

fp = f1;
fs = f2;

thetap = fp/(fa/2);
thetas = fs/(fa/2);

lambdap = 2*tan((thetap*pi)/2);
lambdas = 2*tan((thetas*pi)/2);

Ws = lambdas/lambdap;
Wp = 1;

E = sqrt(10^(0.1*Ap)-1);
n = ceil((log((10^(0.1*As)-1)/E^2))/(2*log(Ws))); % ordem do filtro
k = 1:n;
pk = E^(-1/n)*exp((1j*(2*k+n-1)/(2*n)*pi));

a = real(poly(pk)); % denominador
b = a(end);         % numerador

syms p
Np(p) = poly2sym(b, p);
Dp(p) = poly2sym(a, p);
Hp(p) = Np(p) / Dp(p);
pretty(vpa(collect(Hp), 5)); % collect simplifica ao maximo a funcao

figure(1)
[h, w] = freqs(b,a, logspace(-1, 1, 10000));
semilogx(w, 20*log10(abs(h)));ylim([-60 10]);grid on;hold on;
% Máscara
plot([0.1,Ws,Ws,10],[0,0,-As,-As], '--r')
plot([0.1,Wp,Wp,],[-Ap,-Ap,-80], '--r')
title('H(p)');xlabel('rad/s');ylabel('dB');

syms s
Hs(s) = collect(subs(Hp(p), s/lambdap));
[Ns, Ds] = numden(Hs(s));
pretty(vpa(Hs(s), 3))

bs = sym2poly(Ns);
as = sym2poly(Ds);

an = as(1);
bsn = bs/an;
asn = as/an;

Hsn(s) = poly2sym(bsn,s) / poly2sym(asn,s);
pretty(vpa(Hsn(s),5))

figure(2)
[hs,ws] = freqs(bsn,asn,linspace(0,8,10000));
plot(ws/pi, 20*log10(abs(hs)));ylim([-60 10]);
title('H(s)');xlabel('rad/s');ylabel('dB');
grid on; hold on;
% Fazer a mascara em cima do LAMBDA
plot([0,lambdas/pi,lambdas/pi,2],[0,0,-As,-As], '--r')
plot([0,lambdap/pi,lambdap/pi],[-Ap,-Ap,-80], '--r')

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
[hz, wz] = freqz(bzn, azn, linspace(0, pi, 10000));
plot(wz/pi*fa/2, 20*log10(abs(hz)));ylim([-40 10]);
title_txt = ['H(z) - BP - Filtro IIR - Butterworth - N = ' num2str(n)];
title(title_txt);xlabel('Hz');ylabel('dB');
grid on;hold on;
% Máscara
plot([0.01,fs,fs,2000],[0,0,-As,-As], '--r')
plot([0.01,fp,fp,],[-Ap,-Ap,-80], '--r')

subplot(212)
plot(wz/pi*fa/2, 20*log10(abs(hz)));ylim([-5 5]);xlim([800 1100]);
title_txt = ['H(z) - BP - Filtro IIR - Butterworth - N = ' num2str(n)];
title(title_txt);xlabel('Hz');ylabel('dB');
grid on;hold on;
% Máscara
plot([0.01,fs,fs,2000],[0,0,-As,-As], '--r')
plot([0.01,fp,fp,],[-Ap,-Ap,-80], '--r')

figure(4)
subplot(121)
zplane(bzn, azn);title('Diagrama de pólos e zeros');xlabel('Parte real');ylabel('Parte imaginária');axis([-2 2 -3 3]);
subplot(122)
grpdelay(bzn, azn);title('Atraso de grupo');
xlabel('Frequência normalizada [x\pi rad/amostra]');
ylabel('Atraso de grupo [amostra]');

%% Projeto Filtro FIR - Janela Ajustável

clear all;
clc;

Ap = 2; % Ganho na banda de passagem em dB
As = 30; % Atenuação no stopband em dB
fa = 4000; % Hz
f1 = 1000; % Hz
f2 = 1300; % Hz
GdB = 5; % dB

%G0 = GdB;
G0 = 0; %%%%%%%%%%%%%%%%%teste

f = [1000 1300]; % frequências em Hz
% substituindo de Hz para ômega
w = f/fa*(2*pi); %w = 2*pi*f
wp = w(1);
ws = w(2);
Dw = ws - wp;   %largura da banda de transição
wc = sqrt(wp*ws); % frequência de corte

M =9;

% %Ajuste do ganho
    %levar o pico para abaixo de 0
    ganho = 0.321-GdB; %ganho dB mediddo no plot do filtro
    G0 = ganho;
    
% primeiro ajuste de M (N/2)
    wp1 = 0.5183*pi; ws1 = 0.6641*pi; % valores medidos no gráfico
    Dw1 = ws1 - wp1;
    M1 = ceil(M*Dw1/Dw);
    M = M1; 
    wc = wc - 0.015*pi;

betha = 0.5842*(As-21)^0.4 + 0.07886*(As-21); %janela Kaiser, usa essa formula pois As = 30
N = ceil((As - 8)/(2.285*Dw)+1);
wkaiser = kaiser(N, betha);

%M = N; %como determinar o M ? é arbitrario ?

N = 2*M+1;
wcheb = chebwin(N, As-8)';

k = 1:M;

bi = sin(wc*k)./(pi*k);
b0 = wc/pi;
b = [flip(bi) b0 bi];

m = -M:M;
b = b.*wcheb*10^(-G0/20); % janela de keiser

subplot(3,2,[4 6])
zplane(b, 1);
axis([-2 2 -2 2])
[h, w] = freqz(b, 1, 'whole');
subplot(322)
stem(b), grid on;
subplot(321)
[h, w] = freqz(b, 1, linspace(0,pi,10000));
% plot(w/pi, abs(h)); grid on;
plot(w/pi, 20*log10(abs(h))); grid on; ylim([-80 10]);
hold on;
plot([0,wp,wp]/pi,[-Ap,-Ap,-85]+GdB, '-red')
plot([0,ws/pi,ws/pi,1],[0,0,-As,-As]+GdB, '-red')


subplot(323)
plot(w/pi, unwrap(angle(h))/pi); grid on;
subplot(325)
grpdelay(b, 1)

figure(2)
[h, w] = freqz(b, 1, linspace(0,pi,10000));
% plot(w/pi, abs(h)); grid on;
plot(w/pi, 20*log10(abs(h))); grid on;
ylim([-80 10])
hold on;
plot([0,wp,wp]/pi,[-Ap,-Ap,-85]+GdB, '-red')
plot([0,ws/pi,ws/pi,1],[0,0,-As,-As]+GdB, '-red')



