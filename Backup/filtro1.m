%% Dados
%   LP - (fa = 4000 Hz, f1 = 1000 Hz; f2 = 1300 Hz, Ap = 2 dB, As = 30 dB, GdB = 5 dB)
%   IIR - Butterworth, FIR - Janela Ajust??vel

%% Projeto Filtro IIR - Butterworth

Ap = 2; % Ganho na banda de passagem em dB
As = 30; % Atenua????o no stopband em dB
fa = 4000; % Hz
fp = 1000; % Hz
fs = 1300; % Hz
GdB = 5; % dB

thetap = fp/(fa/2);
thetas = fs/(fa/2);

lambdap = 2*tan((thetap*pi)/2);
lambdas = 2*tan((thetas*pi)/2);

Ws = lambdas/lambdap;
Wp = 1;

E = sqrt(10^(0.1*Ap)-1);
n = ceil((log((10^(0.1*As)-1)/E^2))/(2*log(Ws))); % ordem do filtro
%n = n-1; % teste para ver se ?? a menor ordem poss??vel
k = 1:n;
pk = E^(-1/n)*exp((1j*(2*k+n-1)/(2*n)*pi));

% k2 = real(prod(-pk));
% k3 = 1/E;
% k2, k3 e b produzem o mesmo resultado
a = real(poly(pk)); % denominador
b = a(end);         % numerador

syms p
Np(p) = poly2sym(b, p);
Dp(p) = poly2sym(a, p);
Hp(p) = Np(p) / Dp(p);
pretty(vpa(collect(Hp), 5)); % collect simplifica ao maximo a funcao

figure(1)
[h, w] = freqs(b,a, logspace(-1, 1, 10000));
semilogx(w, 20*log10(abs(h)))
ylim([-80 10])
grid on
hold on
plot([0.1,Ws,Ws,10],[0,0,-As,-As], '--r')
plot([0.1,Wp,Wp,],[-Ap,-Ap,-80], '--r')
title('H(p)')

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
plot(ws, 20*log10(abs(hs)));
ylim([-80 10])
title('H(s)')
grid on
hold on

% Fazer a mascara em cima do LAMBDA
plot([0,lambdas,lambdas,8],[0,0,-As,-As], '--r')
plot([0,lambdap,lambdap],[-Ap,-Ap,-80], '--r')

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
subplot(321)
[hz, wz] = freqz(bzn, azn, linspace(0, pi, 10000));
plot(wz/pi*fa/2, 20*log10(abs(hz)));
ylim([-80 10])
title('Resposta de magnitude de H(z)')
grid on
hold on
plot([0.01,fs,fs,2000],[0,0,-As,-As], '--r')
plot([0.01,fp,fp,],[-Ap,-Ap,-80], '--r')

% subplot(322)
% stem([flip(bi) b0 bi]); grid on;
% title('Resposta ao impulso');

subplot(3,2,[4 6])
zplane(bzn, azn);
xlabel('Parte real');
ylabel('Parte imaginaria');

subplot(323)
plot(w/pi, unwrap(angle(h))/pi); grid on;
title('Resposta de fase de H(z)');

subplot(325)
grpdelay(b, 1);title('Atraso de grupo');
xlabel('Frequencia normalizada [x\pi rad/amostra]');
ylabel('Atraso de grupo [amostra]');

%% Projeto Filtro FIR - Janela Ajustavel

clear all;
close all; %%%%%%%%%%%%%
clc;

Ap = 2; % Ganho na banda de passagem em dB
As = 30; % Atenua��o no stopband em dB
fa = 4000; % Hz
fp = 1000; % Hz
fs = 1300; % Hz
GdB = 5; % dB

ganho = GdB;

wp = fp/(fa)*(2*pi);
ws = fs/(fa)*(2*pi);
dw = (ws - wp);
wc = sqrt(wp*ws); % frequ�ncia de corte
betha = 0.5842*(As-21)^0.4+0.07886*(As-21);
n = ceil((As - 8)/(2.285*dw) +1);

if mod(n,2) == 1 
    %impar
    n = n+1;
end 

n = n-6;

% Ajuste do ganho
    GdB = GdB - 0.039;
% Primeiro ajuste de M
    wp1 = 0.5177*pi; ws1 = 0.6641*pi;
    Dw1 = ws1-wp1;
    n1 = ceil(n*Dw1/dw);
    n = n1;
    wc = wc - 0.015*pi;

Jkaiser = kaiser(n+1,betha);

k = 1:(n/2);
b1 = sin(k*wc)./(k*pi);
b0 = wc/pi; % L'Hospital do b0 acima
b = [flip(b1) b0 b1];
b = b.'; % matriz b transposta
b = b.*Jkaiser*10^(GdB/20)*10^(-0.189/20);

figure(5)
subplot(211)
[h, w] = freqz(b, 1, linspace(0,pi,10000));
plot(w/pi*fa/2, 20*log10(abs(h))); grid on;
ylim([-60 10])
hold on;title_txt = ['BP - Filtro FIR - Janela ajust�vel Kaiser - N = ' num2str(n)];
title(title_txt);xlabel('Hz');ylabel('dB');
fmax = fa/2;
% M�scara
plot([0,fs,fs,fmax],[0,0,-As,-As]+ganho, '--red')
plot([0,fp,fp],[-Ap,-Ap,-140]+ganho, '--red');hold off;

subplot(212)
plot(w/pi*fa/2, 20*log10(abs(h))); grid on;hold on;
title_txt = ['BP - Filtro FIR - Janela ajust�vel Kaiser - N = ' num2str(n)];
title(title_txt);xlabel('Hz');ylabel('dB');xlim([900 1100]); ylim([-5 10]);
% M�scara
plot([0,fs,fs,fmax],[0,0,-As,-As]+ganho, '--red')
plot([0,fp,fp],[-Ap,-Ap,-140]+ganho, '--red');hold off;

figure(6)
subplot(121)
zplane(h,1);title('Diagrama de p�los e zeros');xlabel('Parte real');ylabel('Parte imagin�ria');axis([-2 2 -3 3]);
subplot(122)
grpdelay(h,1);title('Atraso de grupo');
xlabel('Frequ�ncia normalizada [x\pi rad/amostra]');
ylabel('Atraso de grupo [amostra]');