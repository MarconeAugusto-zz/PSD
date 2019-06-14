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

[h, w] = freqs(b,a, logspace(-2, 3, 10000));
semilogx(w, 20*log10(abs(h)))
ylim([-80 10])
grid on
hold on
plot([0,fs,fs,2000],[0,0,-As,-As], 'r')
plot([0,fp,fp,],[-Ap,-Ap,-80], 'r')

syms p
Np(p) = poly2sym(b, p);
Dp(p) = poly2sym(a, p);
Hp(p) = Np(p) / Dp(p);
pretty(vpa(collect(Hp), 5)); % collect simplifica ao maximo a funcao

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

clear h w
figure(2)
[h,w] = freqs(bsn,asn,linspace(0,10,1000));
plot(w,20*log10(abs(h)));
ylim([-80 10])
title('H(s)')
grid on
hold on

% Fazer a mascara em cima do LAMBDA
plot([0,lambdas,lambdas,10],[0,0,-As,-As], 'r')
plot([0,lambdap,lambdap],[-Ap,-Ap,-80], 'r')

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
subplot(121)
[hz, wz] = freqz(bzn, azn, linspace(0, pi, 1000));
plot(wz/pi*fa/2, 20*log10(abs(hz)));
ylim([-80 10])
title('H(z)')
grid on
hold on
plot([0,fs,fs,2000],[0,0,-As,-As], 'r')
plot([0,fp,fp,],[-Ap,-Ap,-80], 'r')

subplot(122)
zplane(bzn, azn);

%% Projeto Filtro FIR - Janela Ajust??vel

Ap = 2; % Ganho na banda de passagem em dB
As = 30; % Atenua????o no stopband em dB
fa = 4000; % Hz
f1 = 1000; % Hz fp
f2 = 1300; % Hz fs
GdB = 5; % dB

wp = f1/(fa);
ws = f2/(fa);
dw = (2*pi)*(ws - wp);

n = ceil((As - 8)/(2.285*dw) +1);
betha = 0.5842*(As-21)^0.4+0.07886*(As-21);

Jkaiser = kaiser(n+1,betha);

wc = sqrt(wp*ws); % frequ??ncia de corte, m??dia das frequ??ncias

k = 1:n;

g0 = 0;

b1 = sin(k*wc)./(k*pi);
% b0 = sin(0*wc)./(0*pi); % sem L'Hospital
b0 = wc/pi; % L'Hospital do b0 acima
b = [flip(b1) b0 b1];
% b = b.*Jkaiser*10^(-g0/20);

% COMO OS AJUSTES FUNCIONAM?
