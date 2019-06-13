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

syms p
Np(p) = poly2sym(b, p);
Dp(p) = poly2sym(a, p);
Hp(p) = Np(p) / Dp(p);
pretty(vpa(collect(Hp), 5)); % collect simplifica ao maximo a funcao

syms s
Hs(s) = collect(subs(Hp, s/Wp));
[Ns, Ds] = numden(Hs(s));
pretty(vpa(Hs(s), 3))

syms z
s = 2*fa*((z-1)/(z+1));


%%

% -------------------- Plotando --------------------
figure(1)
zplane(b,a);

figure(2)
[h, w] = freqs(b,a, logspace(-2, 3, 10000));
semilogx(w, 20*log10(abs(h)))
ylim([-80 10])
grid on
hold off

%% Projeto Filtro FIR - Janela Ajust??vel

Ap = 2; % Ganho na banda de passagem em dB
As = 30; % Atenua????o no stopband em dB
fa = 4000; % Hz
f1 = 1000; % Hz fp
f2 = 1300; % Hz fs
GdB = 5; % dB

w1 = f1/(fa/2);
w2 = f2/(fa/2);
dw = w2 - w1;

n = ceil((As - 8)/(2.285*dw) +1);
betha = 0.5842*(As-21)^0.4+0.07886*(As-21);

Jkaiser = kaiser(n+1,betha);

wc = sqrt(wp*ws); % frequ??ncia de corte, m??dia das frequ??ncias

k = 1:n;

b1 = sin(k*wc)./(k*pi);
% b0 = sin(0*wc)./(0*pi); % sem L'Hospital
b0 = wc/pi; % L'Hospital do b0 acima
b = [flip(b1) b0 b1];

% COMO OS AJUSTES FUNCIONAM?

%%

fcuts = [fp fs];

w = fcuts/fa*(2*pi);
wp = w(1)/pi;
ws = w(2)/pi;
mags = [1 0];
devs = [1-10^(-Ap/20) 10^(-As/20)];

[n,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,fa);
