%% Dados
%   LP - (fa = 4000 Hz, f1 = 1000 Hz; f2 = 1300 Hz, Ap = 2 dB, As = 30 dB, GdB = 5 dB)
%   IIR - Butterworth, FIR - Janela Ajustável

%% Projeto Filtro IIR - Butterworth

Ap = 2; % Ganho na banda de passagem em dB
As = 30; % Atenuação no stopband em dB
fa = 4000; % Hz
f1 = 1000; % Hz
f2 = 1300; % Hz
GdB = 5; % dB

thetap = f1/(fa/2);
thetas = f2/(fa/2);

lambdap = 2*tan((thetap*pi)/2);
lambdas = 2*tan((thetas*pi)/2);

E = sqrt(10^(0.1*Ap)-1);
n = ceil((log(10^(0.1*As)-1)/E^2)/(2*log(Ws))); % ordem do filtro
k = 1:n;
pk = E^(-1/n)*exp((1j*(2*k+n-1)/(2*n)*pi));

% CALCULAR H(P), ACHANDO PRIMEIRO O P E O D(P)

syms p;
%H(p) = 

% -------------------- Plotando --------------------
figure(1)
zplane(b,a);

figure(2)
[h, w] = freqs(b,a, logspace(-2, 3, 10000));
semilogx(w, 20*log10(abs(h)))
ylim([-80 10])
grid on
hold off


%% Projeto Filtro FIR - Janela Ajustável

Ap = 2; % Ganho na banda de passagem em dB
As = 30; % Atenuação no stopband em dB
fa = 4000; % Hz
f1 = 1000; % Hz
f2 = 1300; % Hz
GdB = 5; % dB

