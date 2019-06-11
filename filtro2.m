%% Dados
%   HP - (fa = 4000 Hz, f1 = 1000 Hz; f2 = 1300 Hz, Ap = 0.5 dB, As = 40 dB, GdB = 0 dB)
%   IIR - Chebyshev II, FIR - Janela Fixa

% Projeto Filtro IIR - Chebyshev II

Ap = 0.5; % Ganho na banda de passagem em dB
As = 40; % Atenuação no stopband em dB
fa = 4000; % Hz
f1 = 1000; % Hz
f2 = 1300; % Hz
GdB = 0; % dB

% para Ap = 3dB, E = 1;
E = sqrt(10^(0.1*Ap)-1);

thetap = f1/(fa/2);
thetas = f2/(fa/2);

lambdap = 2*tan((thetap*pi)/2);
lambdas = 2*tan((thetas*pi)/2);

Ws = lambdas/lambdap;
Wp = 1;

n = (acosh(sqrt(10^(0.1*As)-1)/(10^(0.1*Ap)-1))/acosh(Ws))
n = ceil(n); % obtem o proximo valor inteiro, ceil(4.01) = 5

k = 1:n;
fi2 = 1/n * asinh(1/E);
tk = (2*k-1) * pi/(2*n);
pk = -sinh(fi2)*sin(tk) + 1j*cosh(fi2)*cos(tk)

%% Projeto Filtro FIR - Janela Fixa

Ap = 0.5; % Ganho na banda de passagem em dB
As = 40; % Atenua��o no stopband em dB
fa = 4000; % Hz
f1 = 1000; % Hz
f2 = 1300; % Hz
GdB = 0; % dB

g0 = GdB;

f = [1000 1300]; % frequ�ncias em Hz

% substituindo de Hz para �mega
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
    wp1 = 0.4998*pi; ws1 = ws; % valores medidos no gr�fico
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
% % Terceiro ajuste da frequ�ncia de corte
%     Dwp = 0.4*pi - wp2;
%     Dws = 0.6*pi - ws2;
%     wc3 = wc + (Dwp + Dws)/2;
%     wc = wc3;


wc = sqrt(wp*ws);   % frequ�ncia de corte, m�dia das frequ�ncias

k = 1:M;

% Highpass
bi = -sin(k*wc)./(k*pi); 
b0 = 1 - ( wc/pi);
b = [flip(bi) b0 bi];

m = -M : M;
%mk = 0.04; % filtro hamming
mk = 0.00; % filtro hann
wk = (0.5+mk)+(0.5-mk)*cos(2*pi*m/(2*M+1)); %serve para filtro hamming ou hann
%wk = (0.42)+(0.5)*cos(2*pi*m/(2*M+1))+(0.08)*cos(4*pi*m/(2*M+1)); %filtro blackman
b = b.*wk*10^(-g0/20);
%b = b.*10^(-g0/20);     %janela retangular, basta tirar a janela


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
%subplot(2,2,[1 3])
zplane(b, 1)
axis([-2 2 -2 2]);

subplot(323)
plot(w/pi, unwrap(angle(h))/pi); grid on;title('Resposta de fase de H(z)');
subplot(325)
%%subplot(2,2,[2 4])
grpdelay(b, 1);title('Atraso de grupo');

figure(2)
hold on;
plot(w/pi,20*log10(abs(h))); grid on; xlim([0 1]);title('Resposta de magnitude de H(z)');ylim([-80 5]);
plot([0 wp/pi wp/pi 1],[-As -As 0 0], '--red')
plot([ws/pi,ws/pi,1],[-60 -Ap,-Ap], '--red')
hold off;
