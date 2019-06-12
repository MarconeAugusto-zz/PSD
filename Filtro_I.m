%% Dados
%   LP - (fa = 4000 Hz, f1 = 1000 Hz; f2 = 1300 Hz, Ap = 2 dB, As = 30 dB, GdB = 5 dB)
%   IIR - Butterworth, FIR - Janela Ajustável

% Projeto Filtro IIR - Butterworth

Ap = 2; % Ganho na banda de passagem em dB
As = 30; % Atenuação no stopband em dB
fa = 4000; % Hz
f1 = 1000; % Hz
f2 = 1300; % Hz
GdB = 5; % dB

% Projeto Filtro FIR - Janela Ajustável

clear all;
close all;
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



