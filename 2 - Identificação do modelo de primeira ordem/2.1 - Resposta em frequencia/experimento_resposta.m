clear; clc;

freqs = logspace(log10(0.189), log10(18.9), 10)

gjwdb = zeros(1, length(freqs));
phi = zeros(1, length(freqs));

for i = 1:10

    filename = sprintf('dados/FREQ0%d.CSV', i);
    data = readtable(filename);
 
    % Frequência da amostra atual
    f_k = freqs(i);
    T = 1/f_k; % Período da frequência
    
    % Ajuste de tempo, começando de zero
    k = (data.('inS') - data.('inS')(1800));  % Subtrai o valor inicial do ponto de índice 1800
    
    % Filtrar os valores de tempo que são menores ou iguais a T, a partir do índice 1800
    t = k(1800:end);  % Começa do índice 1800 até o final
    t = t(t <= T);  % Filtra para pegar apenas um período
    
    % Filtrar os sinais de entrada e saída com base no tamanho do vetor de tempo
    saida = data.('C2InV')(1800:1800+length(t)-1);  % Ajusta o sinal de saída para o mesmo tamanho do tempo filtrado
    entrada = data.('C1InV')(1800:1800+length(t)-1);  % Ajusta o sinal de entrada para o mesmo tamanho do tempo filtrado

    
    % Cálculo para a entrada
    a0_in_k = trapz(t, entrada) / T;
    an_in_k = (2/T) * trapz(t, cos(2*pi*f_k*t) .* entrada);
    bn_in_k = (2/T) * trapz(t, sin(2*pi*f_k*t) .* entrada);
    zn_in_k = bn_in_k + 1i * an_in_k;
    cn_in_k = abs(zn_in_k);
    phase_in_k = angle(zn_in_k);
    
    
    % Cálculo para a saída (ajuste correto)
    a0_out_k = trapz(t, saida) / T;
    an_out_k = (2/T) * trapz(t, cos(2*pi*f_k*t) .* saida);
    bn_out_k = (2/T) * trapz(t, sin(2*pi*f_k*t) .* saida);
    zn_out_k = bn_out_k + 1i * an_out_k; % Corrigido para saída
    cn_out_k = abs(zn_out_k);
    phase_out_k = angle(zn_out_k);
    
    % Cálculo do ganho (em dB) e defasagem (em graus)
    gjwdb(i) = 20 * log10(abs(cn_out_k) / abs(cn_in_k));

    if (phase_out_k - phase_in_k) * (180/pi) < 0 % Convertendo para graus
        phi(i) = (phase_out_k - phase_in_k) * (180/pi);
    else
        phi(i) = (phase_out_k - phase_in_k) * (180/pi) - 360;
    end
end

freqs_wn = 2 * pi * freqs;

% Plotar o diagrama de Bode
figure;
subplot(2, 1, 1);
semilogx(freqs, gjwdb, '-o');
xlabel('Frequência (rad/s)');
ylabel('Ganho (dB)');
title('Diagrama de Bode - Ganho');
grid on;

subplot(2, 1, 2);
semilogx(freqs, phi, '-o');
xlabel('Frequência (rad/s)');
ylabel('Defasagem (graus)');
title('Diagrama de Bode - Fase');
grid on;

% ##### Método dos mínimos quadrados (3 primeiros pontos) ########
figure

first_freqs_points = freqs(1,1:3);
first_gjwdb_points = gjwdb(1,1:3);

[coeffs, adjusted_gjwdb_points] = minimos_quadrados(first_freqs_points, first_gjwdb_points, 8);


% ##### Método dos mínimos quadrados (todos os pontos) #############

[coeffs, adjusted_gjwdb_points] = minimos_quadrados(freqs, gjwdb, 8);

ww = logspace(log10(0.189), log10(18.9), 200);

g = polyval(coeffs, ww);

% deslocando a curva para passar na origem
    y_zero = adjusted_gjwdb_points(1);
    for i = 1:size(g, 2)
        g(i) = g(i) + y_zero;
    end

KdB = mean(adjusted_gjwdb_points(1,1:2))

semilogx(freqs_wn, [gjwdb; adjusted_gjwdb_points], '-o'); 
hold on

legend('real', 'ajustado')
% semilogx(ww*2*pi, g, '--');

yline(KdB, '--r', sprintf('y = %.2f', KdB));
yline((KdB - 3), '--r', sprintf('y = %.2f', KdB - 3));
% xline ((8.82), '--r', 'y = 8.82');
% xline ((1.4037), '--r', 'y = 1.4037');

K = 10^(KdB/20)

[~, idx] = min(abs(adjusted_gjwdb_points - (KdB-3)));
freq_corte = freqs_wn(idx)
tau = 1/freq_corte

xline(freq_corte, '--r', sprintf('x = %.2f', freq_corte));

xlabel('w');
ylabel('dB');

title('Ajuste por aproximação linear de gjwdb');

grid on;
hold on


% Criando a função de transferência
num = K;
den = [tau 1];
sys = tf(num, den);

% Gerando o Bode plot
[mag, phase, w] = bode(sys, {freqs_wn(1), freqs_wn(10)});
mag_db = 20 * log10(squeeze(mag));
freq_hz = squeeze(w) / (2 * pi);
phase = squeeze(phase);

% Figura para o Bode Plot Modificado
figure;

% Primeiro subplot - Magnitude (Ganho em dB)
subplot(2, 1, 1);
semilogx(freq_hz, mag_db, '-');       
hold on;
semilogx(freqs, gjwdb, 'o');          
hold off;
grid on;
title('Diagrama de Bode - Ganho');
xlabel('Frequência (rad/s)');
ylabel('Magnitude / Ganho (dB)');
legend('modelo', 'real', 'Location', 'best');

% Segundo subplot - Fase
subplot(2, 1, 2);
semilogx(freq_hz, phase, '-');        
hold on;
semilogx(freqs, phi, 'o');           
hold off;
grid on;
title('Diagrama de Bode - Fase');
xlabel('Frequência (rad/s)');
ylabel('Fase (graus)');
legend('modelo', 'real', 'Location', 'best');