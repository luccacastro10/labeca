clear; clc;

freqs = logspace(log10(0.189), log10(18.9), 10);

gjwdb = zeros(1, length(freqs));
phi = zeros(1, length(freqs));

for i = 1:10

    filename = sprintf('dados/FREQ0%d.CSV', i);
    data = readtable(filename);
 
    % Frequência da amostra atual
    f_k = freqs(i);
    T = 1/f_k; % Período da frequência
    
    % Ajuste de tempo, começando de zero
    k = (data.('inS') - data.('inS')(1));  % Subtrai o valor inicial do ponto de índice 1800
    
    % Filtrar os valores de tempo que são menores ou iguais a T, a partir do índice 1800
    t = k(1:end);  % Começa do índice 1800 até o final
    t = t(t <= T);  % Filtra para pegar apenas um período
    
    % Filtrar os sinais de entrada e saída com base no tamanho do vetor de tempo
    saida = data.('C2InV')(1:1+length(t)-1);  % Ajusta o sinal de saída para o mesmo tamanho do tempo filtrado
    entrada = data.('C1InV')(1:1+length(t)-1);  % Ajusta o sinal de entrada para o mesmo tamanho do tempo filtrado

    
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

% Plotar o diagrama de Bode
figure;
subplot(2, 1, 1);
semilogx(freqs, gjwdb, '-o');
xlabel('Frequência (Hz)');
ylabel('Ganho (dB)');
title('Diagrama de Bode - Ganho');
grid on;

subplot(2, 1, 2);
semilogx(freqs, phi, '-o');
xlabel('Frequência (Hz)');
ylabel('Defasagem (graus)');
title('Diagrama de Bode - Fase');
grid on;

% ##### Método dos mínimos quadrados (3 primeiros pontos) #############
figure

first_freqs_points = freqs(1,1:3);
first_gjwdb_points = gjwdb(1,1:3);

[coeffs, adjusted_gjwdb_points] = minimos_quadrados(first_freqs_points, first_gjwdb_points, 1);

plot(first_freqs_points, [first_gjwdb_points; adjusted_gjwdb_points], '-o'); 

xlabel('w');
ylabel('db');

title('Ajuste por aproximação linear de gjwdb');

grid on;

% ##### Método dos mínimos quadrados (todos os pontos) #############
figure

[coeffs, adjusted_gjwdb_points] = minimos_quadrados(freqs, gjwdb, 3);

plot(freqs, [gjwdb; adjusted_gjwdb_points], '-o'); 

xlabel('w');
ylabel('db');

title('Ajuste por aproximação linear de gjwdb');

grid on;
