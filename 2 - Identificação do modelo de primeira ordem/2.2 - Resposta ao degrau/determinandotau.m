% Nome do arquivo
arquivo = 'dados/oficial_.CSV';

% Ler os dados do arquivo CSV
dados = readtable(arquivo);

% Extrair as colunas
t = dados{:, 1}; % Coluna 't in s'
C1 = dados{:, 2};    % Coluna 'C1 in V'
C2 = dados{:, 3};    % Coluna 'C2 in V'

% Determinando y e u inicial
indice_y_zero = (t >= -0.1 & t < 0);
y_zero = mean(C2(indice_y_zero));
u_zero = mean(C1(indice_y_zero));

% Corrigindo o ponto inicial para passar em 0,0
for i = 1:size(C2)
    C2(i) = C2(i) - y_zero;
    C1(i) = C1(i) - u_zero;
end

% Definir o intervalo de t desejado 
ts = 0.4;
tf = 0.6;
indice = (t >= -0.001 & t <= ts);   % Regime transiente
indice_ss = (t >= ts & t <= tf);    %Regime permanente

% Fatiar os vetores com base no intervalo de t
t_ts = t(indice);
C1_ts = C1(indice);
C1_ss = C1(indice_ss);
C2_ts = C2(indice);
C2_ss = C2(indice_ss);

%%
% Calculando K
y_ss = mean(C2_ss);
A = mean(C1_ts);

K = y_ss / A;

Kt = 0.0159060;
Ka = K/Kt;

fprintf('K = %d \n', K)
fprintf('Kt = %d \n', Kt)
fprintf('Ka = %d \n\n', Ka)

% Calculando tau usando método da área
A0 = trapz(t_ts,y_ss - C2_ts);

tau_area = A0/(y_ss);

fprintf('tau_area = %d \n', tau_area);

% Calculando tau usando método do logaritmo neperiano

t_neperiano = t(t >= -0.001 & t <= 0.2);
C2_neperiano = C2(t >= -0.001 & t <= 0.2);
b = log(y_ss ./ (y_ss - C2_neperiano));

a = (t_neperiano' * b) / norm(t_neperiano)^2;

tau_neperiano = (1/a);

fprintf('tau_neperiano = %d \n', tau_neperiano);

%%
% Plotar C1_ts e C2_ts no mesmo gráfico
figure;
plot(t, C1, '-r', 'DisplayName', 'C1 in V'); % Plotar C1_ts em vermelho
hold on;
plot(t, C2, '-b', 'DisplayName', 'C2 in V'); % Plotar C2_ts em azul
hold off;

% Configurações do gráfico
xlabel('t (s)');
ylabel('Voltagem (V)');
title('Gráfico de C1 e C2 no mesmo gráfico');
legend('show');
grid on;

% Plotar C1_ts e C2_ts no mesmo gráfico
figure;
plot(t_ts, C1_ts, '-r', 'DisplayName', 'C1_{ts} in V'); % Plotar C1_ts em vermelho
hold on;
plot(t_ts, C2_ts, '-b', 'DisplayName', 'C2_{ts} in V'); % Plotar C2_ts em azul
hold off;

% Configurações do gráfico
xlabel('t_ts (s)');
ylabel('Voltagem (V)');
title('Gráfico de C1_ts e C2_ts no mesmo gráfico');
yline(y_ss, '--p', 'DisplayName', 'y_{ss} in V');
legend('show');
grid on;

%%
%Validação do Modelo (em relação ao K obtido no primeiro experimento)

K_exp1 = 1.380095;

Erro = ((K_exp1 - K ) / K_exp1) * 100;
fprintf('Erro percentual = %d \n', Erro)