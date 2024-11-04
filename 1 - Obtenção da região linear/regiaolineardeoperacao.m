%% Determinando Região Linear de Operação

% Dados de tensão de entrada (Va), saída (Vt) e velocidade angular (w)
Va = [0, 0.029064, 0.78054, 1.1635, 1.8808, 2.2112, 2.8303, 3.2567, 3.5928, 4.1681, 4.5053, 5.0858, 5.4667, 6.1924, 6.8298, 7.8074, 8.383, 9.320, 10.201, 11.186, 12.386, 13.150, 14.163, 15.418, 16.158, 17.517, 18.532, 20.234]';
Vt = [0, -0.0005908, -0.0043344, -0.0067431, -0.011596, -0.014576, 1.2643, 1.8703, 2.2837, 3.0337, 3.5945, 4.4143, 4.9795, 5.9835, 6.819, 8.2149, 8.9774, 10.404, 11.577, 13.006, 14.684, 15.594, 16.191, 18.692, 19.701, 21.517, 22.922, 25.378]';
w = [0, 0, 0, 0, 0, 33.2, 79.6, 116.5, 144.0, 199.1, 228.1, 274.5, 318.6, 373.1, 438.2, 516.7, 562.5, 652.4, 729.0, 818.0, 916.0, 983.0, 1066.0, 1183.0, 1249.0, 1252.0, 1449, 1612]';

% Definição do vetor b
b = Vt; 

% Definição da Matriz A

n = 6; % Grau do polinômio para ajuste dos dados

A = zeros(length(Va), n);

for i = 1:n
    A(:, i) = Va.^(n-i+1); % Preenche cada coluna com x^(n-i+1)
end

% Calculando vetor x
x = (A' * A) \ (A' * b);

% Fazendp x* = [x; 0] 
coeffs = [x; 0];    % Coeficientes do polinômio ajustado

newVt = polyval(coeffs, Va);    % Novos valores de Vt para o polinômio ajustado

% Definindo a derivada do polinômio ajustado
coeffs_dot = polyder(coeffs);
newVt_dot = polyval(coeffs_dot, Va);

%%%%%%%%%%%%%% Plote dos dados %%%%%%%%%%%%%%%%%%
plot(Va, [newVt, Vt], '-o'); 

xlabel('Va');
ylabel('Vt');
title('Gráfico de Vt em função de Va');

grid on;
figure

%%%%%%%%%%%%%% Plote dos dados %%%%%%%%%%%%%%%%%%
plot(Va, newVt_dot, '-o'); 

xlabel('Va');
ylabel('Vt_{dot}');
title('Gráfico da derivada do polinomio ajustado');

% Adiciona linhas delimitando região linear
xline(3.25, '--r', 'x = 3.25');
xline(20.234, '--r', 'x = 20.234');

grid on;
figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Definindo os limites inferior e superior
limite_inferior = 3.25;
limite_superior = 20.234;

% Filtrando os pontos dentro da região de estabilidade
Va_linear = Va(Va >= limite_inferior & Va <= limite_superior);
Vt_linear = Vt(Va >= limite_inferior & Va <= limite_superior);
newVt_linear = newVt(Va >= limite_inferior & Va <= limite_superior);
w_linear = w(Va >= limite_inferior & Va <= limite_superior);

% Corrigindo o ponto inicial para passar no 0,0
Va_linear_offset = Va_linear(1);
Vt_linear_offset = Vt_linear(1);
newVt_linear_offset = newVt_linear(1);
w_linear_offset = w_linear(1);

for i = 1:size(Vt_linear)
    Vt_linear(i) = Vt_linear(i) - Vt_linear_offset;
    newVt_linear(i) = newVt_linear(i) - newVt_linear_offset;
    Va_linear(i) = Va_linear(i) - Va_linear_offset;
    w_linear(i) = w_linear(i) - w_linear_offset;
end

%%%%%%%%%%%%%% Plote dos dados %%%%%%%%%%%%%%%%%%
plot(Va_linear, [newVt_linear, Vt_linear], '-o'); 

xlabel('Va_{linear}');
ylabel('Vt_{linear}');

title('Gráfico de Vt_{linear} em função de Va_{linear}');

grid on;
figure

fprintf('Região Linear de operação: %d < Va < %d \n\n', limite_inferior, limite_superior);
%% Definindo o parâmetro K

% Definição do vetor b
b_linear = Vt_linear;

% Definição da Matriz A

n = 1; % Grau do polinômio para ajuste dos dados (dessa vez, uma reta)

A = zeros(length(Va_linear), n);

for i = 1:n
    A(:, i) = Va_linear.^(n-i+1); % Preenche cada coluna com x^(n-i+1)
end

% Calculando o vetor x
x = (A' * A) \ (A' * b_linear);

% Fazendp x* = [x; 0] 
coeffs = [x; 0];    % Coeficiente da reta ajustada

newVt_linear = polyval(coeffs, Va_linear);  % Novos valores de Vt para o polinômio ajustado

% Definindo a derivada da reta ajustada
coeffs_dot = polyder(coeffs);
newVt_linear_dot = polyval(coeffs_dot, Va_linear);

% Comparando os valores de K obtidos através da média dos pontos da
% derivada e do coeficiente da reta ajustada utilizando os pontos da
% região linear de operação
fprintf('K = %d (Média dos valores da região linear) \n', mean(newVt_linear_dot));
fprintf('K = %d (Coeficiente da reta Vt/Va) \n', coeffs_dot);

K = mean(newVt_linear_dot);

%% Definindo o parâmetro Kt e Ka

% Realizando o mesmo procedimento para w_linear, aproximando agora por uma reta
% para descobrir Kt, uma vez que Vt_linear = Kt * w_linear

% Definição do vetor b
b_linear = Vt_linear;

% Definição da Matriz A

n = 1; % Grau do polinômio para ajuste dos dados (dessa vez, uma reta)

A = zeros(length(w_linear), n);

for i = 1:n
    A(:, i) = w_linear.^(n-i+1); % Preenche cada coluna com x^(n-i+1)
end

% Calculando vetor x
x = (A' * A) \ (A' * b_linear);

% Fazendp x* = [x; 0]
coeffs = [x; 0];    % Coeficiente da reta ajustada

newVt = polyval(coeffs, w_linear);   % Novos valores de Vt para o polinômio ajustado

% Definindo a derivada do polinômio ajustado
coeffs_dot = polyder(coeffs);
newVt_dot = polyval(coeffs_dot, w_linear);

Kt = coeffs_dot;
Ka = K/Kt;

%%%%%%%%%%%%%% Plote dos dados %%%%%%%%%%%%%%%%%%
plot(w_linear, [newVt, Vt_linear], '-o'); 

xlabel('w');
ylabel('Vt');
title('Gráfico de Vt em função de w (ajustado por uma reta)');

grid on;
figure

%%%%%%%%%%%%%% Plote dos dados %%%%%%%%%%%%%%%%%%
plot(w_linear, newVt_dot, '-o'); 

xlabel('w');
ylabel('Vt_{dot}');
title('Gráfico da derivada do polinomio ajustado');

grid on;

fprintf('Kt = %d V/rpm \n', Kt)
fprintf('Ka = %d rpm/V \n\nc', Ka)

%%
%Validação do Modelo

KaKt = Ka * Kt;

E = (((Ka * Kt) - K ) / (Ka * Kt)) * 100;

fprintf('Ka * Kt = %d \n', KaKt)
fprintf('Erro percentual = %d \n', E)