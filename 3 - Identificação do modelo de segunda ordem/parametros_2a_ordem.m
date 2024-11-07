clear, clc
% Determinando Região Linear de Operação

% Dados de tensão de entrada (Va), saída (Vt) e velocidade angular (w)
Va = [0, 0.029064, 0.78054, 1.1635, 1.8808, 2.2112, 2.8303, 3.2567, 3.5928, 4.1681, 4.5053, 5.0858, 5.4667, 6.1924, 6.8298, 7.8074, 8.383, 9.320, 10.201, 11.186, 12.386, 13.150, 14.163, 15.418, 16.158, 17.517, 18.532, 20.234]';
Vt = [0, -0.0005908, -0.0043344, -0.0067431, -0.011596, -0.014576, 1.2643, 1.8703, 2.2837, 3.0337, 3.5945, 4.4143, 4.9795, 5.9835, 6.819, 8.2149, 8.9774, 10.404, 11.577, 13.006, 14.684, 15.594, 16.191, 18.692, 19.701, 21.517, 22.922, 25.378]';
W = [0, 0, 0, 0, 0, 33.2, 79.6, 116.5, 144.0, 199.1, 228.1, 274.5, 318.6, 373.1, 438.2, 516.7, 562.5, 652.4, 729.0, 818.0, 916.0, 983.0, 1066.0, 1183.0, 1249.0, 1252.0, 1449, 1612]';


% Filtrando os pontos dentro da região de estabilidade
limite_inferior = 3.25;
limite_superior = 20.234;

Va_linear = Va(Va >= limite_inferior & Va <= limite_superior);
Vt_linear = Vt(Va >= limite_inferior & Va <= limite_superior);
W_linear = W(Va >= limite_inferior & Va <= limite_superior);

% Corrigindo o ponto inicial para passar no 0,0
Va_linear_offset = Va_linear(1);
Vt_linear_offset = Vt_linear(1);
W_linear_offset = W_linear(1);

for i = 1:size(Vt_linear)
    Vt_linear(i) = Vt_linear(i) - Vt_linear_offset;
    Va_linear(i) = Va_linear(i) - Va_linear_offset;
    W_linear(i) = W_linear(i) - W_linear_offset;
end

K_barra = ((Va_linear)'*(Vt_linear))/((Va_linear)'*(Va_linear))
K_t = ((W_linear)'*(Vt_linear))/((W_linear)'*(W_linear))
K_g = K_t/K_barra
K_a = K_g

% Carregando os dados do experimento de identificacao dos parametros de 2a
% ordem:

% Nome do arquivo
arquivo = 'dados/dadoslinear.CSV';

% Ler os dados do arquivo CSV
dados = readtable(arquivo);
ganho_sensor_corrente = 20;

% Extrair as colunas
t = dados{2:2:end, 1}; % Coluna 't in s'
Va = dados{2:2:end, 2};    % Coluna 'C1 in V'
Vt = dados{2:2:end, 3};    % Coluna 'C2 in V'
Ia = ganho_sensor_corrente*dados{2:2:end, 4};    % Coluna 'C4 in V'
h = t(2) - t(1);


ue = zeros(size(t));
um = zeros(size(t));

for tk_i = 1:size(t)
    ue(tk_i) = Va(tk_i)-(K_g/K_t)*Vt(tk_i);
    um(tk_i) = K_a*K_t*Ia(tk_i);
end

Me = [Ia(1:end-1), ue(1:end-1)];
Mm = [Vt(1:end-1), um(1:end-1)];


ia_barra = Ia(2:end);
vt_barra = Vt(2:end);

xe = inv(Me'*Me)*Me'*ia_barra
xm = inv(Mm'*Mm)*Mm'*vt_barra

phi_e = xe(1) 
phi_m = xm(1)
gama_e = xe(2)
gama_m = xm(2)

Ra = (1 - phi_e)/gama_e
La = - (Ra*h)/log(phi_e)
f = (1 - phi_m) / gama_m
J = -(f*h) / log(phi_m)

Km = 1/f
Km2 = gama_m/(1-phi_m)





