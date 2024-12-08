% Carrega os dados
load("dados/control.mat");
load("dados/reference.mat");
load("dados/tacometer.mat");

% Extrai os sinais e o tempo
t = control.time;
control_sig = reshape(control.signals.values, [], 1);
reference_sig = reshape(reference.signals.values, [], 1);
tacometer_sig = reshape(tacometer.signals.values, [], 1);

% Define os índices inicial e final para fatiar os dados
start_idx = 800; % Alterar conforme necessário
end_idx = 1300;   % Alterar conforme necessário

% Fatia os dados
new_t = t(start_idx:end_idx);
new_control_sig = control_sig(start_idx:end_idx);
new_reference_sig = reference_sig(start_idx:end_idx);
new_tacometer_sig = tacometer_sig(start_idx:end_idx);

% Plota e salva o gráfico do sinal de controle
figure;
plot(new_t, new_control_sig, 'DisplayName', "Control signal");
xlabel('Time (s)');
ylabel('Amplitude');
title('Control Signal (Sliced)');
grid on;
legend;
saveas(gcf, 'control_signal_sliced.png'); % Salva como imagem

% Plota e salva o gráfico do sinal de referência
figure;
plot(new_t, new_reference_sig, 'DisplayName', "Reference signal");
xlabel('Time (s)');
ylabel('Amplitude');
title('Reference Signal (Sliced)');
grid on;
legend;
saveas(gcf, 'reference_signal_sliced.png'); % Salva como imagem

% Plota e salva o gráfico do sinal do tacômetro
figure;
plot(new_t, new_tacometer_sig, 'DisplayName', "Tacometer signal");
xlabel('Time (s)');
ylabel('Amplitude');
title('Tacometer Signal (Sliced)');
grid on;
legend;
saveas(gcf, 'tacometer_signal_sliced.png'); % Salva como imagem

% Plota e salva o gráfico com todos os sinais fatiados
figure;
plot(new_t, new_reference_sig, 'DisplayName', 'Reference signal'); hold on;
plot(new_t, new_control_sig, 'DisplayName', 'Control signal');
plot(new_t, new_tacometer_sig, 'DisplayName', 'Tacometer signal');
xlabel('Time (s)');
ylabel('Amplitude');
title('All Signals (Sliced)');
grid on;
legend; % 'DisplayName' preenche a legenda automaticamente
saveas(gcf, 'all_signals_sliced.png'); % Salva o gráfico como imagem
