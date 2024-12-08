
k = 1.384063;
tau = 0.08420596;

kp = 1/k; % para erro em regime permanente nulo, kp = 1/k (não sei o porquê)
kt = 0.015906;
ka = k/kt;

% calculo de ki para criticamente amortecido
ki_crit = 1/(4*k*tau); % criticamente amortecido


% calculo de ki para subamortecido com overshoot de 5%
overshoot = 0.05;
zeta = abs(log(overshoot))/(sqrt(pi^2+(log(overshoot)^2)));
wn = 1/(2*tau*zeta);

ki_sub = (wn^2*tau)/(ka*kt);
ts_ki = 4/(zeta*wn); % tempo de acomodação malha fechada
ts_kp = 4*tau; % tempo de acomodação para primeira ordem malha aberta

z = (1/tau) + 5; % selecionando z mais à esquerda que 1/tau por 5 unidades
Ti = 1/z;

wn_ti = sqrt((k*kp)/(tau*Ti));
zeta_ti = ((1+k*kp)/(2*tau)) * (sqrt((tau*Ti)/(k*kp)));

ts_ti = 4/(wn_ti*zeta_ti);

% Plotando root locus sistema
numerador = [12.5 1];     
denominador = [12.5 0]; 
G = tf(numerador, denominador);
rlocus(G);
title('Lugar das Raízes');
xlabel('Parte Real');
ylabel('Parte Imaginária');
