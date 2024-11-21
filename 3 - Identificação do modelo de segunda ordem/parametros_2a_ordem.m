%clear all

filename = 'dados/dadoslinear.CSV'; % Substitua pelo nome do seu arquivo
data = readtable(filename);

data.Properties.VariableNames = {'Tempo', 'Va', 'Vt', 'Ia'};

% t = data.Tempo;
% Va = data.Va;
% Vt = data.Vt;
% Ia = 20* data.Ia;

ganho_corrente = 10;

t = data.Tempo(2:2:end);
Va = (data.Va(2:2:end)+data.Va(1:2:end))/2;
Vt = (data.Vt(2:2:end)+data.Vt(1:2:end))/2;
Ia = (ganho_corrente *data.Ia(2:2:end)+ ganho_corrente*data.Ia(1:2:end))/2;

Va = Va - Va(1);
Vt = Vt - Vt(1);
Ia = Ia - Ia(1);

% t = data.Tempo(1:2:end);
% Va = data.Va(1:2:end);
% Vt = data.Vt(1:2:end);
% Ia = 20* data.Ia(1:2:end);



kt = 0.0159;
k = 1.380;
ka = k/kt;

kg = kt/k;


ue = Va - (kg/kt)*Vt;

um = ka*kt*Ia;

ue_0 = ue(1:end-1);
Ia_0 = Ia(1:end-1);
Ia_n = Ia(2:end);

me = [Ia_0, ue_0];

um_0 = um(1:end-1);
Vt_0 = Vt(1:end-1);
Vt_n = Vt(2:end);

mm = [Vt_0, um_0];


xe = me\Ia_n;

xm = mm\Vt_n;

Ra = (1 - xe(1))/xe(2);
La = - (Ra*(t(2)-t(1)))/log(xe(1));
f = (1-xm(1))/xm(2);
j = -(f*(t(2)-t(1)))/log(xm(1));


A = [ -Ra/La -kg/(kt*La); (ka*kt)/j -f/j];
B = [1/La ; 0];
C1 = [1 0];
C2 = [0 1];
D = 0;

simin_va.time = t;
simin_va.signals.values = Va;
simin_va.signals.dimensions = 1;

simin_vt.time = t;
simin_vt.signals.values = Vt;
simin_vt.signals.dimensions = 1;

simin_ia.time = t;
simin_ia.signals.values = Ia;
simin_a.signals.dimensions = 1;


%Funcao transferencia Vt/Va
[n,d]=ss2tf(A,B,C2,D)







