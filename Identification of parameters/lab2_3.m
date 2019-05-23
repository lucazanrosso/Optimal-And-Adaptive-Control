%% Esercizio 2
clear; close all;

N = 100; %dati totali
x = linspace(0,1,N); %(1 x n) suddivide lo spazio da 0 a 1 in 100 punti
sigma = 1;
ni = randn(N,sigma)*0.2;
y_real = sin(4*x) + 1./(4*x + 1);
y = y_real + ni';

M = 20;
tr = 80;    % numero dati training
vl = N - tr;    % numero dati validation
index_rand = randperm(N); %genera un vettore di 100 numeri a caso
index_tr = index_rand(1:tr); %indici dati di training 
index_vl = index_rand(tr+1:N); %indici dati di validation
 
x_tr = zeros(tr,1); %istanzio un arrray per i dati di training
x_vl = zeros(vl,1); %istanzio un array per i dati di validation
y_tr = zeros(tr,1);
y_vl = zeros(vl,1);
sigma_hat_tr = zeros(M,1);
sigma_hat_vl = zeros(M,1);

for i = 1 : tr           %creo un ciclo for per creare un array di dati di training   
    j = index_tr(i);
    x_tr(i) = x(j);
    y_tr(i) = y(j);    
end

for i = 1 : vl          %creo un ciclo for per creare un array di dati di validation 
    j = index_vl(i);
    x_vl(i) = x(j);
    y_vl(i) = y(j);   
end

for i = 1 : M
    PHI_tr = x_tr.^(0:i);
    PHI_vl = x_vl.^(0:i);
    
    theta = PHI_tr \ y_tr; %calcolo i coefficienti ottimi per il polinomio ((m + 1) x 1)
    
    y_hat_tr = PHI_tr * theta; %calcolo la funzione ottenuta dal polinomio stimatore
    y_hat_vl = PHI_vl * theta; %uso i theta calcolato con il training

    sigma_hat_tr(i) = norm (y_tr - y_hat_tr) / sqrt(tr);  %sigma training
    sigma_hat_vl(i) = norm (y_vl - y_hat_vl) / sqrt(vl);  %sigma validation
end

% Scelgo la soluzione con m = 3
m = 3;
PHI_tr = x_tr.^(0:m);
theta = PHI_tr \ y_tr;
PHI = x'.^(0:m);
y_hat = PHI * theta;
upper_bound = y_hat + 3*sigma_hat_vl(m);
lower_bound = y_hat - 3*sigma_hat_vl(m);

%Grafico y
figure(1);
plot(x, y, '.'); %blu
hold on;
plot(x_tr, y_hat_tr, '.'); %nero
hold on;
plot(x_vl, y_hat_vl, '.'); %nero
grid on;
title('Esercizio 2 MODEL SECTION Data')
xlabel('m')
ylabel('x')
legend('y_{vero}','y_{training}','y_{validation}')
hold off;

%Grafici sigma
figure(2);
plot((1:M),sigma_hat_tr);
hold on;
plot((1:M),sigma_hat_vl);
hold on;
grid on;
title('Esercizio 2 MODEL SECTION \sigma')
xlabel('m')
ylabel('\sigma')
legend('\sigma_{tr}','\sigma_{vl}')
hold off;

%Grafico m = 3
figure(3)
plot(x, y_real); % blu
hold on;
plot(x, y_hat); % blu
hold on;
plot(x, upper_bound, 'g'); % rosso
hold on;
plot(x, lower_bound, 'g'); % rosa
hold off;
grid on;
title('grafico RLS (known noise variance)')
xlabel('numero misurazioni')
ylabel('teta stim e \sigma_n')
legend('y_{real}','y_{hat}','upper bound','lower bound')
hold off