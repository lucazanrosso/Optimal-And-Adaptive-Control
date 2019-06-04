%% Esercizio 2
clear; close all;

N = 100; %dati totali
x = linspace(0,1,N);                % generates N points. The spacing between the points is (1-0)/(N-1)
sigma = 1;
ni = randn(N,sigma)*0.2;            % noise
y_real = sin(4*x) + 1./(4*x + 1);   % real function
y = y_real + ni';                   % measured function

M = 20;                             % model complexity
tr = 80;                            % number of training data
vl = N - tr;                        % number of validation data
index_rand = randperm(N);           % returns a row vector containing a random permutation of the integers from 1 to n inclusive
index_tr = index_rand(1:tr);        % random index of training data
index_vl = index_rand(tr+1:N);      % random index of validation data 
 
x_tr = zeros(tr,1);                 % training data vector
x_vl = zeros(vl,1);                 % validation data vector
y_tr = zeros(tr,1);
y_vl = zeros(vl,1);
sigma_hat_tr = zeros(M,1);          % average prediction error vector
sigma_hat_vl = zeros(M,1);

for i = 1 : tr                      % fill the vectors with their data   
    j = index_tr(i);
    x_tr(i) = x(j);
    y_tr(i) = y(j);    
end

for i = 1 : vl
    j = index_vl(i);
    x_vl(i) = x(j);
    y_vl(i) = y(j);   
end

for i = 1 : M
    PHI_tr = x_tr.^(0:i);           % PHI training
    PHI_vl = x_vl.^(0:i);           % PHI training
    
    theta = PHI_tr \ y_tr;          % theta
    
    y_hat_tr = PHI_tr * theta;      % y_hat training
    y_hat_vl = PHI_vl * theta;      % y_hat validation

    sigma_hat_tr(i) = norm (y_tr - y_hat_tr) / sqrt(tr);    % sigma_hat training
    sigma_hat_vl(i) = norm (y_vl - y_hat_vl) / sqrt(vl);    % signa_hat validation
end

% I choose the solution that on average has a lower sigma (m = 3)
m = 3;
PHI_tr = x_tr.^(0:m);
theta = PHI_tr \ y_tr;
PHI = x'.^(0:m);
y_hat = PHI * theta;
upper_bound = y_hat + 3*sigma_hat_vl(m);
lower_bound = y_hat - 3*sigma_hat_vl(m);

% y plot (m = 20)
figure(1);
plot(x, y, '.');
hold on;
plot(x_tr, y_hat_tr, '.');
hold on;
plot(x_vl, y_hat_vl, '.');
grid on;
title('Estimation of y with m = 20')
xlabel('x')
ylabel('y')
legend('y_{true}','y_{training}','y_{validation}')
hold off;

% sigma plot
figure(2);
plot((1:M),sigma_hat_tr);
hold on;
plot((1:M),sigma_hat_vl);
hold on;
grid on;
title('sigma training and validation trends')
xlabel('m')
ylabel('\sigma')
legend('\sigma_{tr}','\sigma_{vl}')
hold off;

%y plot (m = 3)
figure(3)
plot(x, y_real);
hold on;
plot(x, y_hat);
hold on;
plot(x, upper_bound, 'g'); % green
hold on;
plot(x, lower_bound, 'g'); % green
hold off;
grid on;
title('y trend for the optimal m based on the validation set')
xlabel('x')
ylabel('y')
legend('y_{real}','y_{hat}','upper bound','lower bound')
hold off