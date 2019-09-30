%% Esercizio 2
% clear; close all;

a = 1;
A_tr = zeros(20, a);
A_vl = zeros(20, a);
for k = 1 : a

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
    PHI_vl = x_vl.^(0:i);           % PHI validation
    
    theta = PHI_tr \ y_tr;          % theta
    
    y_hat_tr = PHI_tr * theta;      % y_hat training
    y_hat_vl = PHI_vl * theta;      % y_hat validation

    sigma_hat_tr(i) = norm (y_tr - y_hat_tr) / sqrt(tr);    % sigma_hat training
    sigma_hat_vl(i) = norm (y_vl - y_hat_vl) / sqrt(vl);    % signa_hat validation
end

    A_tr(:, k) = [sigma_hat_tr];
    A_vl(:, k) = [sigma_hat_vl];
end

sigma_tr_av = zeros(M, 1);
sigma_vl_av = zeros(M, 1);
for i = 1 : M
    sigma_tr_av(i) = mean(A_tr(i, :));
    sigma_vl_av(i) = mean(A_vl(i, :)) ;
end

% I choose the solution that on average has a lower sigma (m = 3)
m = 4;
PHI_tr = x_tr.^(0:m);
theta = PHI_tr \ y_tr;
PHI = x'.^(0:m);
y_hat = PHI * theta;
upper_bound = y_hat + 3*sigma_hat_vl(m);
lower_bound = y_hat - 3*sigma_hat_vl(m);

% y plot (m = 20)
% figure(1);
% plot(x, y, '.');
% hold on;
% plot(x_tr, y_hat_tr, '.');
% hold on;
% plot(x_vl, y_hat_vl, '.');
% grid on;
% title('Estimation of y with m = 20')
% xlabel('x')
% ylabel('y')
% legend('y_{true}','y_{training}','y_{validation}')
% hold off;
figure(1)
p1 = plot(x, y_real, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 5);
hold on;
% p2 = plot(x, y_hat, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 5);
% hold on;
p3 = plot(x_tr, y_tr, '.', 'MarkerSize', 20, 'MarkerEdgeColor', [0.9290, 0.6940, 0.1250]);
hold on;
p4 = plot(x_vl, y_vl, '.', 'MarkerSize', 20, 'MarkerEdgeColor', [0.4940, 0.1840, 0.5560]);
hold off;
grid on;
xlabel({'x'}, 'Interpreter', 'latex')
ylabel({'y'}, 'Interpreter', 'latex')
legend([p1 p3 p4], {'$y$ (real model)', '$\hat{y}$ (estimated model)', '$y_{tr}$ (training data)', '$y_{vl}$ (validation data)'}, 'Interpreter', 'latex', 'FontSize', 36, 'Location', 'Southwest');
set(gca,'FontSize',24)

% sigma plot
% figure(2);
% p5 = plot((1:M), sigma_tr_av, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 5);
% hold on;
% p6 = plot((1:M), sigma_vl_av, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 5);
% hold on;
% grid on;
% xlabel({'$m$ number of parameters'}, 'Interpreter', 'latex')
% ylabel({'$\hat{\sigma}$'}, 'Interpreter', 'latex')
% legend([p5 p6], {'$\hat{\sigma}^{tr}_y$ (training)', '$\hat{\sigma}^{vl}_y$ (validation)'}, 'Interpreter', 'latex', 'FontSize', 36, 'Location', 'North');
% set(gca,'FontSize',24)
% hold off;

%y plot (m = 3)
% figure(3)
% plot(x, y_real);
% hold on;
% plot(x, y_hat);
% hold on;
% plot(x, upper_bound, 'g'); % green
% hold on;
% plot(x, lower_bound, 'g'); % green
% hold off;
% grid on;
% title('y trend for the optimal m based on the validation set')
% xlabel('x')
% ylabel('y')
% legend('y_{real}','y_{hat}','upper bound','lower bound')
% hold off

figure(3)
p7 = plot(x, y_real, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 5);
hold on;
p8 = plot(x, y_hat, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 5);
hold on;
p9 = plot(x, upper_bound, 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 5);
hold on;
plot(x, lower_bound, 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 5);
hold off;
grid on;
xlabel({'x'}, 'Interpreter', 'latex')
ylabel({'y'}, 'Interpreter', 'latex')
legend([p7 p8 p9], {'$y$ (real model)', '$\hat{y}$ (estimated model)', '$\hat{y} \pm 3\hat{\sigma}^{vl}_y$'}, 'Interpreter', 'latex', 'FontSize', 36, 'Location', 'Southwest');
set(gca,'FontSize',24)