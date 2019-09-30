%% 4 Identification of discrete time linear systems

%% 4.1 Parameter identification
n = 100;                      % number of samples
u = data.signals.values(:,1);   % input samples
y_f = data.signals.values(:,2); % output samples
t = data.time;
time = t(1:n+1);

I = eye(3);
PHI_0 =  [ -y_f(1) u(2) u(1); -y_f(2) u(3) u(2); -y_f(3) u(4) u(3);];
Y_0 = [y_f(2); y_f(3); y_f(4)];
P = inv(PHI_0'*PHI_0);

theta_hat = P*PHI_0'*Y_0;
a0_hat = zeros(n+1, 1);
sigma_a0_hat = zeros(n+1, 1);

PHI = zeros(n+1, 3);
for i = 5 : n+1    
    phi =  [ -y_f(i-1); u(i); u(i-1) ];
    K = P * phi / (1 + phi'*P*phi);
    P = (I - K*phi')*P;
    theta_hat = theta_hat + K*(y_f(i) - phi'*theta_hat);     
    a0_hat(i) = theta_hat(1);
    
    PHI(i,:) = phi';
    sum = 0;
    for j = 5:i
        sum = sum + (y_f(j) - PHI(j,:)*theta_hat)^2;
    end
    sigma_a0_hat(i) = sqrt(P(1,1) * sum / (i - 4)); % parto da 5, quindi devo dividere per i - 4
end

b_real = ones(n+1, 1) * b;
b_hat = -(J / Ts) * log(-a0_hat);
b_upper_bound = -(J / Ts) * log(-a0_hat + 3 * sigma_a0_hat);
b_lower_bound = -(J / Ts) * log(-a0_hat - 3 * sigma_a0_hat);

figure(1)
N = (1:n+1);
p1 = plot(time, b_real, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 5);
hold on;
p2 = plot(time, b_hat, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 5);
hold on;
p3 = plot(time, b_upper_bound, 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 5);
hold on;
plot(time, b_lower_bound, 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 5);
hold off;
grid on;
xlabel({'t [s]'}, 'Interpreter', 'latex')
ylabel({'b [Nms]'}, 'Interpreter', 'latex')
legend([p1 p2 p3], {'$b$', '$\hat{b}$', '$\hat{b} \pm \hat{b}_{bounds}$'}, 'Interpreter', 'latex', 'FontSize', 36);
set(gca,'FontSize',24)

b_hat = b_hat(n+1)


%% Verifica b_hat
% Z = zeros(3);
% z = zeros(3,1);
% for i = 2 : n+1
%    phi =  [-y_f(i-1); u(i); u(i-1)];
%    Z = Z + phi*phi'; 
%    z = z + phi*y_f(i); 
% end
% theta_hat = inv(Z)*z;
% a_hat = theta_hat(1);
% b_hat = -(J / Ts) * log(- a_hat)
