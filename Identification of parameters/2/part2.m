%% 2 Estimation of a scalar parameter

%% 2.1 Recursive least squares (known noise variance)
n = 1000;
sigma = 1;
ni = randn(n,sigma);
theta_real = ones(n, 1);
y = theta_real + ni;

theta_hat = 0;
theta_hat_vect = zeros(n, 1);
sigma_hat_vect = zeros(n, 1);

for i= 1:n
    theta_hat = ((i-1)/(i))*theta_hat + (1/(i))*y(i);
    theta_hat_vect(i) = theta_hat;
    sigma_hat = sigma/sqrt(i);
    sigma_hat_vect(i) = sigma_hat;
end

upper_bound = theta_hat_vect + 3*sigma_hat_vect;
lower_bound = theta_hat_vect - 3*sigma_hat_vect;
rn = sigma_hat_vect ./ theta_hat_vect;

%% 2.2 Recursive least squares (unknown noise variance)

theta_hat = 0;
theta_hat_vect = zeros(n, 1);
sigma_hat_vect = zeros(n, 1);

for i = 1:n
    theta_hat = ((i-1)/(i))*theta_hat + (1/(i))*y(i);
    theta_hat_vect(i) = theta_hat;
    sum = 0;
    for j = 1:i
        sum = sum + (y(j) - theta_hat)^2;
    end
    sigma_hat = sqrt(sum / i) / sqrt(i);
    sigma_hat_vect(i) = sigma_hat;
end

upper_bound2 = theta_hat_vect + 3*sigma_hat_vect;
lower_bound2 = theta_hat_vect - 3*sigma_hat_vect;
rn_hat = sigma_hat_vect ./ theta_hat_vect;

% phi = ones(n,1);
% phi\y                                % theta_hat verification
% mean((y - theta_stim).^2)/sqrt(n)    % sigma_hat verification

figure(1)
N = (1:1:n);
p1 = plot(N, theta_real, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 4);
hold on;
p2 = plot(N, theta_hat_vect, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 4);
hold on;
p3 = plot(N, upper_bound, 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 4);
hold on;
plot(N, lower_bound, 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 4);
hold on;
p4 = plot(N, upper_bound2, 'Color', [0.4940, 0.1840, 0.5560], 'LineWidth', 4);
hold on;
plot(N, lower_bound2, 'Color', [0.4940, 0.1840, 0.5560], 'LineWidth', 4);
hold off;
grid on;
% title('Recursive least squares (with known and unknown noise variance)', 'FontSize', 24)
% title('Recursive least squares (with known noise variance)', 'FontSize', 24)
% title('Recursive least squares (with unknown noise variance)', 'FontSize', 24)
xlabel('number of measurements', 'FontSize', 24)
ylabel({'$\hat{\theta}$'}, 'Interpreter', 'latex', 'FontSize', 24)
legend([p1 p2 p3 p4], {'$\theta$', '$\hat{\theta}$', '$\hat{\theta} \pm 3\sigma^\theta$', '$\hat{\theta} \pm 3\hat{\sigma}^\theta$'}, 'Interpreter', 'latex', 'FontSize', 36);
% legend([p1 p2 p3], {'$\theta$', '$\hat{\theta}$', '$\hat{\theta} \pm 3\sigma^\theta$'}, 'Interpreter', 'latex', 'FontSize', 36);
% legend([p1 p2 p4], {'$\theta$', '$\hat{\theta}$', '$\hat{\theta} \pm 3\hat{\sigma}^\theta$'}, 'Interpreter', 'latex', 'FontSize', 36);
set(gca,'FontSize',24)
% [229 57 53]/255 red
% [0 150 136]/255 aqua

figure(3)
N = (1:1:n);
p1 = plot(N, theta_real, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 4);
hold on;
p2 = plot(N, y, '.', 'MarkerSize', 20, 'MarkerEdgeColor', [0.8500, 0.3250, 0.0980]);
hold off;
grid on;
xlabel('number of measurements')
ylabel({'$\theta$'}, 'Interpreter', 'latex')
legend([p1 p2], {'$\theta$', '$y_i$'}, 'Interpreter', 'latex', 'FontSize', 36);
set(gca,'FontSize',24)

figure(4)
N = (1:1:n);
p1 = plot(N, rn, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 4);
hold on;
p2 = plot(N, rn_hat, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 4);
hold off;
grid on;
xlabel('number of measurements')
ylabel({'$r_n$'}, 'Interpreter', 'latex')
legend([p1 p2], {'$r_n$', '$\hat{r}_n$'}, 'Interpreter', 'latex', 'FontSize', 36);
set(gca,'FontSize',24)

%% 2.3 Recursive least squares with forgetting factor

n = 1000;
sigma = 1;
ni2 = randn(2*n,sigma);
theta_real = [ones(n, 1); 2*ones(n, 1)];
y = theta_real + ni2;

lambda = 0.98;
Z = 0;
z = 0;

theta_hat_vect = zeros(2*n, 1);
sigma_hat_vect = zeros(2*n, 1);

for i = 1:2*n
    Z = 1 + lambda*Z;
    z = y(i) + lambda*z;
    theta_hat = z/Z;
%     theta_hat = theta_hat + (y(i) - theta_hat)/Z;     % same solution
    theta_hat_vect(i) = theta_hat;
    
    sum = 0;
    for j = 1:i
        sum = sum + (y(j) - theta_hat)^2;
    end
    sigma_hat = sqrt(sum / i) / sqrt(i);
    sigma_hat_vect(i) = sigma_hat;
end

upper_bound = theta_hat_vect + 3*sigma_hat_vect;
lower_bound = theta_hat_vect - 3*sigma_hat_vect;

% figure(3)
% N = (1:1:n);
% p1 = plot(N, theta_real, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 4);
% hold on;
% p3 = plot(N, y, 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 4);
% hold on;
% plot(N, lower_bound, 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 4);
% p2 = plot(N, theta_hat_vect, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 4);
% hold on;
% hold off;
% grid on;
% % title('grafico RLS (known noise variance)')
% xlabel('number of measurements')
% ylabel({'$\hat{\theta}$'}, 'Interpreter', 'latex')
% legend([p1 p2 p3], {'$\theta$', '$\hat{\theta}$', '$\hat{\theta} \pm 3\hat{\sigma}^\theta$'}, 'Interpreter', 'latex', 'FontSize', 36, 'Location','southeast');
% set(gca,'FontSize',24)

% figure(4)
% N = (1:1:2*n);
% p1 = plot(N, theta_real, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 4);
% hold on;
% p3 = plot(N, upper_bound, 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 4);
% hold on;
% plot(N, lower_bound, 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 4);
% p2 = plot(N, theta_hat_vect, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 4);
% hold on;
% hold off;
% grid on;
% xlabel('number of measurements')
% ylabel({'$\hat{\theta}$'}, 'Interpreter', 'latex')
% legend([p1 p2 p3], {'$\theta$', '$\hat{\theta}$', '$\hat{\theta} \pm 3\hat{\sigma}^\theta$'}, 'Interpreter', 'latex', 'FontSize', 36, 'Location','southeast');
% set(gca,'FontSize',24)
