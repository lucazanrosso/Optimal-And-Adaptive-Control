%%Esercizio 1.1 known noise variance
%versione tre e che sia quella buona da pagina 41 a 53 e da 56 a 59
n = 1000;
sigma = 1; %varianza del rumore, siccome vale 1 ho 68.27% di prob che il rumore stia tra +-1  
ni = randn(n,sigma);
theta_real = ones(n, 1);
y = theta_real + ni; % questa è la y vera, vettore con tutte le misure

%riga che fa i Root least square
% theta_stim = phi\y;

theta_hat = 0;
theta_hat_vect = zeros(n, 1);
sigma_hat = 0;
sigma_hat_vect = zeros(n, 1);

for i= 1:n
    theta_hat = ((i-1)/(i))*theta_hat + (1/(i))*y(i);
    theta_hat_vect(i) = theta_hat;
    sigma_hat = sigma/sqrt(i);
    sigma_hat_vect(i) = sigma_hat;
end

upper_bound = theta_hat_vect + 3*sigma_hat_vect;
lower_bound = theta_hat_vect - 3*sigma_hat_vect;


%% 2.2

theta_hat = 0;
theta_hat_vect = zeros(n, 1);
sigma_hat = 0;
sigma_hat_vect = zeros(n, 1);

for i = 1:n
    theta_hat = ((i-1)/(i))*theta_hat + (1/(i))*y(i);
    theta_hat_vect(i) = theta_hat;
    sum = 0;
    for j = 1:i
        sum = sum + (y(j) - theta_hat)^2;
    end
    sigma_hat = sum / i / sqrt(i);
    sigma_hat_vect(i) = sigma_hat;
end

upper_bound2 = theta_hat_vect + 3*sigma_hat_vect;
lower_bound2 = theta_hat_vect - 3*sigma_hat_vect;

theta_hat;
phi= ones(n,1);
theta_stim = phi\y;
sigma_hat;
mean((y - theta_stim).^2)/sqrt(n);

% figure(1)
% N = (1:1:n);
% plot(N, theta_real, 'b'); % blu
% hold on;
% plot(N, theta_hat_vect, 'k'); % blu
% hold on;
% plot(N, upper_bound, 'r'); % rosso
% hold on;
% plot(N, lower_bound, 'r'); % rosa
% hold on;
% plot(N, upper_bound2, 'g'); % rosso
% hold on;
% plot(N, lower_bound2, 'g'); % rosa
% hold off;
% grid on;
% title('grafico RLS (known noise variance)')
% xlabel('numero misurazioni')
% ylabel('teta stim e \sigma_n')
% legend('\theta_{real}','\theta_{hat}','upper bound','lower bound')
 

%% 2.3

n = 1000;
sigma = 1;
ni = randn(2*n,sigma);
theta_real = [ones(n, 1); 2*ones(n, 1)];
y = theta_real + ni

lambda = 0.99;
Z = 0;
z = 0;

theta_hat = 0;
theta_hat_vect = zeros(2*n, 1);
sigma_hat = 0;
sigma_hat_vect = zeros(2*n, 1);

for i = 1:2*n
    Z = 1 + lambda*Z;
    z = y(i) + lambda*z;
    theta_hat = z/Z;
%     theta_hat = theta_hat + (y(i) - theta_hat)/Z;
    theta_hat_vect(i) = theta_hat;
    
    sum = 0;
    for j = 1:i
        sum = sum + (y(j) - theta_hat)^2;
    end
    sigma_hat = sum / i / sqrt(i);
    sigma_hat_vect(i) = sigma_hat;
end

upper_bound = theta_hat_vect + 3*sigma_hat_vect;
lower_bound = theta_hat_vect - 3*sigma_hat_vect;

figure(1)
N = (1:1:2*n);
plot(N, theta_real, 'b'); % blu
hold on;
plot(N, theta_hat_vect, 'k'); % blu
% hold on;
% plot(N, upper_bound, 'g'); % rosso
% hold on;
% plot(N, lower_bound, 'g'); % rosa
hold off;
grid on;
title('grafico RLS (known noise variance)')
xlabel('numero misurazioni')
ylabel('teta stim e \sigma_n')
legend('\theta_{real}','\theta_{hat}','upper bound','lower bound')
   
   
   
   