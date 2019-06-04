%% 4 Identification of discrete time linear systems

%% 4.1 Parameter identification
n = 20000; % samples
u = data.signals.values(:,1);  %salvo i valori degli ingressi
y_f = data.signals.values(:,2); %salvo i valori dell'uscita

I = eye(3);
PHI_0 =  [ -y_f(1) u(2) u(1); -y_f(2) u(3) u(2); -y_f(3) u(4) u(3);];
Y_0 = [y_f(2); y_f(3); y_f(4)];
P = inv(PHI_0'*PHI_0);

theta_hat = P*PHI_0'*Y_0;
a0_hat = zeros(n+1, 1);
sigma_a0_hat = zeros(n+1, 1);

for i = 5 : n+1    
    phi =  [ -y_f(i-1); u(i); u(i-1) ];  
    K = P * phi / (1 + phi'*P*phi);
    P = (I - K*phi')*P;
    theta_hat = theta_hat + K*(y_f(i) - phi'*theta_hat);     
    a0_hat(i) = theta_hat(1);
    
    sum = 0;
    for j = 5:i
        sum = sum + (y_f(j) - phi'*theta_hat)^2;
    end
    sigma_a0_hat(i) = sqrt(P(1,1) * sum / (i - 4)); % parto da 5, quindi devo dividere per i - 4
end

b_real = ones(n+1, 1) * b;
b_hat = -(J / Ts) * log(-a0_hat)
b_upper_bound = (J / Ts) * log(-a0_hat + 3 * sigma_a0_hat);
b_lower_bound = (J / Ts) * log(-a0_hat - 3 * sigma_a0_hat);

figure(1)
N = (1:n+1);
plot(N, b_real, 'b'); % blu
hold on;
plot(N, b_hat, 'k'); % blu
hold on;
plot(N, b_upper_bound, 'g'); % rosso
hold on;
plot(N, b_lower_bound, 'g'); % rosa
hold off;
grid on;
title('grafico RLS (known noise variance)')
xlabel('numero misurazioni')
ylabel('teta stim e \sigma_n')
legend('\theta_{real}','\theta_{hat}','upper bound','lower bound')


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
