%% 4.2 Identification-based Tuning for Control Design

%% System without b_hat
s = tf('s');
P = ( Ki * Kt ) / ( J * s^2 );

m_phi  = 35*pi/180; %valore ricavato dal grafico nel pdf
ts = 0.5; %tempo di settlement da progetto
xi = 0.6; %valore da tabella pdf 
w_co = 4 / ( ts * xi );
[mod_sys_co, phase_sys] = bode(P, w_co);
m_phi_sys = pi + (phase_sys * pi/180)

a = 1/mod_sys_co
alpha = m_phi - m_phi_sys

Kp = a * cos( alpha );
Kd = a * sin( alpha )/w_co;
tau_c = 1/( w_co * 100 ); %polo inserito ad hoc per rendere il sistema fisicamente realizzabile

%Parametri del PD
C = Kp + Kd * ( s / ( 1 + tau_c * s ) );
[mod_contr , phase_contr ] = bode(P , w_co );
W = ( C * P ) / ( 1 + C * P );

% [magP0, phaseP0] = bode(P, 0); 
% [magC0, phaseC0] = bode(C, 0);
% error = 1/(1 + magP0*magC0)*100; %ovvio perchè rimane un polo in 0?

figure(1);
bode(P, 'c'); %sistema
hold on;
bode(C, 'r'); %controllore PD
hold on;
bode(P*C, 'k'); %guadagno d'anello L
hold on;
bode(W, 'm'); %funzione in retroazione
grid on;
title('Bode')
legend('P', 'C', 'L', 'W');
hold off;

% figure(2);
% nyquist(P*C);

%% System with b_hat
b_hat = b_hat(20001);
P_b = ( Ki * Kt ) / (s * ( J * s + b_hat) ); %equazione del sistema

% Specifiche di progetto
m_phi  = 60*pi/180; %valore ricavato dal grafico nel pdf
ts = 0.5; %tempo di settlement da progetto
xi = 0.6; %valore da tabella pdf 
w_co = 4/( ts * xi )
% w_co = 15
[  mod_sys_b, phase_sys_b ] = bode(P_b , w_co );
m_phi_sys_b = pi + phase_sys_b*pi/180

a_b = 1/mod_sys_b
alpha_b = m_phi - m_phi_sys_b

% Kp_b = a_b*1.5
Kp_b = a_b * cos( alpha_b ) 
Kd_b = a_b * sin( alpha_b )/w_co
tau_c_b = 1/( w_co * 100); %polo inserito ad hoc per rendere il sistema fisicamente realizzabile
% tau_c_b = 0.01*Kd_b/Kp_b;

C_b = Kp_b + Kd_b * ( s / ( 1 + tau_c_b * s ) );
[  mod_contr_b , phase_contr_b ] = bode(P_b , w_co );
W_b = ( C_b * P_b ) / ( 1 + C_b * P_b );
% 
% figure(2);
% bode(P_b, 'c'); %sistema
% hold on;
% bode(C_b, 'r'); %controllore PD
% hold on;
% bode(P_b*C_b, 'k'); %guadagno d'anello L
% hold on;
% bode(W_b, 'm'); %funzione in retroazione
% grid on;
% title('Bode')
% legend('P', 'C', 'L', 'W');
% hold off;

% figure(2);
% nyquist(P*C);