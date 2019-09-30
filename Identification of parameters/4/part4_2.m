%% 4.2 Identification-based Tuning for Control Design

%% System without b_hat
s = tf('s');
P = ( Ki * Kt ) / ( J * s^2 );

m_phi  = 60*pi/180; %valore ricavato dal grafico nel pdf
ts = 0.5; %tempo di settlement da progetto
xi = 0.6; %valore da tabella pdf 
w_co = 4 / ( ts * xi );
[mod_sys_co, phase_sys] = bode(P, w_co);
m_phi_sys = pi + (phase_sys * pi/180)

a = 1/mod_sys_co
alpha = m_phi - m_phi_sys

Kp = a * cos( alpha );
Kd = a * sin( alpha )/w_co;
tau_c = 1/( w_co * 10 ); %polo inserito ad hoc per rendere il sistema fisicamente realizzabile

C = Kp + Kd * ( s / ( 1 + tau_c * s ) );            % PD
W = ( C * P ) / ( 1 + C * P );

% [magP0, phaseP0] = bode(P_b, 0); 
% [magC0, phaseC0] = bode(C_b, 0);
% error = 1/(1 + magP0*magC0)*100 % the error is always zero because there is a pole in the origin

% figure(1);
% bode(P, 'c'); %sistema
% hold on;
% bode(C, 'r'); %controllore PD
% hold on;
% bode(P*C, 'k'); %guadagno d'anello L
% hold on;
% bode(W, 'm'); %funzione in retroazione
% grid on;
% title('Bode')
% legend('P', 'C', 'L', 'W');
% hold off;

w = logspace(-1,5,100)
[mag1,phase1,wout1] = bode(P, w); %sistema
[mag2,phase2,wout2] = bode(C, w); %controllore PD
[mag3,phase3,wout3] = bode(P*C, w); %guadagno d'anello L
[mag4,phase4,wout4] = bode(W, w); %funzione in retroazione
figure(3);
p1 = semilogx(wout1,  20*log10(squeeze(mag1)), 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 5);
hold on;
p2 = semilogx(wout2, 20*log10(squeeze(mag2)), 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 5);
hold on;
p3 = semilogx(wout3, 20*log10(squeeze(mag3)), 'Color', [0.3010 0.7450 0.9330], 'LineWidth', 5);
hold on;
p4 = semilogx(wout4, 20*log10(squeeze(mag4)), 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 5);
hold off;
grid on;
xlabel('Frequency [rad/s]')
ylabel('Magnitude [dB]')
legend([p1 p2 p3 p4], {'$P(s)$', '$C(s)$', '$G(s)$', '$W(s)$'}, 'Interpreter', 'latex', 'FontSize', 36);
set(gca,'FontSize',24)
figure(4);
p1 = semilogx(wout1, squeeze(phase1), 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 5);
hold on;
p2 = semilogx(wout2, squeeze(phase2), 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 5);
hold on;
p3 = semilogx(wout3, squeeze(phase3), 'Color', [0.3010 0.7450 0.9330], 'LineWidth', 5);
hold on;
p4 = semilogx(wout4, squeeze(phase4), 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 5);
hold off;
grid on;
xlabel('Frequency [rad/s]')
ylabel('Phase [deg]')
legend([p1 p2 p3 p4], {'$P(s)$', '$C(s)$', '$G(s)$', '$W(s)$'}, 'Interpreter', 'latex', 'FontSize', 36);
set(gca,'FontSize',24)


%% System with b_hat
P_b = ( Ki * Kt ) / (s * ( J * s + b_hat) ); %equazione del sistema

m_phi  = 60*pi/180; %valore ricavato dal grafico nel pdf
ts = 0.5; %tempo di settlement da progetto
xi = 0.6; %valore da tabella pdf 
w_co = 4/( ts * xi )
[  mod_sys_b, phase_sys_b ] = bode(P_b , w_co );
m_phi_sys_b = pi + phase_sys_b*pi/180

a_b = 1/mod_sys_b
alpha_b = m_phi - m_phi_sys_b

Kp_b = a_b * cos( alpha_b ) 
Kd_b = a_b * sin( alpha_b )/w_co
tau_c_b = 1/( w_co * 10); %polo inserito ad hoc per rendere il sistema fisicamente realizzabile
% tau_c_b = 0.01*Kd_b/Kp_b;

C_b = Kp_b + Kd_b * ( s / ( 1 + tau_c_b * s ) );
[  mod_contr_b , phase_contr_b ] = bode(P_b , w_co );
W_b = ( C_b * P_b ) / ( 1 + C_b * P_b );

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

[mag1,phase1,wout1] = bode(P_b, w); %sistema
[mag2,phase2,wout2] = bode(C_b, w); %controllore PD
[mag3,phase3,wout3] = bode(P_b*C_b, w); %guadagno d'anello L
[mag4,phase4,wout4] = bode(W_b, w); %funzione in retroazione
figure(5);
p1 = semilogx(wout1,  20*log10(squeeze(mag1)), 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 5);
hold on;
p2 = semilogx(wout2, 20*log10(squeeze(mag2)), 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 5);
hold on;
p3 = semilogx(wout3, 20*log10(squeeze(mag3)), 'Color', [0.3010 0.7450 0.9330], 'LineWidth', 5);
hold on;
p4 = semilogx(wout4, 20*log10(squeeze(mag4)), 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 5);
hold off;
grid on;
xlabel('Frequency [rad/s]')
ylabel('Magnitude [dB]')
legend([p1 p2 p3 p4], {'$P(s)$', '$C(s)$', '$G(s)$', '$W(s)$'}, 'Interpreter', 'latex', 'FontSize', 36);
set(gca,'FontSize',24)
figure(6);
p1 = semilogx(wout1, squeeze(phase1), 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 5);
hold on;
p2 = semilogx(wout2, squeeze(phase2), 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 5);
hold on;
p3 = semilogx(wout3, squeeze(phase3), 'Color', [0.3010 0.7450 0.9330], 'LineWidth', 5);
hold on;
p4 = semilogx(wout4, squeeze(phase4), 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 5);
hold off;
grid on;
xlabel('Frequency [rad/s]')
ylabel('Phase [deg]')
legend([p1 p2 p3 p4], {'$P(s)$', '$C(s)$', '$G(s)$', '$W(s)$'}, 'Interpreter', 'latex', 'FontSize', 36);
set(gca,'FontSize',24)
