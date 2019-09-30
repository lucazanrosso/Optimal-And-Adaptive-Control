%% 3.6  Evaluate the performance of the systems for r = {1,10,100,1000}
dati = sprintf('%d.mat',20);
load(dati)


%% 3.8 Compare the performance between the LQR and the Pole Placement approach by using the best performing values for r and omega_n 

% theta1_pp = data.signals.values(:,1);
% theta2_pp = data.signals.values(:,2);
% theta1_dot_pp = data.signals.values(:,3);
% theta2_dot_pp  = data.signals.values(:,4);
% t = data.time;
% 
% % theta1_lqr = theta1_pp;
% % theta2_lqr = theta2_pp;
% % theta1_dot_lqr = theta1_dot_pp;
% % theta2_dot_lqr = theta2_dot_pp;
% 
% figure(1)
% p2 = plot(t, theta1_pp, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 3);
% hold on;
% p1 = plot(t, theta1_lqr, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 3);
% hold off;
% grid on;
% xlabel({'time [s]'}, 'Interpreter', 'latex')
% ylabel({'$\theta_1$ [rad]'}, 'Interpreter', 'latex')
% title({'Comparison between lqr and pp for $\theta_1$'}, 'Interpreter', 'latex')
% legend([p1 p2], {'$\theta_1$ (lqr)', '$\theta_1$ (pp)'}, 'Interpreter', 'latex', 'FontSize', 36, 'Location', 'SouthEast');
% set(gca,'FontSize',24)
% 
% figure(2)
% p4 = plot(t, theta2_pp, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 3);
% hold on;
% p3 = plot(t, theta2_lqr, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 3);
% hold off;
% grid on;
% xlabel({'time [s]'}, 'Interpreter', 'latex')
% ylabel({'$\theta_2$ [rad]'}, 'Interpreter', 'latex')
% title({'Comparison between lqr and pp for $\theta_2$'}, 'Interpreter', 'latex')
% legend([p3 p4], {'$\theta_2$ (lqr)', '$\theta_2$ (pp)'}, 'Interpreter', 'latex', 'FontSize', 36, 'Location', 'SouthEast');
% set(gca,'FontSize',24)
% 
% figure(3)
% p2 = plot(t, theta1_dot_pp, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 3);
% hold on;
% p1 = plot(t, theta1_dot_lqr, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 3);
% hold off;
% grid on;
% xlabel({'time [s]'}, 'Interpreter', 'latex')
% ylabel({'$\dot{\theta_1}$ [rad]'}, 'Interpreter', 'latex')
% title({'Comparison between lqr and pp for $\dot{\theta_1}$'}, 'Interpreter', 'latex')
% legend([p1 p2], {'$\dot{\theta_1}$ (lqr)', '$\dot{\theta_1}$ (pp)'}, 'Interpreter', 'latex', 'FontSize', 36, 'Location', 'SouthEast');
% set(gca,'FontSize',24)
% 
% figure(4)
% p2 = plot(t, theta2_dot_pp, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 3);
% hold on;
% p1 = plot(t, theta2_dot_lqr, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 3);
% hold off;
% grid on;
% xlabel({'time [s]'}, 'Interpreter', 'latex')
% ylabel({'$\dot{\theta_2}$ [rad]'}, 'Interpreter', 'latex')
% title({'Comparison between lqr and pp for $\dot{\theta_2}$'}, 'Interpreter', 'latex')
% legend([p1 p2], {'$\dot{\theta_2}$ (lqr)', '$\dot{\theta_2}$ (pp)'}, 'Interpreter', 'latex', 'FontSize', 36, 'Location', 'SouthEast');
% set(gca,'FontSize',24)



%% 4.5 Compare the Pole Placement design with the LQR design for the Dynamic Observer
% theta1_pp = data.signals.values(:,1);
% theta2_pp = data.signals.values(:,2);
% theta1_dot_pp = data.signals.values(:,3);
% theta2_dot_pp  = data.signals.values(:,4);
% t_pp = data.time;
% 
% theta1_pp2 = data2.signals.values(:,1);
% theta2_pp2 = data2.signals.values(:,2);
% theta1_dot_pp2 = data2.signals.values(:,3);
% theta2_dot_pp2  = data2.signals.values(:,4);
% 
% % theta1_lqr = theta1_pp;
% % theta2_lqr = theta2_pp;
% % theta1_dot_lqr = theta1_dot_pp;
% % theta2_dot_lqr = theta2_dot_pp;
% % theta1_lqr2 = theta1_pp2;
% % theta2_lqr2 = theta2_pp2;
% % theta1_dot_lqr2 = theta1_dot_pp2;
% % theta2_dot_lqr2 = theta2_dot_pp2;
% % t_lqr = t_pp;
% 
% figure(1)
% p2 = plot(t_pp, theta1_pp, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 3);
% hold on;
% p1 = plot(t_lqr, theta1_lqr, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 3);
% hold off;
% grid on;
% xlabel({'time [s]'}, 'Interpreter', 'latex')
% ylabel({'$\theta_1$ [rad]'}, 'Interpreter', 'latex')
% title({'Comparison between lqr and pp for $\theta_1$ (with Observer)'}, 'Interpreter', 'latex')
% legend([p1 p2], {'$\theta_1$ (lqr)', '$\theta_1$ (pp)'}, 'Interpreter', 'latex', 'FontSize', 36, 'Location', 'SouthEast');
% set(gca,'FontSize',24)
% 
% figure(2)
% p4 = plot(t_pp, theta2_pp, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 3);
% hold on;
% p3 = plot(t_lqr, theta2_lqr, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 3);
% hold off;
% grid on;
% xlabel({'time [s]'}, 'Interpreter', 'latex')
% ylabel({'$\theta_2$ [rad]'}, 'Interpreter', 'latex')
% title({'Comparison between lqr and pp for $\theta_2$ (with Observer)'}, 'Interpreter', 'latex')
% legend([p3 p4], {'$\theta_2$ (lqr)', '$\theta_2$ (pp)'}, 'Interpreter', 'latex', 'FontSize', 36, 'Location', 'SouthEast');
% set(gca,'FontSize',24)
% 
% figure(3)
% p2 = plot(t_pp, theta1_dot_pp, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 3);
% hold on;
% p1 = plot(t_lqr, theta1_dot_lqr, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 3);
% hold off;
% grid on;
% xlabel({'time [s]'}, 'Interpreter', 'latex')
% ylabel({'$\dot{\theta_1}$ [rad]'}, 'Interpreter', 'latex')
% title({'Comparison between lqr and pp for $\dot{\theta_1}$ (with Observer)'}, 'Interpreter', 'latex')
% legend([p1 p2], {'$\dot{\theta_1}$ (lqr)', '$\dot{\theta_1}$ (pp)'}, 'Interpreter', 'latex', 'FontSize', 36, 'Location', 'SouthEast');
% set(gca,'FontSize',24)
% 
% figure(4)
% p2 = plot(t_pp, theta2_dot_pp, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 3);
% hold on;
% p1 = plot(t_lqr, theta2_dot_lqr, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 3);
% hold off;
% grid on;
% xlabel({'time [s]'}, 'Interpreter', 'latex')
% ylabel({'$\dot{\theta_2}$ [rad]'}, 'Interpreter', 'latex')
% title({'Comparison between lqr and pp for $\dot{\theta_2}$ (with Observer)'}, 'Interpreter', 'latex')
% legend([p1 p2], {'$\dot{\theta_2}$ (lqr)', '$\dot{\theta_2}$ (pp)'}, 'Interpreter', 'latex', 'FontSize', 36, 'Location', 'SouthEast');
% set(gca,'FontSize',24)