n = 4000;                      % number of samples
u = data2.signals.values(:,2);   % input samples
y_f = data2.signals.values(:,1); % output samples
t = data.time;
time = t(1:n+1);

y_f_b = y_f;
overshoot = 1.1*u(n)*ones(n+1, 1);
upper_settling = 1.02*u(n)*ones(n+1, 1);
lower_settling = 0.98*u(n)*ones(n+1, 1);

figure(1)
p1 = plot(time, u, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 5);
hold on;
p4 = plot(time, overshoot, 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 3);
hold on;
p5 = plot(time, upper_settling, 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 3);
hold on;
plot(time, lower_settling, 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 3);
hold on;
p2 = plot(time, y_f, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 5);
hold on;
p3 = plot(time, y_f_b, 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 5);
hold off;
grid on;
xlabel({'time [s]'}, 'Interpreter', 'latex')
ylabel({'position [rad]'}, 'Interpreter', 'latex')
legend([p1 p2 p3 p4 p5], {'reference signal', 'response signal ($b = 0$)', 'response signal ($b = \hat{b}$)','overshoot limit ($10\%$)', 'settiling time limit ($0.5s$)'}, 'Interpreter', 'latex', 'FontSize', 36, 'Location', 'SouthEast');
set(gca,'FontSize',24)
