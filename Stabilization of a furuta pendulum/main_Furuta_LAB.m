% close all; clear all;

%% Define parameters

% mechanical parameters

J1 = 0.000027+0.00257; 
J2 = 4.15e-5;
l1 = 0.16; 
l2 = 0.15/2;
m2 = 0.0206; 
g = 9.81;

% electromechanical parameters for motor and transducers

rad2deg = 180/pi;       %tranformation from radiant to degrees
Kdrv = 2;               % driver gain (voltage-controlled current generator) [A/V]
Kt = 0.071;             % motor current-to-torque factor [Nm/A]
pulse2rad = 2*pi/2000;  % encoder transduction factor: 2000 pulses per revolution
Vmax = 3;               % maximum allowed voltage for the driver
Tderiv = 1/200;         % high frequecy pole for derivator 


% Define the linearized stateâˆ’space model
gam = J1*J2 + J1*m2*l2^2 + J2*m2*l1^2;
a32 = m2^2*l1*l2^2*g/gam; 
a42 = m2*l2*g*(J1 + m2*l1^2)/gam;
b31 = (J2 + m2*l2^2)/gam; 
b32 = m2*l1*l2/gam;
b41 = m2*l1*l2 /gam;
b42 = (J1 + m2*l1^2)/gam;

A=[0 0 1 0;0 0 0 1;0 a32 0 0; 0 a42 0 0];
B=[0 0 ; 0 0; b31 b32; b41 b42 ];
x0 = [0; 3*pi/180; 0; 0];
q0 = x0(1:2);
dq0 = x0(3:4);

%%% Disturbance

t_impulse = 10; % starting time of the impulse
delta_impulse = 0.1; %duration of impulse
tau2_impulse = 0.001; %magnitude of the impulse [Nm]

% Design the infinte horizon LQR
% Q = [1 1 1 1];
% r = 100;

%% 2.9 Bryson’s rule 
% amax = 30*pi/180;
% bmax = 10*pi/180;
% adotmax = amax*2*pi*10;
% bdotmax = bmax*2*pi*10;
% umax = 0.5*0.071*2;
% Q = diag([1/amax^2,1/bmax^2,1/adotmax^2,1/bdotmax^2]); 
% r = 1/umax^2;



%% 3
C1 = [1 0 0 0];
Q = C1'*C1;
r = 10;
[Kinf,Pinf,lam] = lqr(A,B(:,1),Q,r);

%% 3.1  Compute explicitly the transfer function in terms of the systems parameters
I = eye(4);
s = tf('s');
P = C1*inv(s*I-A)*B(:,1);
[z, p] = pzmap(P);
P = minreal(P);

%% 3.2 Compute open loop poles and zeros of such transfer function
P2 = C1*inv(-s*I-A)*B(:,1);
P2 = minreal(P2);

%% 3.3 Compute the Symmetric Root Locus for such transfer function
% rlocus(P);

%% 3.4 3.5 Compute the location of the closed loop poles of the LQR solution for r
G = P * P2;
% rlocus(G);

%% 3.7 Compute the alternative solution based on pole placement
wn = 10;
pp = wn*[ -sqrt(2)/2+1i*sqrt(2)/2  -sqrt(2)/2-1i*sqrt(2)/2 -sqrt(3)/2+1i/2  -sqrt(3)/2-1i/2];
Kinf = place(A,B(:,1),pp);

%% 3.8 Compare the performance between the LQR and the Pole Placement approach by using the best performing values for r and omega_n 



%% 4
O = obsv(A,C1);
det(O);

%% 4.1 Design a Dynamic Observer based on Pole placement
L = place(A',C1',pp)';

%% 4.2 Design a Dynamic Observer by using LQR approach
F = A';
G = C1';
H = B(:,1)';
Q = H' * H;

%% 4.3 Compute explicitly the transfer function in terms of the systems parameters
P = H*inv(s*I-F)*G;
[z, p] = pzmap(P);

%% 4.4  Compute the LQR gain for this system for different values of r
[Tau,Pinf,lam] = lqr(F,G,Q,r);
% L = Tau';

%% 4.5 Compare the Pole Placement design with the LQR design for the Dynamic Observer


%prestaz migliori r = 10 e wn = 10