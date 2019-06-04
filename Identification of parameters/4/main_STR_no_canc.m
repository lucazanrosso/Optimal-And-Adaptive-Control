%%  Nominal model definitions

%   mechanical parameters
Jm  = 0.000027;         %   motor rotor inertia
Jd1 = 1.08e-5;          %   small disk inertia
Jd2 = 0.000149;         %   large disk inertia
J   = Jm+Jd1;           %   tot load inertia
b   = 2.4299e-004;      %   viscous damping coeff [not known!!]

%   electrical parameters
Ki  = 2;                % current loop transconductance [A/V]
Kt  = 0.071;            % motor torque const [Nm/A]

%   nominal model (continuous time)
numPc = Ki*Kt;
denPc = [J b 0];
sysPc = tf(numPc, denPc);

%   nominal models (discrete time)  
Ts = 0.001;                         %   sampling time [s]
sysP = c2d(sysPc, Ts, 'zoh');       %   discretization
[numP, denP] = tfdata(sysP, 'v');
numP = numP(2:end);                 %   remove leading zeros coefficient


%%  Additional simulink model definitions    

%   plant I/O prefilters
numHu = 1;
denHu = [1 0];

numHy = [1 -1];
denHy = denHu;

sysHu = tf(numHu, denHu, Ts);
sysHy = tf(numHy, denHy, Ts);


%% Inizialization of PD parameters
% comment on all parameters after running simulink.

Kp = 1;
Kd = 1;
tau_c = 1;

Kp_b = 1;
Kd_b = 1;
tau_c_b = 1;