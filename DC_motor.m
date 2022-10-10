clc;
clear all;
close all; 
warning off 

% parameters 
Ra = 1.203;                         %Resistance
L = 0.005584;                       %Henry _ Inductance
Ke = 0.08574;                       %Motor velocity constant
Kt = 1.0005*Ke;                     %Torque constant
J = 0.00014166;                     %Inertia
fr = 0.000245;                      %friction coefficent
Ts = 1e-3;                          %Sampling time
N = 10;                             % Prediction horizon
I0 = 0.2;                           %initial current
n0 = 70.3;                          %initial speed
x0 = [I0 n0]';
xref = [0 0]';

% uncertainty input 
Bw = [-0.0085 -0.0006
    -0.0603  0.0002];
W = [1
    1];
%State Space Model of System
As = [-Ra/L -Ke/L;Kt/J -fr/J];
Bs = [1/L;0];
Cs = [1 0;0 1];
Ds = [0 0]';
sys = ss(As,Bs,Cs,Ds);
Gs = tf(sys);

%Discretized system 
dis_sys = c2d(sys,Ts);
Ad = dis_sys.A;
Bd = dis_sys.B;
Cd = dis_sys.C;
Dd = dis_sys.D;
B_noise = [Bd Bw];
D_noise = [Dd zeros(size(Dd,1),2)];
Gz = tf(dis_sys)
%%
%Constraints of variables
Imin = -1.2;                % minimum current
Imax = 2;                   % maximum current
nmin=-150;                  % minimum speed
nmax=150;                   % maximum speed
xmin=[Imin,nmin]';
xmax=[Imax,nmax]';

Vmin=0;                     % minimum input voltage
sat_min = Vmin;
Vmax=12;                    % maximum input voltage
sat_max = Vmax;
Ax=[1 0;0 1;-1 0;0 -1];
bx=[Imax;nmax;Imax;nmax];
Au=[1;-1];
bu=[12;0];

%% -------------------------------
Q=1000*eye(2);
R=1;
[Kdlqr,P,cp]=dlqr(Ad,Bd,Q,R);
Pr=P;
feedbackSystem = Ad-Bd*Kdlqr;
poles_feedback = eig(feedbackSystem);

QMPC=Q; %optimal Q
R=R;
S=P;
P1 = zeros(size(S));

B_noise1=[Bd Bw];
D_noise1=[Dd Dd Dd];
%%
B_noise=[Bd ones(size(Bd,1),1) Bw];
D_noise=[Dd zeros(size(Dd,1),1) Dd Dd ];
%set power of v_x=0.01 in simulink block
power_vx=0;
%set power of v_y=0.01 in simulink block
power_vy=0.1;
%%

Cnew=[0 1];
Qk=eye(2);          % variance of v_x
Rk=Qk(1)*1000000;               % variance of v_y
NK=0;               % covariance v_x,v_y



open('DC_motor_EX1')
sim('DC_motor_EX1')

