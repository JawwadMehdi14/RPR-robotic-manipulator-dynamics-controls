%sine parameters
A = 1;
Fc = 1;
     c1 = 0.4;
c2 = 0.3;
c3 = 0.16;
r1 = 0.15;

%lengths
lenghts_values = [[0.02, 0.1256, 0.4]; 
               [0.03, 0.03, 0.46]; 
               [0.01, 0.01, 0.15]];
density = 2700; %Almunium density
mass_values = [(pi * lenghts_values(1,1)^2*lenghts_values(1,3)*density), lenghts_values(2,1)*lenghts_values(2,2)*lenghts_values(2,3)*density, lenghts_values(3,1)*lenghts_values(3,2)*lenghts_values(3,3)*density];
         

%masses
params.m = mass_values;

%gravity
params.g = 9.806;

%Friction
params.Fv = [0,0,0];

q1 = pi/3;
q2 = -0.1;
q3 = 0.1;


    r1 = 0.4;
    a1 = 0.02;
    a2 = 0.03;
    a3 = 0.01;
    b1 = 0.1256;
    b2 = 0.03;
    b3 = 0.01;
    c1 = 0.4;
    c2 = 0.46;
    c3 = 0.15;
    
    H = [[0, -sin(q1), cos(q1),              cos(q1)*(c1 + c3 + q3)]
[1,        0,       0,                           - c2 - q2]
[0,  cos(q1), sin(q1), r1 + sin(q1)*(c3 + q3) + c1*sin(q1)]
[0,        0,       0,                                   1]];


%operational space desired values
xd = H(1:3,4)';

%PD control values
Kp = diag([200,200,200]);
Kd = diag([50,50,50]);
Md = diag([0.12,0.12,0.12]);

%force control parameters
fd = [0,0,3];

Kf = diag([50,50,50]);
Ki = diag([50,50,50]);

%enviroment stiffness
parms.Ke = diag([0,0,4]);

%plane position
parms.plane_pos = [0,0,0.3]';
%parms.plane_pos = [0,0,9.0,0,0,0]';

%plane axis
parms.plane_axis = 3; %z-axis
