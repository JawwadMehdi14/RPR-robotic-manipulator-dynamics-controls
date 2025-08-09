clear;
clc;

%sine parameters
A = 1;
Fc = 1;

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


q1 = 0.5629;
q2 = -0.6600;
q3 = 0.7250;

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


    H = [ [0, -sin(q1), cos(q1),              cos(q1)*(c1 + c3 + q3)]
[1,        0,       0,                           - c2 - q2]
[0,  cos(q1), sin(q1), r1 + sin(q1)*(c3 + q3) + c1*sin(q1)]
[0,        0,       0,                                   1]
];


eulers = rotm2eul(H(1:3,1:3),"ZYZ");


%operational space desired values
xd = [H(1:3,4); eulers'];

%PD control values
Kp = diag([200,200,200,50,50,50]);

Kd = diag([150,150,150,20,20,20]);

