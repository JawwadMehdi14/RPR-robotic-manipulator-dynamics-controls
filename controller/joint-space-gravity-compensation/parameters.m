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

%Joint variables values
qd = [pi/3, 0.1, 0.2];

%PD control values
Kp = [[90,0,0];
      [0,90,0];
      [0,0,90]];

Kd = [[60,0,0];
      [0,60,0];
      [0,0,60]];