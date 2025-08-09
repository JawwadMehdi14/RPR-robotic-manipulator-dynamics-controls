function [sys, x, str, ts] = Manipulator(t,x, u, flag, parameters)
    if flag == 0 %initialization
        x = zeros(6,1); % Initializes the state vector
        str = []; % empty state string
        ts = [0 0]; % sample times    

        sizes = simsizes;
        sizes.NumContStates = 6; %number of continuous states^
        sizes.NumDiscStates = 0; %number of discrete states
        sizes.NumOutputs = 6; % number of outputs^
        sizes.NumInputs = 3; %number of inputs
        sizes.DirFeedthrough = 0; % Does u depends on the output
        sizes.NumSampleTimes  = 1; %Number of function to run in time t

        sys = simsizes(sizes);

    elseif flag == 1   % ... (Calculate q_dot_dot using the inverse of inertia matrix, coriolis matrix, and gravity matrix)

        q_dot = x(1:3);
        q = x(4:6);
        tau = u(1:3);
        parameters.q = q;
        parameters.q_dot = q_dot;
        B = double(inertiaMatrix(parameters));
        C = double(corolisMatrix(parameters));
        G = double(gravityMatrix(parameters));
        
        q_dot_dot = pinv(B)*(tau - C*q_dot - G);

        sys = [q_dot_dot; q_dot];

    elseif flag == 3 %generate the output
        sys = x;
    end
end

function B = inertiaMatrix(parameters) %This function calculates the inertia matrix (B) based on the joint angles (q).
    q1 = parameters.q(1);
    q2 = parameters.q(2);
    q3 = parameters.q(3);

    m1 = parameters.m(1);
    m2 = parameters.m(2);
    m3 = parameters.m(3);

    a1 = 0.02;
    a2 = 0.03;
    a3 = 0.01;
    b1 = 0.1256;
    b2 = 0.03;
    b3 = 0.01;
    c1 = 0.4;
    c2 = 0.46;
    c3 = 0.15;

B = [[(a2^2*m2)/12 + (3*a1^4*m1)/2 + (a3^2*m3)/12 + (3*b1^4*m1)/2 + (b3^2*m3)/12 + c1^2*m1 + c1^2*m2 + c1^2*m3 + (7*c2^2*m2)/12 + (c3^2*m3)/2 + m3*q3^2 + 3*a1^2*b1^2*m1 - c1*c2*m2 + c1*c3*m3 + 2*c1*m3*q3 + c3*m3*q3,       0,  0]
[                                                                                                                                                                                                               0, m2 + m3,  0]
[                                                                                                                                                                                                               0,       0, m3]];
 
end

function C = corolisMatrix(parameters) %This function calculates the Coriolis matrix (C) based on joint angles and their velocities (q and dq).
    q1 = parameters.q(1);
    q2 = parameters.q(2);
    q3 = parameters.q(3);

    q1_dot = parameters.q_dot(1);
    q2_dot = parameters.q_dot(2);
    q3_dot = parameters.q_dot(3);

    m1 = parameters.m(1);
    m2 = parameters.m(2);
    m3 = parameters.m(3);

    a1 = 0.02;
    a2 = 0.03;
    a3 = 0.01;
    b1 = 0.1256;
    b2 = 0.03;
    b3 = 0.01;
    c1 = 0.4;
    c2 = 0.46;
    c3 = 0.15;
                                                                                     
C = [[ (m3*q3_dot*(2*c1 + c3 + 2*q3))/2, 0, (m3*q1_dot*(2*c1 + c3 + 2*q3))/2]
[                                0, 0,                                0]
[-(m3*q1_dot*(2*c1 + c3 + 2*q3))/2, 0,                                0]];
 
end

function G = gravityMatrix(parameters)
    q1 = parameters.q(1);
    q2 = parameters.q(2);
    q3 = parameters.q(3);
    
    g = parameters.g;

    m1 = parameters.m(1);
    m2 = parameters.m(2);
    m3 = parameters.m(3);

    a1 = 0.02;
    a2 = 0.03;
    a3 = 0.01;
    b1 = 0.1256;
    b2 = 0.03;
    b3 = 0.01;
    c1 = 0.4;
    c2 = 0.46;
    c3 = 0.15;

    G = [(c1*g*m1*cos(q1))/2 + g*m2*cos(q1)*(c1 - c2/2) + g*m3*cos(q1)*(c1 + c3/2 + q3)
        0;
        g*m3*sin(q1)];
end