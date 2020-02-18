%Example of CAZI 
clear;
clc;
%Data of the problem
N  = 100; %number of iterations
max_segments = 4; %max number of segments forming the zonotopes
sigma = 0.03;
initial_thetas = [0.9;0.8];
gamma = [0.05;0.05]; %expansion vector
H = [-0.2,0.1;0.1,0.2]; %initial generators

%generation of the time-varying parameter theta1
theta1 = zeros(1,N); %initialisation of theta_1
for i = 1:N
    thetap = 0.9;
    if (i-1) > 0 
        thetap = theta1(i-1);
    end
    theta1(i) = thetap + 0.05 - i * 0.0008;
end

%generation of the time-varying parameter theta2
theta2 = zeros(1,N); %initialisation of theta_2
for i = 1:N
    thetap = 0.8;
    if (i-1) > 0 
        thetap = theta2(i-1);
    end
    theta2(i) = thetap + 0.05 - i * 0.001;
end

%creating an array with the parameters
parameters = cell(1,2);
parameters{1} = theta1;
parameters{2} = theta2;

%generation of the bounded error
e = zeros(1,N);
for i = 1:N
    e(i) = -1 + 2 * rand(1);
    e(i) = e(i) / 5;
end

%calculation of the output vector
x = zeros(1,N); %instantiation of the output vector
for i = 1:N
    xp1 = 0.1; %initial condition
    xp2 = 0.1; %initial condition
    if (i-2) > 0
        xp2 = x(i-2);
    end
    if (i-1) > 0
        xp1 = x(i-1);
    end
    x(i) = theta1(i) * xp1 - theta2(i) * xp2 +  e(i); 
end

B1 = cell(1,N);
B2 = cell(1,N);
Phi_u = cell(1,N);
Phi_l = cell(1,N);

k = 1;
while k <= N
    xp1 = 0.1;
    xp2 = 0.1;
    if (k-2) > 0
        xp2 = x(k-2);
    end
    if (k-1) > 0
        xp1 = x(k-1);
    end
    
    b1 = x(k) - e(k);
    B1{k} = b1;
    
    b2 = x(k) - e(k);
    B2{k} = b2;

    phi_u = ones(1,2);
    phi_u(1) = xp1 + abs(sigma * xp1);
    phi_u(2) = -xp2 + abs(sigma * xp2);
    Phi_u{k} = phi_u;
    
    phi_l = ones(1,2);
    phi_l(1) = xp1 - abs(sigma * xp1);
    phi_l(2) = -xp2 - abs(sigma * xp2);
    Phi_l{k} = phi_l;
    k = k + 1;
end

%function call
tic %start clock
CAZI(initial_thetas, H, gamma, B1, B2, Phi_u, Phi_l, max_segments, N, parameters)
toc %stop clock