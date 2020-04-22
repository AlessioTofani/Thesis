%CAZI and PAZI algorithm with the example from BRAVO

clear;
clc;
%Data of the problem
N  = 100; %number of iterations
max_segments = 40; %max number of segments forming the zonotopes
sigma = 0.05;
initial_thetas = [0.9;0.8;1;0.8;0]; %initial parameters vector
H = 1*eye(5);
gamma = [0.01,0.01,0.01,0.01,0.05]; %expansion factors
time_var = 0;

%input
f = zeros(1,N); %instantiation of the input vector 
f = -5 + (5+5)*rand(1,N); %random bounded input |u| < 5
f(1) = 0;

%disturbance
omega = zeros(1,N); %instantiation of the disturbance vector 
for i = 1: N
    variation = -0.05 + (0.05+0.05) * rand(1);
    omegaprec = 0;
    if (i-1) > 0
        omegaprec = omega(i-1);
    end
    omega(i) = omegaprec + variation;
end
omega(1) = 0;

%generation of the bounded error
e = zeros(1,N);
for i = 1: N
    e(i) = -0.05 + (0.05+0.05) * rand(1); 
end

%calculation of the output vector
x = zeros(1,N); %instantiation of the output vector
for i = 1:N
    %xp1 = 0; %initial condition
    %xp2 = 0; %initial condition
    fp1 = 0; %initial condition
    xp1 = 0.1; %initial condition
    xp2 = 0.1; %initial condition
    if (i-2) > 0
        xp2 = x(i-2);
    end
    if (i-1) > 0
        xp1 = x(i-1);
        fp1 = f(i-1);
    end
    x(i) = 0.9*xp1 -0.8*xp2 - f(i) -0.8*fp1 + omega(i) + e(i);
end

B1 = cell(1,N);
B2 = cell(1,N);
Phi_u = cell(1,N);
Phi_l = cell(1,N);

k = 1;
while k <= N
    xp1 = x(2);
    xp2 = x(1);
    fp1 = 0;
    if (k-2) > 0
        xp2 = x(k-2);
    end
    if (k-1) > 0
        xp1 = x(k-1);
        fp1 = f(k-1);
    end
    
    b1 = x(k) - e(k);
    B1{k} = b1;
    
    b2 = x(k) - e(k);
    B2{k} = b2;

    phi_u = ones(1,5);
    phi_u(1) = xp1 + abs(sigma * xp1);
    phi_u(2) = - xp2 + abs(sigma * xp2);
    phi_u(3) = - f(k) + abs(sigma * f(k));
    phi_u(4) = - fp1 + abs(sigma * fp1);
    phi_u(5) = 1 + abs(sigma * 1);
    Phi_u{k} = phi_u;
    
    phi_l = ones(1,5);
    phi_l(1) = xp1 - abs(sigma * xp1);
    phi_l(2) = - xp2 - abs(sigma * xp2);
    phi_l(3) = - f(k) - abs(sigma * f(k));
    phi_l(4) = - fp1 - abs(sigma * fp1);
    phi_l(5) = 1 - abs(sigma * 1);
    Phi_l{k} = phi_l;
    k = k + 1;
end

%function call
tic %start clock
%[centers, generators] = CAZI(initial_thetas, H, gamma, B1, B2, Phi_u, Phi_l, max_segments, N)
[centers,generators] = PAZI(initial_thetas, H, gamma, B1, B2, Phi_u, Phi_l, max_segments, N)
toc %stop clock

%plotting the results
Plot(centers, generators, length(centers), time_var)