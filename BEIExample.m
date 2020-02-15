%Example of Bounded Error Identification
clc;
clear;
%Data of the problem
theta_c = [0.9;-0.8;-1;-0.8;0]; %initial parameters vector
H = 0.2*ones(5); %initial generators
sigma = 0.05; 
gamma = [0.01,0.01,0.01,0.01,0.05]; %expansion factors
max_segments = 40; %max number of segments forming the zonotopes

N = 100; %number of iterations
u = zeros(1,N); %instantiation of the input vector 
u = -5 + (5+5)*rand(1,N); %random bounded input |u| < 5

%initial conditions
u(1) = 0; 

%generation of the bounded error
e = zeros(1,N);
for i = 1: N
    e(i) = -0.05 + (0.05+0.05) * rand(1); 
end

%genearation of the perturbation
omega = zeros(1,N);
for i = 1: N
    variation = -0.05 + (0.05+0.05) * rand(1);
    omegaprec = 0;
    if (i-1) > 0
        omegaprec = omega(i-1);
    end
    omega(i) = omegaprec + variation;
end
omega(1) = 0;

%calculation of the output vector
y = zeros(1,N); %instantiation of the output vector
for i = 1:N
    xp1 = 0; %initial condition
    xp2 = 0; %initial condition
    up1 = 0; %initial condition
    if (i-2) > 0
        xp2 = y(i-2);
    end
    if (i-1) > 0
        xp1 = y(i-1);
        up1 = u(i-1);
    end
    y(i) = 0.9*xp1 -0.8*xp2 - u(i) -0.8*up1 + omega(i) + e(i); 
end

%calculation of the regression vector
regressor = cell(1,N); %instantiation of the regresion vector
for i = 1:N
    yp1 = y(2);
    yp2 = y(1);
    up1 = 0;
    if (i-2) > 0
        yp2 = y(i-2);
    end
    if (i-1) > 0
        yp1 = y(i-1);
        up1 = u(i-1);
    end
    regressor{i} = [yp1;yp2;u(i);up1;1];
end

%function call
tic %start clock
BoundedIdentification(theta_c, H, sigma, gamma, y, regressor, max_segments, N);
toc %stop clock