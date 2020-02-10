%Example of CAZI 

%Data of the problem
N  = 100; %number of iterations
max_order = 10; %max order for the zonotope
sigma = 0.03;
initial_thetas = [0.9;0.8];
theta_1 = zeros(1, N); %initialisation of theta_1
theta_1(1) = initial_thetas(1,:);
theta_2 = zeros(1,N); %initialisation of theta_2
theta_2(1) = initial_thetas(2,:);
gamma = [0.05;0.05]; %expansion vector
m = 1; %number of measurements

initial_center = [0.9;0.8]; %initial search space center
H = [-0.2,0.1;0.1,0.2]; %initial generators

%generation of the time-varying parameter theta1
for i = 2:N + 2
    theta_1(i) = theta_1(i-1) + 0.05 - i * 0.0008;
end

%generation of the time-varying parameter theta2
for i = 2:N + 2
    theta_2(i) = theta_2(i-1) + 0.05 - i * 0.001;
end

%generation of the bounded error (check if it's enough)
e = zeros(1,N+2);
for i = 3: N + 2
    e(i) = -1 + (1+1) * rand(1);
end

x = zeros(1,N+2); %instantiation of the output vector
%initial conditions
x(1) = 0.1;
x(2) = 0.1;

%calculate the output
for i = 3:N + 2
    x(i) = theta_1(i) * x(i-1) - theta_2(i) * x(i-2) +  e(i); 
end

%calculate the regression vector
regressor = cell(1,N + 2);
for i = 3:N+2
    regressor{i} = [x(:,i-1);x(:,i-2)];
end

%function call
tic %start clock
CAZI(initial_thetas, H, sigma, gamma, x , max_order, N, e, theta_1, theta_2, m)
toc %stop clock