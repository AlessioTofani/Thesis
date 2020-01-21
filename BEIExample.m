%Example of Bounded Error Identification

%Data of the problem
theta_c = [0.9;-0.8;-1;-0.8;0]; %initial parameters vector
H = 2*ones(5); %initial generators
sigma = 0.05; 
gamma = [0.01,0.01,0.01,0.01,0.05]; %expansion factors
max_order = 40; %max order for the zonotopes

k = 100; %number of iterations
u = zeros(1,k+2); %instantiation of the input vector 
u = -5 + (5+5)*rand(1,k+2); %random bounded input |u| < 5

%initial conditions
u(1) = 0; 
u(2) = 0;
%generation of a bounded error
e = -0.05 + (0.05+0.05) *rand(1); 

%generate the perturbation
omega_0 = 0;
omega = zeros(1,k+2);
omega(1) = omega_0;
for i = 2: k + 2
    variation = -0.05 + (0.05+0.05) *rand(1);
    omega(i) = omega(i-1) + variation;
end

%calculate the output
y = zeros(1,k+2); %instantiation of the output vector
for i = 3:k + 2
    y(i) = 0.9*y(:,i-1) -0.8*y(:,i-2) - u(i) -0.8*u(:,i-1) + omega(i) + e; 
end

%calculate the regression vector
regressor = cell(1,k+2);
for i = 3:k+2
    regressor{i} = [y(:,i-1);y(:,i-2);u(i);u(i-1);1];
end

%function call
BoundedIdentification(theta_c, H, sigma, gamma, y, regressor, max_order, k);