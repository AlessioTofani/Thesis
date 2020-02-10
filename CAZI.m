function CAZI(initial_thetas, H, sigma, gamma, x, regressor, max_segments, N, e, parameters, m)
%Function for guaranteed system identification
%Based on the paper "Bounded Error Identification of Systems With Time-Varying Parameters"

%Parameters
%initial_thetas - initial values
%H - initial matrix of the search zonotope
%sigma - value of the error
%gamma - vector of the expansion factors
%x - vector of the system output
%regressor - regression vector
%max_segments - maximum number of generators of the computed zonotopes
%N - number of iterations
%e - error
%parameters - real values of the parameters of the system
%m - number of measurements

[dimension1, dimension2] = size(initial_thetas); %get the number of parameters to be estimated
parameters_number = dimension1; 

p = initial_thetas; %initialization of the zonotope's center 
[nrows,ncolumns] = size(H); 
order = ncolumns; %extraction of the order of the zonotope
Tbest = cell(1,N+2); %instantiation of the matrix containing the T matrixes
vbest = cell(1,N+2); %instantiation of the matrix containing the v vectors
Gamma = diag(gamma); %diagonal matrix of the expansion factors
max_order = max_segments / parameters_number; %calculating the maximum order of the zonotopes

k = 3;
while k <= N + 2 %iteration over the number of iterations
    i = 1;
    matrix = horzcat(p,H);
    Z = zonotope(matrix);
    Zhalf = halfspace(Z); %creation of the half space representation
    half = Zhalf.halfspace; %extracing the values
    A = half.H;
    b = half.K;
    
    while i <= m %iteration over the number of measurements
        b1 = x(k) - e(k);
        b2 = x(k) - e(k);
        phi_u = ones(1,2);
        phi_u(1) = x(k-1) + abs(sigma * x(k-1));
        phi_u(2) = -x(k-2) + abs(sigma * x(k-2));
        phi_l = ones(1,2);
        phi_l(1) = x(k-1) - abs(sigma * x(k-1));
        phi_l(2) = -x(k-2) - abs(sigma * x(k-2));
        
        %first strip
        Abase = vertcat(A,phi_l);
        bbase = vertcat(b,b2);
        Aeq = []; %equality constraints
        beq = []; %equality constraints
        x0 = [1,1]; %initial guess
        lb = [0,0]; %lower bounds
        up = []; %upper bounds
        options = optimoptions('fmincon','Display','off'); 

        theta1 = fmincon(@objective, x0, Abase, bbase, Aeq, beq, lb, up, @nlcon, options, phi_u, phi_u, b1);
        delta1_star = objective(theta1,phi_u,phi_u,b1);
        %strip variables
        c1 = phi_u.';
        d1 = b1 + delta1_star/2 * norm(phi_u);
        sigma1 = delta1_star/2 * norm(phi_u);
        
        %intersection of the zonotope and the first strip
        [T_set_1,v_set_1,volumes_list_1] = intersection(order,p,d1,c1,H,sigma1);
        
        %second strip
        theta2 = fmincon(@objective, x0, Abase, bbase, Aeq, beq, lb, up, @nlcon, options, phi_l, phi_u, b2);
        delta2_star = objective(theta2,phi_l,phi_u,b2);

        %strip variables
        c2 = phi_l.';
        d2 = b2 + delta2_star/2 * norm(phi_l);
        sigma2 = delta2_star/2 * norm(phi_l);
        
        %intersection of the zonotope and the second strip
        [T_set_2,v_set_2,volumes_list_2] = intersection(order,p,d2,c2,H,sigma2);
        T_set =[T_set_1,T_set_2];
        v_set = [v_set_1,v_set_2];
        volumes_list = [volumes_list_1,volumes_list_2];
        [min_volume, jstar] = min(volumes_list); %get the smallest volume and its corresponding index (j*)
        H_new = T_set{jstar};
        v_new = v_set{jstar};
        i = i +1;
    end
    H_new = horzcat(H_new,Gamma); %set expansion
    if order >= max_order %zonotope reduction 
        zono_matrix = horzcat(v_new,H_new);
        z = zonotope(zono_matrix);
        z_reduced = reduce(z,'girard',max_order);
        H_new = generators(z_reduced);
        v_new = center(z_reduced);
    end
    H = H_new;
    p = v_new;
    Tbest{k} = H;
    vbest{k} = p;
    [nrows,ncolumns] = size(H);
    order = ncolumns;
    k = k + 1;
end

Tbest = Tbest(1,4:102); %cut the first empty values
vbest = vbest(1,4:102); %cut the first empty values
steps = N - 1; %number of steps for the graphs

%calculations of the limits of the zonotopes at every instant k
bounds = cell(1,steps);
for i = 1:steps
    bounds{i} = sum(abs(Tbest{i}),2);
end

%visualization of the parameters
for i = 1:parameters_number
    figure();
    upper = zeros(1,steps);
    lower = zeros(1,steps);
    centers = zeros(1,steps);
    for j = 1:steps
        current_bound = bounds{j};
        current_bound = current_bound(i,:);
        actual_value = vbest{j};
        center_temp = actual_value(i,:);
        centers(j) = center_temp;
        upper(j) = current_bound + center_temp;
        lower(j) = center_temp - current_bound;
    end
    hold on;
    title("θ_" + i);
    plot(centers,'g');
    plot(upper,'b');
    plot(lower,'b');
    plot(parameters{i}, 'k');
    xlabel('k');
end