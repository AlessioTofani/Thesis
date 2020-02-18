function CAZI(initial_thetas, H, gamma, B1, B2, Phi_u, Phi_l, max_segments, N, parameters)
%Function perform the CAZI algorithm
%Based on the paper "Zonotope-based recursive estimation of the feasible solution 
%set for linear static systems with additive and multiplicative uncertainties"

%Parameters
%initial_thetas - initial values
%H - initial matrix of the search zonotope
%gamma - vector of the expansion factors
%B1 - y - u;
%B2 - y - u;
%Phi_u - upper bounds;
%Phi_l - lower bounds;
%max_segments - maximum number of generators of the computed zonotopes
%N - number of iterations
%parameters - real values of the parameters of the system

[dimension1, dimension2] = size(initial_thetas); %get the number of parameters to be estimated
parameters_number = dimension1; 

p = initial_thetas; %initialization of the zonotope's center 
[nrows,ncolumns] = size(H); 
order = ncolumns; %extraction of the order of the zonotope
Tbest = cell(1,N); %instantiation of the matrix containing the T matrixes
vbest = cell(1,N); %instantiation of the matrix containing the v vectors
Gamma = diag(gamma); %diagonal matrix of the expansion factors
max_order = max_segments / parameters_number; %calculating the maximum order of the zonotopes

k = 1;
while k <= N %iteration over the number of iterations
    matrix = horzcat(p,H);
    Z = zonotope(matrix);
    Zhalf = halfspace(Z); %creation of the half space representation
    half = Zhalf.halfspace; %extracing the values
    A = half.H;
    b = half.K;
    
    b1 = B1{k};
    b2 = B2{k};
    phi_u = Phi_u{k};
    phi_l = Phi_l{k};
    
    %first strip
    theta = zeros(length(phi_u),1);
    [maximum,index] = max(abs(phi_u));
    theta(index) = b1 / phi_u(index);
    f=phi_u/norm(phi_u); %objective function
    Ap=[A;-phi_u;phi_l]; %parameters for the linear programming problem
    bp=[b;-b1;b2]; %parameters for the linear programming problem
    [X,fval] = linprog(f,Ap,bp); %calling the linear progragming function
    delta1=fval-theta'*f';
    [X,fval] = linprog(-f,Ap,bp); %calling the linear progragming function
    delta2=fval+theta'*f';
    delta1_star=max([abs(delta1),abs(delta2)]); %getting the best delta

    %strip variables
    c1 = phi_u.';
    d1 = b1 + delta1_star/2 * norm(phi_u);
    sigma1 = delta1_star/2 * norm(phi_u);

    %intersection of the zonotope and the first strip
    [T_set_1,v_set_1,volumes_list_1] = intersection(order,p,d1,c1,H,sigma1);

    %second strip
    theta = zeros(length(phi_l),1);
    [maximum,index] = max(abs(phi_l));
    theta(index) = b2 / phi_l(index);
    f=phi_l/norm(phi_l); %objective function
    Ap=[A;-phi_u;phi_l]; %parameters for the linear programming problem
    bp=[b;-b1;b2]; %parameters for the linear programming problem
    [X,fval] = linprog(f,Ap,bp); %calling the linear progragming function
    delta1=fval-theta'*f';
    [X,fval] = linprog(-f,Ap,bp); %calling the linear progragming function
    delta2=fval+theta'*f';
    delta2_star=max([abs(delta1),abs(delta2)]); %getting the best delta

    %strip variables
    c2 = phi_l.';
    d2 = b2 + delta2_star/2 * norm(phi_l);
    sigma2 = delta2_star/2 * norm(phi_l);

    %intersection of the zonotope and the second strip
    [T_set_2,v_set_2,volumes_list_2] = intersection(order,p,d2,c2,H,sigma2);
    T_set =[T_set_1,T_set_2]; %merging the 2 T_sets
    v_set = [v_set_1,v_set_2]; %merging the 2 v_sets
    volumes_list = [volumes_list_1,volumes_list_2];
    [min_volume, jstar] = min(volumes_list); %get the smallest volume and its corresponding index (j*)
    H_new = T_set{jstar};
    v_new = v_set{jstar};
    
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

%calculations of the limits of the zonotopes at every instant k
bounds = cell(1,N);
for i = 1:N
    bounds{i} = sum(abs(Tbest{i}),2);
end

%visualization of the parameters
for i = 1:parameters_number
    figure();
    upper = zeros(1,N);
    lower = zeros(1,N);
    centers = zeros(1,N);
    for j = 1:N
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
    plot(centers,'g','LineWidth',1.5);
    plot(upper,'b', 'LineWidth',1.5);
    plot(lower,'b', 'LineWidth',1.5);
    para = parameters{i};
    para = para(1,1:N);
    plot(para, 'k', 'LineWidth',1.5);
    xlabel('k');
end

%alternative way to visualize the paramters with their bounds as zonotopes
figure();
cc = lines; %color map for the tight strips
for j = 1:N
    if mod(j,25) == 0 || (j == 1)
        hold on;
        zono_matrix = horzcat(vbest{j}, Tbest{j});
        z = zonotope(zono_matrix);
        plot(z, [1 2],'color',cc(j+1,:), 'LineWidth',1.5);
        plot(parameters{1}(j), parameters{2}(j),'color',cc(j+1,:), 'Marker', '*', 'MarkerSize', 10);
    end
end
plot(parameters{1}, parameters{2}, 'k', 'LineWidth',1.5);


end