function BoundedIdentification(theta_c, H, sigma, gamma, y, regressor, max_segments, N)
%Function for guaranteed system identification
%Based on the paper "Bounded Error Identification of Systems With Time-Varying Parameters"

%Parameters
%theta_c - initial values
%H - initial matrix of the search zonotope
%sigma - value of the error
%gamma - vector of the expansion factors
%y - vector of the system output
%regressor - regression vector
%max_order - maximum order of the computed zonotopes
%N - number of iterations

[dimension1, dimension2] = size(theta_c); %get the number of parameters to be estimated
parameters_number = dimension1;

p = theta_c; %initialisation of the zonotope's center 
[nrows,ncolumns] = size(H); 
order = ncolumns; %extraction of the order of the zonotope
Tbest = cell(1,N); %instantiation of the matrix containing the T matrixes
vbest = cell(1,N); %instantiation of the matrix containing the v vectors
Gamma = diag(gamma); %diagonal matrix of the expansion factors
max_order = max_segments / parameters_number; %calculating the maximum order of the zonotopes

for index = 1:N
    T_set = cell(1,order + 1); %list of the matrixes T
    v_set = cell(1,order + 1); %list of vectors v
    volumes_list = []; %list of the volumes of the zonotopes for finding the minimum
    c = regressor{index};
    d = y(index); 
    [T_set,v_set,volumes_list] = intersection(order,p,d,c,H,sigma); %calling the function for the intersection of a zonotope and a strip
    [min_volume, jstar] = min(volumes_list); %get the smallest volume and its corresponding index (j*)
    T_star = T_set{jstar}; %get the T(j*)
    v_star = v_set{jstar}; %get the v(j*)
    T_new = horzcat(T_star,Gamma); %set expansion
    if order >= max_order %zonotope reduction 
        zono_matrix = horzcat(v_star,T_new);
        z = zonotope(zono_matrix);
        z_reduced = reduce(z,'girard',max_order);
        T_new = generators(z_reduced);
        v_star = center(z_reduced);
    end
    H = T_new;
    p = v_star;
    Tbest{index} = H;
    vbest{index} = p;
    [nrows,ncolumns] = size(H);
    order = ncolumns;
end

%calculations of the limits of the zonotopes at every instant k
bounds = cell(1,N);
for i = 1:N
    bounds{i} = sum(abs(Tbest{i}),2);
end

%visualisation of the parameters
for i = 1:parameters_number
    figure();
    upper = zeros(1,N); %upper_bounds
    lower = zeros(1,N); %lower_bounds
    centers = zeros(1,N); %center
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
    plot(centers,'g', 'LineWidth',1.5);
    plot(upper,'b','LineWidth',1.5);
    plot(lower,'b','LineWidth',1.5);
    xlabel('k');
end

end