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
time_var = 1; %specify if use also the alternative visualization

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
    [nrows,ncolumns] = size(T_new); 
    new_order = ncolumns; %extraction of the order of the zonotope
    if new_order >= max_order %zonotope reduction 
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

%visualization of the parameters
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
    set(gca,'FontSize',18);
    plot(centers,'g', 'LineWidth',1.5);
    plot(upper,'b','LineWidth',1.5);
    plot(lower,'b','LineWidth',1.5);
    xlabel('k');
    ylabel("θ_" + i);
end

%alternative way to visualize the paramters with their bounds as zonotopes
if time_var == 1
    iteration_count = 20; %specify every many iterations to display the zonotope
    for i = 1:parameters_number
        figure();
        cc = lines; %color map for the tight strips
        centers = zeros(1,N); %center
        for j = 1:N
            actual_value = vbest{j};
            center_temp = actual_value(i,:);
            centers(j) = center_temp;
            if mod(j,iteration_count) == 0 || (j == 1) 
                hold on;
                v = zeros(parameters_number,1);
                v(1) = j;
                v(2) = center_temp(1);
                zono_matrix = horzcat(v, Tbest{j});
                z = zonotope(zono_matrix);
                plot(z, [1 2],'color',cc(j+1,:), 'LineWidth',1.5);
            end
        end
        plot(centers,'g', 'LineWidth',1.5);
        xlabel('k');
        ylabel("θ_" + i);
    end
end

end