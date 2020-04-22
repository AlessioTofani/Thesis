function Plot(centers, generators, N, time_var)
%Function for plotting the results of the parameters identification

%Parameters
%centers - cell array containing the centers of the zonotopes
%generators - cell array containing the generators of the zonotopes
%N - number of iterations
%time_var - flag for specifing if using also the alternative way to
%visualize the results

bounds = cell(1,N);
for i = 1:N
    bounds{i} = sum(abs(generators{i}),2);
end

parameters_number = size(centers{1},1);

%plot for the N iterations
for i = 1:parameters_number
    figure();
    upper = zeros(1,N);
    lower = zeros(1,N);
    center = zeros(1,N);
    for j = 1:N
        current_bound = bounds{j};
        current_bound = current_bound(i,:);
        actual_value = centers{j};
        center_temp = actual_value(i,:);
        center(j) = center_temp;
        upper(j) = current_bound + center_temp;
        lower(j) = center_temp - current_bound;
    end
    hold on;
    plot(center,'g','LineWidth',1.5);
    plot(upper,'b', 'LineWidth',1.5);
    plot(lower,'b', 'LineWidth',1.5);
    xlabel('k');
    ylabel("θ_" + i);
end

%alternative way to visualize the parameters with their bounds as zonotopes
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
