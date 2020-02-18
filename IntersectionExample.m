%Bounded Error Identification of Systems With Time-Varying Parameters

%Example of the intersection of a Zonotope and a Strip 
%Data of the problem
Z = zonotope([0,0.2812,0.1968,0.4235;0,0.0186,-0.2063,-0.2267]); 
p = [0;0];
H = [0.2812,0.1968,0.4235;0.0186,-0.2063,-0.2267];
d = 0;
c = [1;-1];
sigma = 0.1;
r = 3; %order of the zonotope

figure; hold on;
grid on;
set(gca,'FontSize',18);
bounds = [-1,1];
plotStrip(c, d, sigma, bounds); %plot the strip
cc = lines; %color map for the tight strips
[T_set,v_set,volumes_list] = intersection(r,p,d,c,H,sigma);
for i = 1:4
    v = v_set{i};
    T = T_set{i};
    inter = zonotope([v(1,:), T(1,1), T(1,2), T(1,3); v(2,:), T(2,1), T(2,2), T(2,3)]);
    plot(inter, [1 2],'color',cc(i+1,:), 'LineWidth',2); %plot the tight strip
end
plot(Z, [1 2], 'LineWidth',2, 'color', 'b'); %plot the zonotope
[min_volume, jstar] = min(volumes_list); %get the smallest volume and its corresponding index (j*)
