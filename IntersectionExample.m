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

plotStrip(c, d, sigma, bounds); %plot the strip
plot(Z); %plot the zonotope
ylim([-0.8,0.6]);
T = []; %initialize T
T_set = cell(1,r + 1); %list of the matrixes T
v_set = cell(1,r + 1 );
volumes_list = []; %list of the volumes of the zonotopes
cc = lines; %color map for the tight strips
for j = 0:r
    if j == 0
        v = p;
        T = H;
    else
        v = p + ((d - c.' * p) / (c.'* H(:,j))) * H(:,j);
        for i = 1:r
            if i == j 
                 Tji = (sigma / (c.' * H(:,j))) * H(:,j);
            else
                Tji = H(:,i) - (c.' * H(:,i) / (c.' * H(:,j))) * H(:,j);
            end
            T = horzcat(T, Tji);
        end
    end
    T_set{j + 1} = T;
    v_set{j + 1} = v;
    intersection = zonotope([p(1,:), T(1,1), T(1,2), T(1,3); p(2,:), T(2,1), T(2,2), T(2,3)]);
    volume = det(T * T.');
    volumes_list = horzcat(volumes_list, volume);
    plot(intersection, [1 2],'color',cc(j+1,:)); %plot the tight strip
    T = []; %reinitialize T
end

[min_volume, jstar] = min(volumes_list); %get the smallest volume and its corresponding index (j*)