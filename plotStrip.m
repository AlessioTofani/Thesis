function plotStrip(c,d,sigma,x)
%Function that plot a strip
%Parameters:
%c - vector c for constructing the strip
%d - value d for constructing the strip
%sigma - interval
%x -bounds of the plotting
%the strip is defined as: 
%S = {x:|c.T * x - d| <= sigma}

%General formula of a straight line:
%y = mx + q

m = - (c(1,:) / c(2,:));
q = d / c(2,:);
hold on;
y = m*x + q + (sigma / c(2,:)); 
plot(x,y,'r:');
y = m*x + q - (sigma / c(2,:)); 
plot(x,y,'r:');

end