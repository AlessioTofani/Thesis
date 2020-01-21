function plotStrip(c, d, sigma, bounds)
%Function that plot a strip
%Parameters:
%c
%d
%sigma
%the strip is defined as: 
%S = {x:|c.T * x - d| <= sigma}

%General formula of a straight line:
%y = mx + q

x = bounds;
q = sigma;
for i = 1:size(c)
    m = c(i);
    m = abs(m);
    y=m*x+q;
    hold on;
    plot(x,y,'r:');
    y=m*x - q;
    plot(x,y,'r:');
end
end
