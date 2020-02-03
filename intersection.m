function [T_set,v_set,volumes_list] = intersection(order,p,d,c,H,sigma)
%Function that calculates the intersection between a zonotope and a strip
%The zonotope is given in the form Z = p + H*B^r where + stands for 
%Minkowski Sum
%The strip is given in the form S = {x: |c.T* x - d| <= sigma}
%Parameters:
%order - order of the considered zonotope
%p - centre of the zonotope
%d - value d for constructing the strip
%c - vector c for constructing the strip
%H - matrix H (generators of the zonotope)
%sigma - value of sigma

%output
%T_set - set of the matrices T calculated
%v_set - set of the vector v calculated
%volumes_list - list of the computed volumes

volumes_list = []; %list of the volumes of the zonotopes
for j = 0:order
    T = []; %initialise T
        if j > 0 
            ctHj = abs(c.' * H(:,j));
        end
        if (j >= 1 & j <= order) & ctHj ~= 0 
            v = p + ((d - c.' * p) / (c.' * H(:,j))) * H(:,j);
            for i = 1:order
                if i == j 
                    Tji = (sigma / (c.' * H(:,j))) * H(:,j);
                else
                    Tji = H(:,i) - (c.' * H(:,i) / (c.' * H(:,j))) * H(:,j);
                end
                T = horzcat(T, Tji);
            end
        else
            v = p;
            T = H;  
        end
        T_set{j + 1} = T;
        v_set{j + 1} = v;
        volume = det(T * T.');
        volumes_list = horzcat(volumes_list, volume);
end