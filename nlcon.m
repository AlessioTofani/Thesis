function [c,ceq] = nlcon(theta,phi_calc,phi_cons,b12)
    %non linear constraint for the maximization problem
    c = -(phi_cons * [theta(1);theta(2)] - b12); %constraint phi.'*thetai >= bi
    ceq = 0;
end