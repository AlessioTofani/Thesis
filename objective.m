function result = objective(theta,phi_calc,phi_cons,b1)
%objective function that has to be maximized
% ((theta_hat - theta_tilde)^T * phi) / norm(phi)
    theta_tilde = pinv(phi_calc) * b1;
    result = - abs((([theta(1) - theta_tilde(1); theta(2) - theta_tilde(2)]).' * phi_calc.') / norm(phi_calc));
end