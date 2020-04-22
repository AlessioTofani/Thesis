function[Delta,P] = PRadC_BS_SMA(Phi, Gamma, Sigma)
%- Use of Bisection Algorithm to obtain beta

%-------START OF CODE-----%
ns = size(Gamma,1);
nh = size(Sigma,1);
try
P = sdpvar(ns,ns,'symmetric'); % The vector in state dimensions
X = sdpvar(ns,nh,'full'); %The gain matrix
tau = sdpvar(1,1); % The factor to be maximized for decreasing the P radius
cond_temp = max(max(Gamma*Gamma'));
%------Bisection Algorithm-----%
beta_up = 1;  % Max value of beta
beta_lo = 0; % Min value of beta

beta_tol = 0.1; % 0.1

beta_wr = beta_lo; % Current working value of beta
%------Start of loop
while(beta_up-beta_lo)> beta_tol
    beta_tst = (beta_up+beta_lo)/2;
    GAx = blkvar;
    GAx(1,1) = beta_tst*P;
    GAx(1,2) = 0;
    GAx(1,3) = 0;
    GAx(1,4) = P-Phi'*X';
    GAx(2,2) = Gamma.'*Gamma;
    GAx(2,3) = 0;
    GAx(2,4)= Gamma'*P-Gamma'*Phi'*X';
    GAx(3,3) = Sigma.'*Sigma; %
    GAx(3,4) = Sigma*X';%
    GAx(4,4) = P;
    GAx = sdpvar(GAx); %make the matrix symmetric
    cond2 = ((1-beta_tst)*P)/(norm(Sigma)+cond_temp);
    % Define the optimization criterion
    crit = -tau;
    %% Solve the problem
    % -----------------
    options_sdp=sdpsettings;
    options_sdp.solver='mosek'; % sedumi %
    options_sdp.shift=1e-5; % next two lines: numerical parameters
    options_sdp.verbose=0;  % =0: suppress screenoutput, =1: allow screen output
    % LMI problem to be solved
    pblmi =  [(P>=0) , (GAx>=0) ,(cond2>=eye(ns)*tau) ,(tau>=0) ];
    
    % Solve LMI conditions
    solpb = optimize(pblmi,crit,options_sdp);
    % Check if LMI is feasible
    if solpb.problem == 1
        disp('LMIs are infeasible');
        beta_up = beta_tst;
    else
        %     disp('LMIs are FEASIBLE');
        beta_lo = beta_tst;
        beta_wr = beta_tst;
    end
end
P = value(P);
X = value(X);
tau = value(tau);
Delta = inv(P)*X;

end

