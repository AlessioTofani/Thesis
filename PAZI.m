function [vbest,Tbest] = PAZI(initial_thetas, H, gamma, B1, B2, Phi_u, Phi_l, max_segments, N)
%Function perform the PAZI algorithm
%Based on the paper "Zonotope-based recursive estimation of the feasible solution 
%set for linear static systems with additive and multiplicative uncertainties"

%Parameters
%initial_thetas - initial values
%H - initial matrix of the search zonotope
%gamma - vector of the expansion factors
%B1 - y - u;
%B2 - y - u;
%Phi_u - upper bounds;
%Phi_l - lower bounds;
%max_segments - maximum number of generators of the computed zonotopes
%N - number of iterations

[dimension1, dimension2] = size(initial_thetas); %get the number of parameters to be estimated
parameters_number = dimension1; 
batch_dimension = 4;
p = initial_thetas; %initialization of the zonotope's center 
[nrows,ncolumns] = size(H); 
order = ncolumns; %extraction of the order of the zonotope
%error message
if  mod(N,batch_dimension) > 0
    display('ERROR FROM PAZI: Number of iterations not divisible by batch dimension');
    return
end

vbest = cell(1,N / batch_dimension); %instantiation of the matrix containing the v vectors
Tbest = cell(1,N / batch_dimension); %instantiation of the matrix containing the T matrixes
Gamma = diag(gamma); %diagonal matrix of the expansion factors
max_order = max_segments / parameters_number; %calculating the maximum order of the zonotopes
vbest{1} = initial_thetas;
Tbest{1} = H;
nosol = 0;
i_plot = 0; %index for plotting the results
m_index = 0; %index of the batch
batch_C = [];
batch_D = [];
batch_Sigma = [];

for k=1:N %iteration over the number of iterations
    fprintf('ITERATION %d\n',k);
    nosol = 0;
    matrix = horzcat(p,H);
    Z = zonotope(matrix);
    Zhalf = halfspace(Z); %creation of the half space representation
    half = Zhalf.halfspace; %extracing the values
    A = half.H;
    b = half.K;
   
    b1 = B1{k};
    b2 = B2{k};
    phi_u = Phi_u{k};
    phi_l = Phi_l{k};

    %first strip
    theta = zeros(length(phi_u),1);
    [maximum,index] = max(abs(phi_u));
    theta(index) = b1 / phi_u(index);
    f=phi_u/norm(phi_u); %objective function
    Ap=[A;-phi_u;phi_l]; %parameters for the linear programming problem
    bp=[b;-b1;b2]; %parameters for the linear programming problem
    [X,fval] = linprog(f,Ap,bp); %calling the linear progragming function
    delta1=fval-theta'*f';
    [X,fval] = linprog(-f,Ap,bp); %calling the linear progragming function
    delta2=fval+theta'*f';
    delta1_star=max([abs(delta1),abs(delta2)]); %getting the best delta

    %strip variables
    c1 = phi_u.';
    d1 = b1 + delta1_star/2 * norm(phi_u);
    sigma1 = delta1_star/2 * norm(phi_u);

    batch_C = horzcat(batch_C,c1);
    batch_D = vertcat(batch_D,d1);
    batch_Sigma = horzcat(batch_Sigma,sigma1);

    %second strip
    theta = zeros(length(phi_l),1);
    [maximum,index] = max(abs(phi_l));
    theta(index) = b2 / phi_l(index);
    f=phi_l/norm(phi_l); %objective function
    Ap=[A;-phi_u;phi_l]; %parameters for the linear programming problem
    bp=[b;-b1;b2]; %parameters for the linear programming problem
    [X,fval] = linprog(f,Ap,bp); %calling the linear progragming function
    delta1=fval-theta'*f';
    [X,fval] = linprog(-f,Ap,bp); %calling the linear progragming function
    delta2=fval+theta'*f';
    delta2_star=max([abs(delta1),abs(delta2)]); %getting the best delta

    %strip variables
    c2 = phi_l.';
    d2 = b2 + delta2_star/2 * norm(phi_l);
    sigma2 = delta2_star/2 * norm(phi_l);
    
    batch_C = horzcat(batch_C,c2);
    batch_D = vertcat(batch_D,d2);
    batch_Sigma = horzcat(batch_Sigma,sigma2);
    
    m_index = m_index + 1;
    if m_index < batch_dimension
       continue; 
    end

    PHI = batch_C;
    C = PHI';
    D = batch_D;
    sig = batch_Sigma;
    SIG = diag(sig);
    
    %Solving the LMI problem
    [Lambda,P] = PRadC_BS_SMA(C, Gamma, SIG);
    
    if nosol == 0
         v_new = p+Lambda*(D-PHI'*p); %center of the candidate zonotope    
         H_new = [(eye(parameters_number)-Lambda*PHI')*H Lambda*SIG]; %generators of the candidate zonotope
         logsol(k) = 1;
    end

    H_new = horzcat(H_new,Gamma); %set expansion
    if order >= max_order %zonotope reduction 
        zono_matrix = horzcat(v_new,H_new);
        z = zonotope(zono_matrix);
        z_reduced = reduce(z,'girard',max_order);
        H_new = generators(z_reduced);
        v_new = center(z_reduced);
    end
    i_plot = i_plot + 1;
    H = H_new;
    p = v_new;
    Tbest{i_plot + 1} = H;
    vbest{i_plot + 1} = p;
    [nrows,ncolumns] = size(H);
    order = ncolumns;

    %reset the batches
    batch_C = [];
    batch_D = [];
    batch_Sigma = [];
    m_index = 0;
end
end