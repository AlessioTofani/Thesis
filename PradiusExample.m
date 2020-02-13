%This example shows an alternative way ob building candidate zonotopes

%Initial zonotope
Z = zonotope([0.1,0.1,0.2,0.3;-0.5,0.3,0.2,0.1]);
p0 = center(Z); %zonotope center
H0 = generators(Z); %zonotope generators

%plot the initial zonotope
plot(Z,[1 2],'LineWidth',2);
hold on;
set(gca,'FontSize',18);

%first strip
d1 = -0.1163;
c1 = [5;1];
sigma1 = 0.2;
bounds = [-0.6,0.6]; %limits to be displayed
plotStrip(c1,d1,sigma1,bounds);

%second strip
d2 = -0.2935;
c2 = [-4;1];
sigma2 = 0.2;
plotStrip(c2,d2,sigma2,bounds);

%third strip
d3 = -0.6928;
c3 = [1;2];
sigma3 = 0.3;
plotStrip(c3,d3,sigma3,bounds);

%polyhedron
%S = {theta : |Phi.'*theta - D| <= sig}
Phi = [5,-4,1;1,1,2];
C = Phi';
D = [-0.1163;-0.2935;-0.6928];
sig = [0.2;0.2;0.3];
SIG = diag(sig);

beta=0.01;
epsilon=sig'*sig; 
n_lmi=9; %dimensions of the LMI matrix
P=[1 0;0 1]; %user defined matrix P
Gamma=zeros(2,2); %Expansion matrix

%LMI Problem
cvx_begin
    variable X(2,3);
    variable tau;
    maximize tau
    subject to
        (1-beta)*P/epsilon-tau*eye(2) == semidefinite(2);
        tau>=0;
    A11 = blkdiag(beta*P,Gamma'*Gamma,SIG'*SIG);
    A12 = [P-C'*X';Gamma'*P-Gamma'*C'*X'; SIG*X'];
    A21 = A12';
    A22 = P;
    f_matrix = [A11 A12; A21 P];
    f_matrix == semidefinite(n_lmi);
cvx_end

Lambda = inv(P)*X;

%Resulting zonotope
pz = p0+Lambda*(D-Phi'*p0); %center of the candidate zonotope    
Hz = [(eye(2)-Lambda*Phi')*H0 Lambda*SIG]; %generators of the candidate zonotope
zono_matrix = horzcat(pz,Hz);
Z2 = zonotope(zono_matrix);
plot(Z2,[1 2],'g','LineWidth',2);
ylim([-1.5 0.5]) %limits for the y axis
grid on;