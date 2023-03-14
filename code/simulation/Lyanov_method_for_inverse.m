% test the static output feedbak by Lyapunov
rng(121)
A = eye(3);
B = rand(3,2);
D = eye(2);
M = [A B; B' D];
%
eig(M)
%% networked control system exercise
% Define the system
A = [0.9053 0.0928; 0.0098 0.9512];
% set yalmip options
ops = sdpsettings('solver','mosek');
ops = sdpsettings(ops,'verbose',1);
% Define unknowns, parameters and constraints
P = sdpvar(2,2); % Unknown 2x2 symmetric matrix
Q = 1 * eye(2,2);
CONS1=[A'*P*A-P<=-Q]; % Constraint 1
CONS2=[P>=0]; % Constraint 2
CONS =[CONS1 ,CONS2]; % Combine all constraints
% Solving for P
infosol = optimize(CONS,[],ops);
Psol = double(P) % Converts to standard matrix format
% Check if the solution is OK
infosol.info
% Double-check the solution verifies the constraints
% Psol must be >0. Chech that all eigenvalues of Psol are > 0
eigP=eig(Psol)
% the matrix -(A'*Psol*A-Psol+Q) must be >0. Chech that all its eigenvalues are > 0
eigcons0=eig(-(A'*Psol*A-Psol+Q))
%% static output feedback: CT W method
A = [zeros(2,4); 2 0 1 2; 0 2 2 1];
B = [eye(2); zeros(2,2)];
C = [4 -3 3 0; -2 2 -1 2];

ops = sdpsettings('solver','mosek');
ops = sdpsettings(ops,'verbose',1);
% Define unknowns, parameters and constraints
P = sdpvar(4,4); % Unknown 2x2 symmetric matrix
M = sdpvar(2,2,'full');
N = sdpvar(2,2,'full');
Q = 1 * eye(4);
CONS=[A*P + P*A' - B*N*C - C'*N'*B'<= -Q]; % Constraint 1
CONS=[CONS, P>=0]; % Constraint 2
CONS=[CONS, M*C == C*P]; % Constraint 3
CONS = [CONS, N(1,2)==0, N(2,1)==0, M(2,1)==0, M(1,2)==0];
% Solving for P
infosol = optimize(CONS,[],ops);
Psol = double(P) % Converts to standard matrix format

eigP=eig(Psol)
F = double(N)*inv(double(M))
eig(A-B*F*C)
%% static output feedback: DT P method
A = [2 0.3 2; 1 0 1; 0.3 0.6 -0.6];
B = [1 0; 0 1; 1 0];
C = [1 1 0];

ops = sdpsettings('solver','mosek');
ops = sdpsettings(ops,'verbose',1);
% Define unknowns, parameters and constraints
P = sdpvar(3,3); % Unknown 2x2 symmetric matrix
M = sdpvar(2,2,'full');
N = sdpvar(2,1,'full');
Q = 1 * eye(6);
CONS=[ [ P  A'*P-C'*N'*B';
    P*A-B*N*C P] >= 0]; % Constraint 1
CONS=[CONS, P>=0]; % Constraint 2
CONS=[CONS, B*M == P*B]; % Constraint 3
% CONS = [CONS, N(1,2)==0, N(2,1)==0, M(2,1)==0, M(1,2)==0];
% Solving for P
infosol = optimize(CONS,[],ops);

Psol = double(P) % Converts to standard matrix format

eigP=eig(Psol)
F = inv(double(M))*double(N)
abs(eig(A-B*F*C))
infosol.info
%% W method
P = sdpvar(3,3); % Unknown 2x2 symmetric matrix
M = sdpvar(1,1,'full');
N = sdpvar(2,1,'full');
Q = 1 * eye(6);
CONS=[ [ P P*A' - C'*N'*B';
   A*P-B*N*C  P] >= Q]; % Constraint 1
CONS=[CONS, P>=0.1*eye(3)]; % Constraint 2
CONS=[CONS, M*C == C*P]; % Constraint 3
% CONS = [CONS, N(1,2)==0, N(2,1)==0, M(2,1)==0, M(1,2)==0];
% Solving for P
infosol = optimize(CONS,[],ops);

Psol = double(P) % Converts to standard matrix format

eigP=eig(Psol)
F = double(N)*inv(double(M))
abs(eig(A-B*F*C))
infosol.info


%% static output feedback: DT W method C not full row rank
A = [2 0.3 2; 1 0 1; 0.3 0.6 -0.6];
B = [1 0; 0 1; 1 0];
C_ori = rand(4,3);%[eye(3);0.5 0 0];
T = pinv(C_ori);
C = T*C_ori; 

ops = sdpsettings('solver','mosek');
ops = sdpsettings(ops,'verbose',1);
W = sdpvar(3,3); % Unknown 2x2 symmetric matrix
M = sdpvar(3,3,'full');
N = sdpvar(2,3,'full');
Q = 1 * eye(6);
CONS=[ [ W W*A' - C'*N'*B';
   A*W-B*N*C  W] >= Q]; % Constraint 1
CONS=[CONS, W>=0]; % Constraint 2
CONS=[CONS, M*C == C*W]; % Constraint 3
CONS = [CONS, N(1,1:3)==0];%N(2,1)==0, M(2,1)==0, M(1,2)==0];
% Solving for P
infosol = optimize(CONS,[],ops);

Psol = double(W) % Converts to standard matrix format

eigP=eig(Psol)
F = double(N)*pinv(double(M))*T
abs(eig(A-B*F*C_ori))
infosol.info