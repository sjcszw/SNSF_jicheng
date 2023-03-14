% polydome_hankel_prediction_mystruct
clc
clear
rng(353) % 72:N1,20,N2, 0
%%
% system definition
% for i = 1:1000
% A = 2*rand(3,3); eig(A)
% B = 2*rand(3,2);
% C = 2*rand(2,3);
n1 = 20; n2 = 0;
A = round(n1*rand(3,3)-n2)/10; eig(A)
B = round(n1*rand(3,2)-n2)/10; %
C = round(n1*rand(2,3)-n2)/10;
D =  0*round(n1*rand(2,2)-n2)/10;
Ob = obsv(A,C);Co = ctrb(A,B);
% rank(Ob) 
% rank(Co)
nu = size(B, 2);
nx = size(A, 1);
ny = size(C, 1);
fx = @(x,u) A*x + B*u; % state dynamics
fy = @(x) C*x; % output dynamics
a = tzero(A,B,C,zeros(ny,nu),eye(nx))
% if rank(B) > rank(C*B)
%     1
% end
% end
%%
N_ini = 5;
N_pred = 1;
T = 50; % T >= (nu + 1)*(N_ini + N_pred + nx) - 1
T_tot = 100; % total steps for simulation

% Simulate a trajectory by random input
x_cl = zeros(nx,T_tot); y_cl = zeros(ny,T_tot);
Q = eye(nx);R = eye(nu);
[K,~,~] = dlqr(A,B,Q,R);
K = -K;
u_cl = 3*rand(nu,T_tot)-1.5;
x_cl(:, 1) = 0;
y_cl(:, 1) = C*x_cl(:, 1);
for t = 1:T_tot-1
    u_cl(:,t) = u_cl(:,t) + K*x_cl(:,t);
    x_cl(:, t+1) = fx(x_cl(:,t), u_cl(:,t));
    y_cl(:, t+1) = fy(x_cl(:,t+1));
end

%%

h=figure;
hold on
% set(h,'Units','normalized','Position',[0 0 1 .5]); 
yyaxis left
plot(u_cl(1,:),'b','LineWidth',1)
yyaxis right
% ylim([5 30])
%     plot(tExp,exp.setpoint_winter.value/10.0,'r','LineWidth',1)
plot(C(1,:)*x_cl,'r','LineWidth',1)
plot(y_cl(1,:),'g','LineWidth',1)
% plot(time,exp.supply_temp/10.0,'y','LineWidth',1.5)
% title('Meteorological data','FontSize',18)
legend({'u','x1','y'},'FontSize',18)
% ylabel('[°C]','FontSize',18)

%%

% Define the weights for g vector: L2 norm or adaptive term
% diag_vector = 0.02:- 0.018/(T-N_ini-N_pred):0.002;
% diag_vector = 2:- 0.1/(T-N_ini-N_pred):1.9;
diag_vector = 0.02:- 0.001/(T-N_ini-N_pred):0.019;
Qg = diag(diag_vector);

% build the Henkal matrix 
disp('Computing Hankel...')
N_total = N_ini + N_pred;
Hankel_col = T - N_total + 1;
Hankel_U = compute_Hankel_matrix(u_cl(:,1:T),nu,N_total, Hankel_col);
% Hankel_W_elec = compute_Hankel_matrix(exp_elec_iden.w,nw,N_total, Hankel_col);

N_total = N_ini + N_pred + 1;
Hankel_col = T - N_total + 2;
Hankel_Y = compute_Hankel_matrix(y_cl(:,1:T+1),ny,N_total, Hankel_col);
Hankel_U_past = Hankel_U(1:nu*N_ini, :);
Hankel_U_future = Hankel_U(nu*N_ini+1:end, :);
Hankel_Y_past = Hankel_Y(1:ny*(N_ini+1), :);
Hankel_Y_future = Hankel_Y(ny*(N_ini+1)+1:end, :);
%%
N_total = N_ini + N_pred + nx;
Hankel_col = T - N_total +1 ;

% Need: u is persistently exciting of order is N_ini + N_pred + nx
Hankel_UW_check = compute_Hankel_matrix(u_cl(:,1:T), nu, N_total, Hankel_col);
% Hankel_UW_check = compute_Hankel_matrix([exp_elec_iden.u; exp_elec_iden.w], nu+nw, N_total, Hankel_col);

% Check if the past u trajectory is persistently exciting of order N_ini + N_pred + nx
if rank(Hankel_UW_check)== nu*(N_ini + N_pred + nx) 
    disp('Hankel rank is ok')
else
    disp('Exciting of order of Hu is samller than N_ini + N_pred + nx')
end

%% check AS
tem = pinv([Hankel_U_past;Hankel_Y]); % Moore-Penrose pseudoinverse
U1 = Hankel_U_future*tem(:,1:nu*N_ini);
Y1 = Hankel_U_future*tem(:,nu*N_ini+1:end);
% if N_pred==1
U_check = [zeros(N_ini*nu-nu,nu),eye(N_ini*nu-nu)];
U_check = [U_check; U1];

abs(eig(U_check))
%%
tem3 = [Hankel_U_past;Hankel_Y]; % Hankel for u_init, y_init, y_{t,t+1}
[U,S,V] = svd(tem3); % SVD to get different G
% [Uf,Sf,Vf] = svd(Hankel_U_future); % SVD 
G = V*[inv(S(1:15,1:15)) rand(15,9);rand(30,24)]*U';
% G = V*[inv(S(1:15,1:15)) zeros(15,9);10000*rand(30,24)]*U'; % U1 are almost the same
% G = V*[inv(S(1:15,1:15)) -0.1+0.2*rand(15,9);zeros(30,24)]*U';

% if N_pred==1 build the dynamics matrix
U1 = Hankel_U_future*G(:,1:nu*N_ini);
U_check = [zeros(N_ini*nu-nu,nu),eye(N_ini*nu-nu)];
U_check = [U_check; U1];
abs(eig(U_check))

%% Use matrix addition and multiplication to get U_check
tem3 = [Hankel_U_past;Hankel_Y]; % Hankel for u_init, y_init, y_{t,t+1}
[U,S,V] = svd(tem3); % SVD to get different G
G = V*[inv(S(1:15,1:15)) ones(15,9);ones(30,24)]*U';
U1 = Hankel_U_future*G(:,1:nu*N_ini);
U_check = [zeros(N_ini*nu-nu,nu),eye(N_ini*nu-nu)];
U_check = [U_check; U1];

G_bar = V*[inv(S(1:15,1:15)) zeros(15,9); zeros(30,24)]*U';
G_var = V*[zeros(15,15) ones(15,9);ones(30,24)]*U';
% [0 I;0] + [0; I]*Hu*G_bar*[I;0;0]
U_bar = [zeros(N_ini*nu-nu,nu),eye(N_ini*nu-nu); zeros(nu,N_ini*nu)] + ...
    [zeros(N_ini*nu-nu,nu); eye(nu)]*Hankel_U_future*G_bar*[eye(N_ini*nu);  zeros((N_ini+2)*ny,N_ini*nu)];

% U_var = [0; I]*Hu*Gvar*[I;0;0]
U_var = [zeros(N_ini*nu-nu,nu); eye(nu)]*Hankel_U_future*G_var*[eye(N_ini*nu);  zeros((N_ini+2)*ny,N_ini*nu)];

%%
Au = U_bar;
Bu = [zeros(N_ini*nu-nu,nu); eye(nu)]*Hankel_U_future*V;
Cu_ori = U'*[eye(N_ini*nu);  zeros((N_ini+2)*ny,N_ini*nu)];

%% only last several columns to change
T_1 = [zeros(9,15),eye(9)];
Cu = T_1*Cu_ori; 
rank(Cu)
ops = sdpsettings('solver','mosek');
ops = sdpsettings(ops,'verbose',1);
W = sdpvar(10,10); % Unknown 2x2 symmetric matrix
M = sdpvar(9,9,'full');
N = sdpvar(45,9,'full');
Q = 1 * eye(20);
CONS=[ [ W W*Au'- Cu'*N'*Bu';
   Au*W-Bu*N*Cu  W] >= 0]; % Constraint 1
CONS=[CONS, W>=0]; % Constraint 2
CONS=[CONS, M*Cu == Cu*W]; % Constraint 3
% CONS = [CONS, N(1:15,:)==0];
% Solving for P
infosol = optimize(CONS,[],ops);

Psol = double(W) % Converts to standard matrix format

eigP=eig(Psol)
F = double(N)*pinv(double(M))*T_1
abs(eig(Au-Bu*F*Cu_ori))
infosol.info

% check generalized inverse
G_bar = V*[inv(S(1:15,1:15)) zeros(15,9); zeros(30,24)]*U';
G_var = V*(-F)*U';
G = G_bar + G_var;
% check inverse condition;check A_UIE
% tem3 - tem3*(G_bar+G_var)*tem3
% U1 = Hankel_U_future*G(:,1:nu*N_ini);
% U_check = [zeros(N_ini*nu-nu,nu),eye(N_ini*nu-nu)];
% U_check = [U_check; U1];
% U_check - (Au-Bu*F*Cu_ori)
%% prediction input
% y_pred = Hg_pred*[y_ini;u_ini;u_pred;w_ini;w_pred];
% Hg_pred = Hg_problem_u(Hankel_U,Hankel_Y,nu,ny,N_ini,Qg);
%
num = size(u_cl,2);
t = T+1;
% u_pred(:,1:T) = u_cl(:,1:T);
% u_pred_basic(:,1:T) = u_cl(:,1:T);
% u_pred(:,1:T) = 10*ones(nu,T)-1;
u_pred_basic(:,1:T) = 0*ones(nu,T)-1;

while(t<=num-1)
    if t+N_pred<=num
%         ui = reshape(u_cl(:,t-N_ini:t-1), [N_ini*nu,1]);
%         ui = reshape(u_pred(:,t-N_ini:t-1), [N_ini*nu,1]);
%         wip = reshape(w_cl(:,t-N_ini:t+N_pred-1), [(N_ini+N_pred)*nw,1]);
        yi = reshape(y_cl(:,t-N_ini:t), [(N_ini+1)*ny,1]); 
        yp = reshape(y_cl(:,t+1:t+N_pred), [N_pred*ny,1]); 
%         u_pred(:,t:t+N_pred-1) = reshape(Hg_pred*[yi;ui;yp],[nu, N_pred]);
%         y_pred(t+1:t+N_pred) = Hg_pred*[yi;uip;wip];
        
        ui = reshape(u_pred_basic(:,t-N_ini:t-1), [N_ini*nu,1]);
%         g = [Hankel_U_past;Hankel_Y]\[ui;yi;yp];
%         g = tem*[ui;yi;yp];
        g = G*[ui;yi;yp];
    %     [Hankel_U;Hankel_y_past]*g-[u_ini_pred;y_ini] % Almost 0
        u_pred_basic(:,t:t+N_pred-1) = reshape(Hankel_U_future*g,[nu, N_pred]);

        t = t+N_pred;
    else
        break;
    end
end
%%
num = 2;
figure()
hold on; grid on;
start = T+1; last = size(u_pred_basic,2);
plot( u_cl(num,start:last), 'Color',[0.4660 0.6740 0.188],'LineWidth',3)
plot( u_pred_basic(num,start:last), '*r--','LineWidth',1)
% plot( u_pred_basic(num,start:last), '*r--','LineWidth',1)
% plot( u_pred_basic(num,start:last), 'Color',[0.8500 0.3250 0.0980],'LineWidth',2,'LineStyle','--')
scatter( [(start:N_pred:last)-start+1],u_cl(num,start:N_pred:last), 100,'blue')
% legend('u1 real','u1 estimate','Estimation point','FontSize',18)
legend('u2 real','u2 estimate','Estimation point','FontSize',18)
% legend('u real','u bi-level','u basic','prediction point','FontSize',18)
title('u estimation','fontsize',20)
%%
% difference
figure(); hold on
plot( u_cl(1,start:last)-u_pred_basic(1,start:last), 'sr-','LineWidth',1)
plot( u_cl(2,start:last)-u_pred_basic(2,start:last), 'sb-','LineWidth',1)
legend('Error u1','Error u2','FontSize',18)
title('Estimation error $u_{closed \ loop} - u_{estimation}$','interpreter', 'latex','fontsize',20)
%%
start = T+1-5; last = size(u_pred_basic,2);
u2_p1 = u_pred_basic(:,start:last); % u2: D=0; p1: LMI u
u2_r1 = u_cl(:,start:last); % r1: real u
save('./data_for_latex/u2_p1.mat','u2_p1','u2_r1');
%% compute hankel matrix
function Hankel = compute_Hankel_matrix(x,nx,Hankel_row, Hankel_col)
    % Build Hankel matrix by x vector sequences
    if Hankel_row+Hankel_col-1 ~= size(x,2)
        error("The length of history data is wrong!")
    end
    Hankel = zeros(nx*Hankel_row, Hankel_col);
    for i = 1:Hankel_row
        Hankel((i-1)*nx+1:i*nx,:) = x(:,i:i+Hankel_col-1);
    end
end