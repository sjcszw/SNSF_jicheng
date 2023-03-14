% polydome_hankel_prediction_mystruct
clc
clear
rng(123)
%%
% system definition
A = 0.9;
B = 1;%[0.5; 1]; 
C = 1;
Q = 1;R = 1;
nu = size(B, 2);
nx = size(A, 1);
ny = size(C, 1);
fx = @(x,u) A*x + B*u; % state dynamics
fy = @(x) C*x; % output dynamics
%%
N_ini = 3;
N_pred = 10;
T = 50; % T >= (nu + 1)*(N_ini + N_pred + nx) - 1
T_tot = 100; % total steps for simulation

% Simulate a trajectory by random input
x_cl = zeros(nx,T_tot); y_cl = zeros(ny,T_tot);
u_cl = 3*rand(nu,T_tot)-1.5;
x_cl(:, 1) = 0;
y_cl(:, 1) = C*x_cl(:, 1);
for t = 1:T_tot-1
    x_cl(:, t+1) = fx(x_cl(:,t), u_cl(:,t));
    y_cl(:, t+1) = fy(x_cl(:,t+1));
end

%%

h=figure;
hold on
% set(h,'Units','normalized','Position',[0 0 1 .5]); 
yyaxis left
plot(u_cl,'b','LineWidth',1)
yyaxis right
% ylim([5 30])
%     plot(tExp,exp.setpoint_winter.value/10.0,'r','LineWidth',1)
plot(x_cl(1,:),'r','LineWidth',1)
plot(y_cl,'g','LineWidth',1)
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

%% prediction input
% y_pred = Hg_pred*[y_ini;u_ini;u_pred;w_ini;w_pred];
Hg_pred = Hg_problem_u(Hankel_U,Hankel_Y,nu,ny,N_ini,Qg);
%
num = size(u_cl,2);
t = T+1;
while(t<=num-1)
    if t+N_pred<=num
        ui = reshape(u_cl(:,t-N_ini:t-1), [N_ini*nu,1]);
%         wip = reshape(w_cl(:,t-N_ini:t+N_pred-1), [(N_ini+N_pred)*nw,1]);
        yi = reshape(y_cl(:,t-N_ini:t), [(N_ini+1)*ny,1]); 
        yp = reshape(y_cl(:,t+1:t+N_pred), [N_pred*ny,1]); 
        u_pred(:,t:t+N_pred-1) = reshape(Hg_pred*[yi;ui;yp],[nu, N_pred]);
%         y_pred(t+1:t+N_pred) = Hg_pred*[yi;uip;wip];

        g = [Hankel_U_past;Hankel_Y]\[ui;yi;yp];
    %     [Hankel_U;Hankel_y_past]*g-[u_ini_pred;y_ini] % Almost 0
        u_pred_basic(:,t:t+N_pred-1) = reshape(Hankel_U_future*g,[nu, N_pred]);
        t = t+N_pred;
    else
        break;
    end
end
%%
figure()
hold on; grid on;
start = T+1; last = size(u_pred,2);
plot( u_cl(:,start:last), 'Color',[0 0.4470 0.7410],'LineWidth',1)
plot( u_pred(:,start:last), 'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
plot( u_pred_basic(:,start:last), 'Color',[0.4660 0.6740 0.188],'LineWidth',1)
scatter( [(start:N_pred:last)-start+1],u_pred(:,start:N_pred:last), 100,'red')
legend('u real','u bi-level','u basic','prediction point','FontSize',18)

title('u','fontsize',20, 'interpreter', 'latex')
%% prediction output
% % y_pred = Hg_pred*[y_ini;u_ini;u_pred;w_ini;w_pred];
% Hg_pred = Hg_problem_y(Hankel_U,Hankel_Y,ny,N_ini,Qg);
% %
% num = size(u_cl,2);
% t = T+1;
% while(t<=num-1)
%     if t+N_pred<=num
%         uip = reshape(u_cl(:,t-N_ini:t+N_pred-1), [(N_ini+N_pred)*nu,1]);
% %         wip = reshape(w_cl(:,t-N_ini:t+N_pred-1), [(N_ini+N_pred)*nw,1]);
%         yi = reshape(y_cl(:,t-N_ini:t), [(N_ini+1)*ny,1]); 
%         y_pred(:,t+1:t+N_pred) = reshape(Hg_pred*[yi;uip],[ny, N_pred]);
% %         y_pred(t+1:t+N_pred) = Hg_pred*[yi;uip;wip];
% 
%         g = [Hankel_U;Hankel_Y_past]\[uip;yi];
%     %     [Hankel_U;Hankel_y_past]*g-[u_ini_pred;y_ini] % Almost 0
%         y_pred_basic(:,t+1:t+N_pred) = reshape(Hankel_Y_future*g,[ny, N_pred]);
%         t = t+N_pred;
%     else
%         break;
%     end
% end
%% plot
% figure()
% hold on; grid on;
% start = T+1; last = size(y_pred,2);
% plot( y_cl(:,start:last), 'Color',[0 0.4470 0.7410],'LineWidth',1)
% plot( y_pred(:,start:last), 'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
% plot( y_pred_basic(:,start:last), 'Color',[0.4660 0.6740 0.188],'LineWidth',1)
% scatter( [(start:N_pred:last)-start+1],y_pred(:,start:N_pred:last), 100,'red')
% legend('y real','y bi-level','y basic','prediction point','FontSize',18)
% 
% xlabel('Step: 1 step = 15 minutes', 'interpreter', 'latex','fontsize',20);
% ylabel('Temperature [$^{\circ}$C] ', 'interpreter', 'latex','fontsize',20);
% title('Prediction by bi-level deepc: $N_{pred}$ = 12 = 3 hours, T = 5 days','fontsize',20, 'interpreter', 'latex')
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
%% compute the prediction matrix
function Hg_pred = Hg_problem(Hankel_U,Hankel_W,Hankel_Y,ny,N_ini,Qg)
    % y_pred = Hg_pred*[y_ini;u_ini;u_pred;w_ini;w_pred];
            Hy_init = Hankel_Y(1:ny*(N_ini+1), :);
            Hy_pred = Hankel_Y(ny*(N_ini+1)+1:end, :);           
            H = [Hankel_U;Hankel_W];
            
            % ====== Compute g_norm from KKT condition
            temp = inv(Qg);
            % desired inversion
            temp = temp - temp*Hy_init'*((eye((N_ini+1)*ny)... %/(2*w_s)...
                +Hy_init*temp*Hy_init')\Hy_init)*temp;
            temp_inv = temp*H'/(H*temp*H');
         
            Hg_pred = Hy_pred*[(temp-temp*H'*((H*temp*H')\H)*temp)*Hy_init', temp_inv];
end

function Hg_pred = Hg_problem_y(Hankel_U,Hankel_Y,ny,N_ini,Qg)
    % u_pred = Hg_pred*[y_ini;u_ini;u_pred;w_ini;w_pred];
            Hy_init = Hankel_Y(1:ny*(N_ini+1), :);
            Hy_pred = Hankel_Y(ny*(N_ini+1)+1:end, :);           
            H = Hankel_U;
            
            % ====== Compute g_norm from KKT condition
            temp = inv(Qg);
            % desired inversion
            temp = temp - temp*Hy_init'*((eye((N_ini+1)*ny)... %/(2*w_s)...
                +Hy_init*temp*Hy_init')\Hy_init)*temp;
            temp_inv = temp*H'/(H*temp*H');
            Hg_pred = Hy_pred*[(temp-temp*H'*((H*temp*H')\H)*temp)*Hy_init', temp_inv];
end

function Hg_pred = Hg_problem_u(Hankel_U,Hankel_Y,nu,ny,N_ini,Qg)
    % u_pred = Hg_pred*[y_ini;u_ini;y_pred];
            Hy_init = Hankel_Y(1:ny*(N_ini+1), :);
            Hy_pred = Hankel_Y(ny*(N_ini+1)+1:end, :); 
            Hu_init = Hankel_U(1:nu*N_ini, :);
            Hu_pred = Hankel_U(nu*N_ini+1:end, :);           
            H = [Hu_init;Hy_pred];
            
            % ====== Compute g_norm from KKT condition
            temp = inv(Qg);
            % desired inversion
            temp = temp - temp*Hy_init'*((eye((N_ini+1)*ny)... %/(2*w_s)...
                +Hy_init*temp*Hy_init')\Hy_init)*temp;
            temp_inv = temp*H'/(H*temp*H');
            Hg_pred = Hu_pred*[(temp-temp*H'*((H*temp*H')\H)*temp)*Hy_init', temp_inv];
end