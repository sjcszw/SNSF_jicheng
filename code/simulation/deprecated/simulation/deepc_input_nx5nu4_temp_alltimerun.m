% polydome_hankel_prediction_mystruct
clc
clear
rng(125)
%%
% system definition
% A = [1 1; 0 1];
% B = [1; 1];%[0.5; 1]; 
% A = diag([0.9,0.8]);
% B = eye(2);
% C = eye(2);
% rank(Ob) 
% rank(Co)
nu = 1;
nw = 3;
ny = 2;
%%
load('raw_06-12-2021_19-12-2021.mat')
%
N_ini = 5;
N_pred = 1;

% y_cl = mean([exp.sensor_temp_1,exp.sensor_temp_2,exp.sensor_temp_3,exp.sensor_temp_4],2)';
y_cl = [mean([exp.sensor_temp_1,exp.sensor_temp_2,exp.sensor_temp_3,exp.sensor_temp_4],2)';
        mean([exp.co2_1,exp.co2_2,exp.co2_3,exp.co2_4],2)'/100];
% y_cl = [ mean([exp.co2_1,exp.co2_2,exp.co2_3,exp.co2_4],2)'/100];

% No calibration at night; co2 and temp improve it
% co2, N_ini=10 raw 6.9910, clear 25.44 N_ini=4 raw = 6.3313
% temp: N_ini=10 raw 14.1492 clear:55.3524 N_ini=4 raw 11.4557
% temp + co2: N_ini=10 raw 6.8423, clear 17.257 N_ini=4 raw 5.3618, clear 26.4552
%       N_ini=5 raw 5.5563 N_ini=6 raw 7.0524 N_ini=7 raw 6.8705
% detrend 
    y_cl_env = mean([exp.co2_1,exp.co2_2,exp.co2_3,exp.co2_4]')/100;
    y_cl_env = [y_cl_env(:,96*7+1:96*7+32), y_cl_env(:,96*6+33:96*7)];
    y_cl_env = kron(ones(ny,14),y_cl_env);
    y_cl_env = y_cl_env(1:size(y_cl,2));
    y_cl = y_cl - y_cl_env;

u_cl = exp.people';
w_cl = [exp.power-2.35,exp.weather_temp,exp.weather_rad/100]';

T = 96*4+50; % T >= (nu + 1)*(N_ini + N_pred + nx) - 1
T_tot = size(y_cl,2); % total steps for simulation
%%

h=figure;
hold on
% set(h,'Units','normalized','Position',[0 0 1 .5]); 
yyaxis left
plot(u_cl(1,:),'b','LineWidth',1)
yyaxis right
% ylim([5 30])
%     plot(tExp,exp.setpoint_winter.value/10.0,'r','LineWidth',1)
% plot(C(1,:)*x_cl,'r','LineWidth',1)
plot(y_cl(1,:),'g','LineWidth',1)
% plot(tem1(1,:),'r','LineWidth',1)
% plot(time,exp.supply_temp/10.0,'y','LineWidth',1.5)
% title('Meteorological data','FontSize',18)
legend({'u','y'},'FontSize',18)
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
Hankel_W = compute_Hankel_matrix(w_cl(:,1:T),nw,N_total, Hankel_col);
% Hankel_W_elec = compute_Hankel_matrix(exp_elec_iden.w,nw,N_total, Hankel_col);

N_total = N_ini + N_pred + 1;
Hankel_col = T - N_total + 2;
Hankel_Y = compute_Hankel_matrix(y_cl(:,1:T+1),ny,N_total, Hankel_col);
Hankel_U_past = Hankel_U(1:nu*N_ini, :);
Hankel_U_future = Hankel_U(nu*N_ini+1:end, :);
Hankel_Y_past = Hankel_Y(1:ny*(N_ini+1), :);
Hankel_Y_future = Hankel_Y(ny*(N_ini+1)+1:end, :);
%%
nx = 5;
N_total = N_ini + N_pred + nx;
Hankel_col = T - N_total +1 ;

% Need: u is persistently exciting of order is N_ini + N_pred + nx
Hankel_UW_check = compute_Hankel_matrix([u_cl(:,1:T);w_cl(:,1:T)], nu+nw, N_total, Hankel_col);
% Hankel_UW_check = compute_Hankel_matrix([exp_elec_iden.u; exp_elec_iden.w], nu+nw, N_total, Hankel_col);

% Check if the past u trajectory is persistently exciting of order N_ini + N_pred + nx
if rank(Hankel_UW_check)== (nu+nw)*(N_ini + N_pred + nx) 
    disp('Hankel rank is ok')
else
    disp('Exciting of order of Hu is samller than N_ini + N_pred + nx')
end

%% prediction input
% y_pred = Hg_pred*[y_ini;u_ini;u_pred;w_ini;w_pred];
Hg_pred = Hg_problem_u(Hankel_U,Hankel_Y,Hankel_W, nu,ny,N_ini,Qg);
%
num = size(u_cl,2);
t = T+1;
u_pred = zeros(nu,T_tot);


for t = 7*96+7*4+1:T_tot
    if t+N_pred<=T_tot
        ui = reshape(u_pred(:,t-N_ini:t-1), [N_ini*nu,1]);
        wip = reshape(w_cl(:,t-N_ini:t+N_pred-1), [(N_ini+N_pred)*nw,1]);
        yi = reshape(y_cl(:,t-N_ini:t), [(N_ini+1)*ny,1]); 
        yp = reshape(y_cl(:,t+1:t+N_pred), [N_pred*ny,1]); 
        u_pred(:,t:t+N_pred-1) = reshape(Hg_pred*[yi;ui;yp;wip],[nu, N_pred]);
%         y_pred(t+1:t+N_pred) = Hg_pred*[yi;uip;wip];

%             g = [Hankel_U_past;Hankel_Y;Hankel_W]\[ui;yi;yp;wip];
    %     [Hankel_U;Hankel_y_past]*g-[u_ini_pred;y_ini] % Almost 0
%             u_pred_basic(:,t:t+N_pred-1) = reshape(Hankel_U_future*g,[nu, N_pred]);
%             t = t+N_pred;

    else
        break;
    end
end
%%
num = 1;
figure()
hold on; grid on;
start = 7*96+1; last = size(u_pred,2)-96*2;
% plot( u_cl(num,start:last), 'Color',[0 0.4470 0.7410],'LineWidth',1,'Marker','s')
% plot( u_pred(num,start:last), 'Color',[0.8500 0.3250 0.0980],'LineWidth',1,'Marker','o')
plot(u_cl(num,start:last),'b--');
plot(u_pred(num,start:last),'r','LineWidth',1.5);
% scatter( [(start:N_pred:last)-start+1],u_pred(num,start:N_pred:last), 100,'red')
legend('Real people','Estimated people','FontSize',18) % 'Estimation point',

title('Estimate occpancy number in Polydome','fontsize',20, 'interpreter', 'latex')
mean(abs(u_pred(num,start:last)-u_cl(num,start:last)))
%% The training dataset

for t = N_ini+1:T
    if t+N_pred<=T
        ui = reshape(u_pred(:,t-N_ini:t-1), [N_ini*nu,1]);
        wip = reshape(w_cl(:,t-N_ini:t+N_pred-1), [(N_ini+N_pred)*nw,1]);
        yi = reshape(y_cl(:,t-N_ini:t), [(N_ini+1)*ny,1]); 
        yp = reshape(y_cl(:,t+1:t+N_pred), [N_pred*ny,1]); 
        u_pred(:,t:t+N_pred-1) = reshape(Hg_pred*[yi;ui;yp;wip],[nu, N_pred]);
    else
        break;
    end
end
figure();hold on;
start = N_ini+1; last = T;
plot(u_cl(num,start:last),'b--');
plot(u_pred(num,start:last),'r','LineWidth',1.5);
xlabel('time');
ylabel('y');
legend('Data','Predictions');
hold off
%%
tzero(A,B,C,zeros(ny,nu),eye(nx))
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

function Hg_pred = Hg_problem_u(Hankel_U,Hankel_Y,Hankel_W,nu,ny,N_ini,Qg)
    % u_pred = Hg_pred*[y_ini;u_ini;y_pred];
            Hy_init = Hankel_Y(1:ny*(N_ini+1), :);
            Hy_pred = Hankel_Y(ny*(N_ini+1)+1:end, :); 
            Hu_init = Hankel_U(1:nu*N_ini, :);
            Hu_pred = Hankel_U(nu*N_ini+1:end, :);           
            H = [Hu_init;Hy_pred;Hankel_W];
            
            % ====== Compute g_norm from KKT condition
            temp = inv(Qg);
            % desired inversion
            temp = temp - temp*Hy_init'*((eye((N_ini+1)*ny)... %/(2*w_s)...
                +Hy_init*temp*Hy_init')\Hy_init)*temp;
            temp_inv = temp*H'/(H*temp*H');
            Hg_pred = Hu_pred*[(temp-temp*H'*((H*temp*H')\H)*temp)*Hy_init', temp_inv];
end