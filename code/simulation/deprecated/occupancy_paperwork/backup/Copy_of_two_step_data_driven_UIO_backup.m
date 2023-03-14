% This version: when start to build it, just finished the first Kalman
% filter
%% 16052022
% 1. System id: check model with training and validation data
% 2. Use model estimate the x (in training set),  
% 3. [Yp Xp Yf Xf] estimate x then compute u by x and model 
%% ====== 1. System id ======
clear
%% --- 1.1 Load and prepare the data
T_train_end = 96*4+50;
T_test_start = 96*7+1;
T_tot = 1332; % size(exp.co2_1, 1); % total steps for simulation

ifDetrend = 1;
ifTemp = 0; % using temperature related information
ifExp = 1; % 0: simulation 
ifU = 0; % In simulation, whether use known u
ifIno = 0; % If the state space is in the inovation form

if ifExp
    load('raw_06-12-2021_19-12-2021.mat')
    if ifTemp
        y_cl = [mean([exp.sensor_temp_1,exp.sensor_temp_2,exp.sensor_temp_3,exp.sensor_temp_4],2)'; ...
            mean([exp.co2_1,exp.co2_2,exp.co2_3,exp.co2_4],2)'/100];
    else
        y_cl = mean([exp.co2_1,exp.co2_2,exp.co2_3,exp.co2_4],2)'/100;
    end
    
    u_cl = exp.people';
    w_cl = [exp.power-2.35,exp.weather_temp,exp.weather_rad/100]';

    % detrend
    if ifDetrend==1
        ny = size(y_cl,1);
        y_cl_env = mean([exp.co2_1,exp.co2_2,exp.co2_3,exp.co2_4]')/100;
        y_cl_env = [y_cl_env(:,96*7+1:96*7+32), y_cl_env(:,96*6+33:96*7)];    
        % wrong: 
%              y_cl_env = kron(ones(ny,14),y_cl_env); 
%              y_cl_env = y_cl_env(1:size(y_cl,2));
         % ======== 16052022
            y_cl_env = kron(ones(1,14),y_cl_env); 
            if ifTemp
                % detrend the co2
                y_cl_env = [zeros(1,size(y_cl,2)); y_cl_env(1,1:size(y_cl,2))];
                % detrend all
%                 y_cl_env = [mean(y_cl(1,:))*ones(1,size(y_cl,2)); y_cl_env(1,1:size(y_cl,2))];
%                 w_cl = w_cl - mean(w_cl,2);
            else
                y_cl_env = y_cl_env(1,1:size(y_cl,2));
            end
        % ======== 16052022
        
        y_cl = y_cl - y_cl_env;
    end


else
    rng(128)
    A = 2*rand(3,3); eig(A)
    B = 2*rand(3,2);
    C = 2*rand(2,3);
    D = zeros(2,2);
    Ob = obsv(A,C);Co = ctrb(A,B);
    % rank(Ob) 
    % rank(Co)
    if ifU
        Bw = 2*rand(3,2); % for known u
        nw = size(Bw, 2);
        w_cl = 2*rand(nw,T_tot)-1;
    else
        Bw = zeros(3,2);
        nw = size(Bw, 2);
        w_cl = zeros(nw,T_tot)-1;        
    end
    nu = size(B, 2);
    nx = size(A, 1);
    ny = size(C, 1);
    fx = @(x,u, w) A*x + B*u + Bw*w; % state dynamics
    fy = @(x) C*x; % output dynamics
    a = tzero(A,B,C,D,eye(nx))    
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
        x_cl(:, t+1) = fx(x_cl(:,t), u_cl(:,t), w_cl(:,t));
        y_cl(:, t+1) = fy(x_cl(:,t+1));
    end    
end
%% iden and test data structure for ssest (PEM)
T_samp = 900;

if ifTemp
    start = 1; last = T_train_end;
    exp_iden = iddata(	y_cl(:,start:last)',		...
                [u_cl(:,start:last); w_cl(:,start:last)]', T_samp); %,
    start = T_test_start; last = T_tot;
    exp_test = iddata(	y_cl(:,start:last)',		...
                [u_cl(:,start:last); w_cl(:,start:last)]', T_samp); %,    
else
    start = 1; last = T_train_end;
    exp_iden = iddata(	y_cl(:,start:last)',		...
                [u_cl(:,start:last)]', T_samp); %,
    start = T_test_start; last = T_tot;
    exp_test = iddata(	y_cl(:,start:last)',		...
                [u_cl(:,start:last)]', T_samp); %,
end
%     exp_elec.InputName	= {'electrical power', 'air_temp', 'solar_rad','people'};
%     exp_elec.OutputName = {'sensor_temp'};
%% --- 1.2 model selection
% nn = struc(1:7, 1:7,1);
% V = arxstruc(exp_iden,exp_test,nn);
% selstruc(V); %invoke selstruc in an interactive mode
%
% nx = 1:10;
% sys = ssest(exp_iden,nx);

%% 1.3 --- Comparation of different order of ss model ---
opt = ssestOptions('SearchMethod','auto');
opt.Focus = 'prediction'; 
opt.N4Horizon = 1*3600/T_samp;

nx = 2;
[m1, X0] = ssest(exp_iden, nx, 'Ts', T_samp, opt);
nx = 3;
[m2, X0] = ssest(exp_iden, nx, 'Ts', T_samp, opt);
nx = 1;
[m3, X0] = ssest(exp_iden, nx, 'Ts', T_samp, opt);

pred_step = 4;
[~,fit1,~] = compare(exp_iden,m1,m2,m3,pred_step)
[~,fit,~] = compare(exp_test,m1,m2,m3,pred_step)
if ifTemp || ~ifExp
fit1{1}
fit1{2}
fit1{3}
fit{1}
fit{2}
fit{3}
end
%%
% num = 1;
% figure
% compare(exp_test(:,:,:,num), m1,m2,m3, 1);
% figure
% compare(exp_test(:,:,:,num), m1,m2,m3,4);
% figure
% compare(exp_test(1:96,:,:,num), m1,m2,m3);

%
%% ====== 3. Design UIO ======
% 1. Load parameters
if ifExp
    model = m2;
    A = model.A
    B = model.B
    C = model.C
    D = model.D
    ny = size(C, 1);
    nx = size(A, 1);
    nu = 1;
    nw = size(B, 2) - nu;
    Bw = B(:, nu+1:end); 
    B = B(:,1:nu);    
    D = D(:,1:nu); 
else
    1
end
%
%% Luenberger observer and Kalman filter for state estimation
if ifExp
    % Kalman filter
    W = model.K*model.K'*model.NoiseVariance; 
    V = model.NoiseVariance;
    M = model.K;
    x_est(:,1) = zeros(nx,1);
    P_post = W;    
    for i = 1:T_tot-1
 		% Time Update
        x_pri	= A * x_est(:,i) + B * u_cl(:, i);
        P_pri	= A * P_post * A' + W;        
        if ifIno
            % Measurements Update
            k_post	= (P_pri * C' + M) * (C * P_pri * C' + V + C*M + M'*C' )^-1;
            x_est(:, i+1)	= x_pri + k_post * (y_cl(:,i) - C * x_pri);
            P_post	= (eye(nx) - k_post*C) * P_pri - k_post*M';           
        else
            % Measurements Update
            k_post	= P_pri * C'* (C * P_pri * C' + V)^-1;
            x_est(:, i+1)	= x_pri + k_post * (y_cl(:,i) - C * x_pri);
            P_post	= (eye(nx) - k_post*C) * P_pri;
        end
%         L_kal = A*qvar_kal*C'*(C*qvar_kal*C' + V)^(-1);
%         x_est(:,i+1) = A*x_est(:,i) + B*u_cl(:,i) + L_kal*(y_cl(:,i)-C*x_est(:,i));
%         qvar_kal = qvar_kal - qvar_kal*C'*(C*qvar_kal*C' + V)^(-1)*C*qvar_kal;
%         qvar_kal = A*qvar_kal*A' + W;
    end   
    figure()
    hold on
    plot(C*x_est(:,50:T_tot) - y_cl(:,50:T_tot), 'r')
    mean(abs(C*x_est(:,50:T_tot) - y_cl(:,50:T_tot)))

else
    Q = eye(nx); 
    R = 1;
    [K,~,~] = dlqr(A',C',Q,R);
    K = -K;                             
    eig(A'+ C'*K)
    % Luenberger observer
    clear x_elec_est y_pred_elec
    x_est(:,1) = zeros(nx,1);

    for i = 1:T_tot-1
        x_est(:,i+1) = A*x_est(:,i) + B*u_cl(:,i) - K'*(y_cl(:,i)-C*x_est(:,i));
    end
    figure()
    hold on
    plot(C*x_est(:,50:T_tot) - y_cl(:,50:T_tot), 'r')
    mean(abs(C*x_est(:,50:T_tot) - y_cl(:,50:T_tot)))
    
end

%% --- 2.2 Check conditions
rng(123)
% K = eye(ny);
K = rand(ny,ny);
fprintf('====== Condition 1 \n')
rank([B;D]) == nu
rank([K*C*B; D]) == nu


fprintf('====== U nonsingular and inverse \n')
% if ifExp==1
%     U = [zeros(2*ny-nu,nu); eye(nu)] / [K*C*B; D]
%     U(1,2)  = 10*rand(1); U(2,2) = 10*rand(1);
%     det(U)
%     U*[K*C*B; D] - [zeros(2*ny-nu,nu); eye(nu)]
% else
    U = [zeros(2*ny-nu,nu); eye(nu)] / [K*C*B; D]
    tem = null([K*C*B; D]'); j =1;
    for i =1:size(U,1)
        if sum(U(i,:))==0
            U(i,:) = tem(:,j)';
            j = j+1;
        end
    end
    det(U)
    U*[K*C*B; D] - [zeros(2*ny-nu,nu); eye(nu)]   
% end
%%
fprintf('====== Condition 2 \n')
tem = -U* [K*C*A; C];
T1 = tem(1:2*ny-nu,:);
T2 = tem(2*ny-nu+1:end,:);
Ob = obsv(A+B*T2, T1);
rank(Ob) - nx
fprintf('====== Condition 2 by zeros \n')
a = tzero(A,B,C,zeros(ny,nu),eye(nx))

% Compute one X1
eig(A+B*T2)
Q = eye(nx); R = eye(2*ny-nu);
[K1,~,~] = dlqr((A+B*T2)',T1',Q,R);
X1 = -K1';
% X1 = -0*K1';
eig(A+B*T2 + X1*T1)

% Compute N E L
N = A+B*T2 + X1*T1;
tem = [X1 B]*U;
ME = tem(:,1:ny);
Y = tem(:,ny+1:end);
E = ME*K; 
L = ((eye(nx)-ME*K*C)*A - Y*C)*ME*K + Y;
%% --- 2.3 estimate input
u_pred = zeros(nu,T_tot);
z_pred = zeros(nx,T_tot);
x_pred = zeros(nx,T_tot);
x_pred(:,1) = z_pred(:,1) + E*y_cl(:,1); 
for t = 2:T_tot
        if ifTemp
        else
            if ifU || ifTemp
                z_pred(:,t) = N*z_pred(:,t-1) + L*y_cl(:,t-1) + (eye(nx)-E*C)*Bw*w_cl(:,t-1);
            else
                z_pred(:,t) = N*z_pred(:,t-1) + L*y_cl(:,t-1);
            end
            x_pred(:,t) = z_pred(:,t) + E*y_cl(:,t); 
        end
            
end

for t = 2:T_tot
        if ifTemp || ifU
            u_pred(:,t-1) = B\(x_pred(:,t)-A*x_pred(:,t-1) - Bw*w_cl(:,t-1));
        else
            u_pred(:,t-1) = B\(x_pred(:,t)-A*x_pred(:,t-1)); %inv(B)*
%             x_pred(:,t)-A*x_pred(:,t-1) - B*u_pred(t-1)
        end
end
%%
num = 1;
figure()
hold on; grid on;
% start = 97; last = T_train_end;
start = 7*96+1; last = 12*96;
% last = size(u_pred,2)-96*2;
% plot( u_cl(num,start:last), 'Color',[0 0.4470 0.7410],'LineWidth',1,'Marker','s')
% plot( u_pred(num,start:last), 'Color',[0.8500 0.3250 0.0980],'LineWidth',1,'Marker','o')
plot(u_cl(num,start:last),'b--');
plot(u_pred(num,start:last),'r','LineWidth',1.5);
% plot(u_pred_basic(num,start:last),'g','LineWidth',1.5);
% scatter( [(start:N_pred:last)-start+1],u_pred(num,start:N_pred:last), 100,'red')
legend('Real people','Estimated people','FontSize',18) % 'Estimation point',

title('Estimate occpancy number in Polydome','fontsize',20, 'interpreter', 'latex')
mean(abs(u_pred(num,start:last)-u_cl(num,start:last)))
% mean(abs(u_pred_basic(num,start:last)-u_cl(num,start:last)))
% CO2 m2
%   No de: 12.2048 0night, 9.840
%   De, right, 9.936 0night, 8.39
% CO2+temp
%   No de: m2 8.4665 0night, 10.879
%   De: very bad
1
%%
u_pred = zeros(nu,T_tot);
z_pred = zeros(nx,T_tot);
x_pred = zeros(nx,T_tot);
x_pred(:,1) = z_pred(:,1) + E*y_cl(:,1); 
for i = 0:14
    z_pred(:,i*96+1)=0;
    for t = i*96+1:(i+1)*96
        if t+2>T_tot
            break;
        end
        if t<=i*96+27 || t>=(i+1)*96-4
            x_pred(:,t+1) = A*x_pred(:,t);
            z_pred(:,t+1) = x_pred(:,t+1) - E*y_cl(:,t+1);
        else
            z_pred(:,t+1) = N*z_pred(:,t) + L*y_cl(:,t);
            x_pred(:,t+1) = z_pred(:,t+1) + E*y_cl(:,t+1);
        end
           
    end
end
for t = 2:T_tot
        if ifTemp
            u_pred(:,t-1) = B\(x_pred(:,t)-A*x_pred(:,t-1) - Bw*w_cl(:,t-1));
        else
            u_pred(t-1) = B\(x_pred(:,t)-A*x_pred(:,t-1)); % inv(B)*
%             x_pred(:,t)-A*x_pred(:,t-1) - B*u_pred(t-1)
        end
end
mean(abs(u_pred(num,start:last)-u_cl(num,start:last)))
% %%
% u_pred = zeros(nu,T_tot);
% z_pred = zeros(nx,T_tot);
% x_pred = zeros(nx,T_tot);
% x_pred(:,1) = z_pred(:,1) + E*y_cl(:,1); 
% for i = 0:14
%     z_pred(:,i*96+1)=0;
%     for t = i*96+1:(i+1)*96
%         if t+1<=T_tot
%             z_pred(:,t+1) = N*z_pred(:,t) + L*y_cl(:,t);
%             x_pred(:,t+1) = z_pred(:,t+1) + E*y_cl(:,t+1);            
%         else
%            
%             break;
%         end
%     end
% end
%