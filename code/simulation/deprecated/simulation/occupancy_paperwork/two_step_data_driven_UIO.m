%% 16052022
% 1. System id: check model with training and validation data
% 2. Use model estimate the x (in training set),  
% 3. [Yp Xp Yf Xf] estimate x then compute u by x and model
%% Save u1 D~=0 u2 D=0 p4 data-driven UIO
% start = T_test_start-5; last = T_tot-1;
% u2_p4 = u_pred(:,start:last); 
% u2_r4 = u_cl(:,start:last); % r1: real u
% save('./u2_p4_1805.mat','u2_p4','u2_r4');
% start = T_test_start-5; last = T_tot;
% u1_p4 = u_pred(:,start:last); 
% u1_r4 = u_cl(:,start:last); % r1: real u
% save('./u1_p4_1805.mat','u1_p4','u1_r4');
% start = 7*96+1; last = 12*96;
% ur_DUIO = u_cl(num,start:last);
% u_DUIO = u_pred(num,start:last);
% save('./DUIO.mat','ur_DUIO','u_DUIO')
%% ====== 1. System id ======
clear
%% --- 1.1 Load and prepare the data
T_train_end = 96*4+50;
T_test_start = 96*7+1;
T_tot = 1332; % size(exp.co2_1, 1); % total steps for simulation

ifDetrend = 1;
ifTemp = 0; % using temperature related information
ifExp = 1; % 0: simulation 
ifIno = 1; % If the state space is in the inovation form
ifU = 0; % In simulation, whether use known u
ifD = 0;

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
%                 y_cl_env = [zeros(1,size(y_cl,2)); y_cl_env(1,1:size(y_cl,2))];
                % detrend all
                y_cl_env = [mean(y_cl(1,:))*ones(1,size(y_cl,2)); y_cl_env(1,1:size(y_cl,2))];
                w_cl = w_cl - mean(w_cl,2);
            else
                y_cl_env = y_cl_env(1,1:size(y_cl,2));
            end
        % ======== 16052022
        
        y_cl = y_cl - y_cl_env;
    end


else
    rng(353)
    T_train_end = 50;
    T_test_start = 51;
    T_tot = 100;  
    n1 = 20; n2 = 0;
    A = round(n1*rand(3,3)-n2)/10; eig(A)
    B = round(n1*rand(3,2)-n2)/10; %
    C = round(n1*rand(2,3)-n2)/10;
    nu = size(B, 2);
    nx = size(A, 1);
    ny = size(C, 1);    
    if ifD
        for i=1:1000
            D =  round(n1*rand(2,2)-n2)/10; %[1 0; 1 0];%zeros(ny,nu); 
            a = tzero(A,B,C,D,eye(nx));
            if max(abs(a))<0.8
                abs(a)
                break
            end
        end
    else
        D =  0*round(n1*rand(2,2)-n2)/10;    
    end
    a = tzero(A,B,C,D,eye(nx))
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
        w_cl = zeros(nw,T_tot);        
    end
    fx = @(x,u, w) A*x + B*u + Bw*w; % state dynamics
    fy = @(x,u) C*x + D*u; % output dynamics
    
    % Simulate a trajectory by random input
    x_cl = zeros(nx,T_tot); y_cl = zeros(ny,T_tot);
    Q = eye(nx);R = eye(nu);
    [K,~,~] = dlqr(A,B,Q,R);
    K = -K;
    u_cl = 3*rand(nu,T_tot)-1.5;
    x_cl(:, 1) = 0;
    u_cl(:,1) = u_cl(:,1) + K*x_cl(:,1);
    y_cl(:, 1) =  fy(x_cl(:,1),u_cl(:,1));
    for t = 1:T_tot-1
        x_cl(:, t+1) = fx(x_cl(:,t), u_cl(:,t), w_cl(:,t));
        u_cl(:,t+1) = u_cl(:,t+1) + K*x_cl(:,t+1);
        y_cl(:, t+1) = fy(x_cl(:,t+1),u_cl(:,t+1));
    end    
end
%% iden and test data structure for ssest (PEM)
T_samp = 900;

if ifTemp || ifU
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

if ifExp
    pred_step = 4;
    opt.Focus = 'prediction'; 
    opt.N4Horizon = 1*3600/T_samp;
     order = 7;
    [m1, X0] = ssest(exp_iden, order, 'Ts', T_samp, opt);
    order = 3;
    [m2, X0] = ssest(exp_iden, order, 'Ts', T_samp, opt);
    order = 2;
    [m3, X0] = ssest(exp_iden, order, 'Ts', T_samp, opt);

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
else
    pred_step = 1;
    order = nx;
    [model, X0] = ssest(exp_iden, order, 'Ts', T_samp, opt, 'Feedthrough',ifD);

    [~,fit1,~] = compare(exp_iden,model,pred_step)
    [~,fit,~] = compare(exp_test,model,pred_step)   
%     opt.N4Horizon = 1;
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
%% ====== 2.  State estimation ======
% 1. Load parameters
if ifExp
    model = m1;
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
    A = model.A
    B = model.B
    C = model.C
    D = model.D    
end
%
%% Luenberger observer and Kalman filter for state estimation
if ifExp
    % Kalman filter
    % ifIno, 1~=3, steady K, 0 slowly updated K
    W = model.K*model.NoiseVariance*model.K'; 
    V = model.NoiseVariance;
    M = model.K;
    x_est(:,1) = zeros(nx,1);
    P_post = W;    
    for i = 1:T_tot
        if ifIno==1 
            % === Luca
            W = eye(nx);
            S = M*V;
            Dsr = S*V^-1;
            Phi = A-Dsr*C;
            Q_mod = W - S*V^-1*S';
            % Time Update
            if ifTemp
                x_pri	= Phi  * x_est(:,i) + B * u_cl(:, i) + Bw*w_cl(:,i) + Dsr * y_cl(:,i);
            else
                x_pri	= Phi  * x_est(:,i) + B * u_cl(:, i) + Dsr * y_cl(:,i);
            end
            P_pri	= Phi * P_post * Phi' + Q_mod;     
             % Measurements Update
            k_post	= P_pri * C'* (C * P_pri * C' + V)^-1;
            x_est(:, i+1)	= x_pri + k_post * (y_cl(:,i) - C * x_pri);
            P_post	= (eye(nx) - k_post*C) * P_pri;    
        elseif  ifIno==2 % use the one from book
            % === ifIno
            % Time Update
            x_pri	= A * x_est(:,i) + B * u_cl(:, i);
            P_pri	= A * P_post * A' + W;     
            % Measurements Update
            k_post	= (P_pri * C' + M) * (C * P_pri * C' + V + C*M + M'*C' )^-1;
            x_est(:, i+1)	= x_pri + k_post * (y_cl(:,i) - C * x_pri);
            P_post	= (eye(nx) - k_post*C) * P_pri - k_post*M';     
        elseif ifIno==3 % use model.K as gain
%             for i = 1:T_tot-1
                if ifTemp
                    x_est(:,i+1) = A*x_est(:,i) + B*u_cl(:,i) + Bw*w_cl(:,i) + M*(y_cl(:,i)-C*x_est(:,i));
                else
                    x_est(:,i+1) = A*x_est(:,i) + B*u_cl(:,i) + M*(y_cl(:,i)-C*x_est(:,i));
                end
%             end
            
        elseif ifIno==0 % basic KF
            % Time Update
            if ifTemp
                x_pri = A * x_est(:,i) + B * u_cl(:, i) + Bw*w_cl(:,i);
            else           
                x_pri	= A * x_est(:,i) + B * u_cl(:, i);
            end
            P_pri	= A * P_post * A' + W;                     
            % Measurements Update
            k_post	= P_pri * C'* (C * P_pri * C' + V)^-1;
            x_est(:, i+1)	= x_pri + k_post * (y_cl(:,i) - C * x_pri);
            P_post	= (eye(nx) - k_post*C) * P_pri;
        end
    end   
    figure()
    hold on
    num = 1;
    plot(C(num,:)*x_est(:,1:T_tot) - y_cl(num,1:T_tot), 'r')
    mean(abs(C*x_est(:,50:T_tot) - y_cl(:,50:T_tot)),2)
% x_est(:,1:5)
else
    Q = eye(nx); 
    R = eye(nu);
    [K,~,~] = dlqr(A',C',Q,R);
    K = -K;                             
    eig(A'+ C'*K)
    % Luenberger observer
    clear x_elec_est y_pred_elec
    x_est(:,1) = ones(nx,1);

    for i = 1:T_tot
        if ifU
            x_est(:,i+1) = A*x_est(:,i) + B*u_cl(:,i) + Bw*w_cl(:,i) - K'*(y_cl(:,i)-C*x_est(:,i));
        else
            if ifD
                x_est(:,i+1) = A*x_est(:,i) + B*u_cl(:,i) - K'*(y_cl(:,i) - C*x_est(:,i) - D*u_cl(:,i));
            else
                x_est(:,i+1) = A*x_est(:,i) + B*u_cl(:,i) - K'*(y_cl(:,i)-C*x_est(:,i));
            end
        end
    end
    figure()
    hold on
    num = 1;
    if ~ifD
        plot(C(num,:)*x_est(:,1:T_tot) - y_cl(num,1:T_tot), 'r')
        mean(abs(C*x_est(:,1:T_tot) - y_cl(:,1:T_tot)),2)
    else
        plot(C(num,:)*x_est(:,1:T_tot)+D(num,:)*u_cl(:,1:T_tot) - y_cl(num,1:T_tot), 'r')
        mean(abs(C*x_est(:,10:T_tot)+D*u_cl(:,10:T_tot)  - y_cl(:,10:T_tot)),2)
    end
end

for t = 2:T_tot+1
        if ifTemp || ifU
            u_pred(:,t-1) = B\(x_est(:,t)-A*x_est(:,t-1) - Bw*w_cl(:,t-1));
        else
            u_pred(:, t-1) = B\(x_est(:,t)-A*x_est(:,t-1)); % inv(B)*
%             x_pred(:,t)-A*x_pred(:,t-1) - B*u_pred(t-1)
        end
end
start = 0*96+10; last = T_tot;
mean(abs(u_pred(1,start:last)-u_cl(1,start:last)))
% for i = start:last
%     u_pred(1,i) = max(u_pred(1,i),0);
%     if u_pred(1,i)>80
%         u_pred(1,i) = mean(u_pred(1,i-1:i));
%     end
% end
% mean(abs(u_pred(1,start:last)-u_cl(num,start:last)))
figure()
hold on; grid on;
plot(u_cl(1,start:last),'b--');
plot(u_pred(1,start:last),'r','LineWidth',1.5);


%% ====== 3. Design UIO ======
% build the Hankel matrix 
T = T_train_end;
T_start = 10;
disp('Computing Hankel...')
N_total = 2;
Hankel_col = T - T_start +1 - N_total + 1;
Hankel_Y = compute_Hankel_matrix(y_cl(:,T_start:T),ny,N_total, Hankel_col);
Hankel_X = compute_Hankel_matrix(x_est(:,T_start:T),nx,N_total, Hankel_col);

Hankel_X_past = Hankel_X(1:nx, :);
Hankel_X_future = Hankel_X(nx+1:end, :);

if ifTemp || ifU
    Hankel_W = compute_Hankel_matrix(w_cl(:,T_start:T),nw,N_total, Hankel_col);
end
%% 
if ifTemp || ifU
    G = pinv([Hankel_X_past; Hankel_Y; Hankel_W]); % Moore-Penrose pseudoinverse
else
    G = pinv([Hankel_X_past;Hankel_Y]); %;Hankel_W]); % Moore-Penrose pseudoinverse
end
X1 = Hankel_X_future*G(:,1:nx);
abs(eig(X1))   

% tem3 = [Hankel_X_past;Hankel_Y];
% tem3*tem*tem3-tem3

% %% SVD first
% [U,S,V] = svd(tem3);
% % G = V*[inv(S(1:15,1:15)) rand(15,9);rand(30,24)]*U';
% num1 = 7;
% G = V*[inv(S(1:num1,1:num1)); zeros(424-num1,num1)]*U';
% % G = V*[inv(S(1:16,1:16)) zeros(16,37-16); zeros(429-16,37)]*U'; G = V*S'*U';
% a = tem3*G*tem3-tem3;
%% prediction input
u_pred = zeros(nu,T_tot);
z_pred = zeros(nx,T_tot);
x_pred = zeros(nx,T_tot);

if ifExp
    tem = T_test_start;
    for t = tem:T_tot-1
                if ifU || ifTemp
                    yip = reshape(y_cl(:,t:t+1), [2*ny,1]); 
                    xi = reshape(x_pred(:,t), [nx,1]);
                    wip = reshape(w_cl(:,t:t+1), [2*nw,1]); 
                    x_pred(:,t+1) = Hankel_X_future*G*[xi;yip;wip];
                else
                    yip = reshape(y_cl(:,t:t+1), [2*ny,1]); 
                    xi = reshape(x_pred(:,t), [nx,1]);
        %             g = G*[xi;yip];
                    x_pred(:,t+1) = Hankel_X_future*G*[xi;yip];
                end      
    end
    for t = tem:T_tot
            if ifTemp
                u_pred(:,t-1) = B\(x_pred(:,t)-A*x_pred(:,t-1) - Bw*w_cl(:,t-1));
            else
                u_pred(:,t-1) = B\(x_pred(:,t)-A*x_pred(:,t-1)); %inv(B)*
    %             x_pred(:,t)-A*x_pred(:,t-1) - B*u_pred(t-1)
            end
    end
else
    for t = T_test_start-1:T_tot-1
                if ifU
                    yip = reshape(y_cl(:,t:t+1), [2*ny,1]); 
                    xi = reshape(x_pred(:,t), [nx,1]);
                    wip = reshape(w_cl(:,t:t+1), [2*nw,1]); 
                    x_pred(:,t+1) = Hankel_X_future*G*[xi;yip;wip];
                else
                    yip = reshape(y_cl(:,t:t+1), [2*ny,1]); 
                    xi = reshape(x_pred(:,t), [nx,1]);
        %             g = G*[xi;yip];
                    x_pred(:,t+1) = Hankel_X_future*G*[xi;yip];
                end      
    end
    if ifD
        for t = T_test_start:T_tot
            u_pred(:,t) = D\(y_cl(:,t)-C*x_pred(:,t));
        end
    else
        for t = T_test_start+1:T_tot
            if ifU
                u_pred(:,t-1) = B\(x_pred(:,t)-A*x_pred(:,t-1) - Bw*w_cl(:,t-1));
            else
                u_pred(:,t-1) = B\(x_pred(:,t)-A*x_pred(:,t-1)); %inv(B)*
    %             x_pred(:,t)-A*x_pred(:,t-1) - B*u_pred(t-1)
            end     
        end
    end
end
%Comupte the input

%%
num = 1;
figure()
hold on; grid on;
if ifExp
    % start = 97; last = T_train_end;
    start = 7*96+1; last = 12*96;
else
    start = T_test_start; last = T_tot;
end
plot(u_cl(num,start:last),'b--');
plot(u_pred(num,start:last),'r','LineWidth',1.5);
legend('Real people','Estimated people','FontSize',18) % 'Estimation point',
title('Estimate occpancy number in Polydome','fontsize',20, 'interpreter', 'latex')
mean(abs(u_pred(num,start:last)-u_cl(num,start:last)))

for i = start:last
    if mod(i,96)<=28 || mod(i,96)>=93  % not work time
        u_pred(1,i) = 0;
    end    
%     u_pred(1,i) = max(u_pred(1,i),0);
end
mean(abs(u_pred(1,start:last)-u_cl(1,start:last)))

% mean(abs(u_pred_basic(num,start:last)-u_cl(num,start:last)))
% CO2  '()': force 0 at night
%   m2=3 No de: Luca(wrong, Dsr * y_cl) 8.30 (5.62)  0night, 5.45 noIno 12.93 (8.53)  0night, 8.68
%               ifIno=1 14.6(10.0)11.5 ifIno=2 bad
%               m=7 ifIno=3 8.7(6.16)5.98 
%   m2 =3 De, right, Luca(wrong, Dsr * y_cl) 5.1576 (4.894)  0night, 4.875; ifIno=0 7.22(6.30) 6.27
%               ifIno=1 7.43(6.60)6.61 ifIno=2 bad
%                m=6 ifIno=3 5.98(5.51)5.53  m=6 ifIno=0 5.9(5.4)5.4
%               m=6 ifIno=1(W=eye) 5.4(5.0)5.0
%               m=7 ifIno=3 6.2(5.59)5.59  m=7 ifIno=1(W=eye) 5.1(4.7)4.7
%                m=8 ifIno=3 7.5(5.9)6.0
% CO2+temp
%   No de: m2 
 %     ifIno=3 13.6(10.4)11.4 ifIno=0 8.1(6.3)6.8 wrong at night
 %      m=2 ifIno=0 7.7(6.0)6.5 m=6 ifIno=0 9.4(7.3)7.4 
 %      m=6 ifIno=3 8.2(5.7)7.2
%   DeCo2:
%       m=3 ifIno=3 8.3(6.2)6.0  m=3ifIno=0 8.2(5.8)6.3
%       m=6 ifIno=0 8.6(6.6)6.7
%   DeAll:
%       m=3 ifIno=0 (6.2) m=3ifIno=1(W=eye) (6.1)
1
%%
u_pred = zeros(nu,T_tot);
z_pred = zeros(nx,T_tot);
x_pred = zeros(nx,T_tot); 
for i = 0:14
    for t = i*96+1:(i+1)*96
        if t+2>T_tot
            break;
        end
        if t<=i*96+28 || t>=(i+1)*96-4
            if ifTemp
                x_pred(:,t+1) = A*x_pred(:,t) + Bw*w_cl(:,t);
            else
                x_pred(:,t+1) = A*x_pred(:,t) ;
            end

        else
            if ifTemp
                yip = reshape(y_cl(:,t:t+1), [2*ny,1]); 
                xi = reshape(x_pred(:,t), [nx,1]);
                wip = reshape(w_cl(:,t:t+1), [2*nw,1]); 
                x_pred(:,t+1) = Hankel_X_future*G*[xi;yip;wip];                
            else
                yip = reshape(y_cl(:,t:t+1), [2*ny,1]); 
                xi = reshape(x_pred(:,t), [nx,1]);
                x_pred(:,t+1) = Hankel_X_future*G*[xi;yip];
            end

        end
           
    end
end

for t = 2:T_tot
        if ifTemp
            u_pred(:,t-1) = B\(x_pred(:,t)-A*x_pred(:,t-1) - Bw*w_cl(:,t-1));
        else
            u_pred(:, t-1) = B\(x_pred(:,t)-A*x_pred(:,t-1)); % inv(B)*
%             x_pred(:,t)-A*x_pred(:,t-1) - B*u_pred(t-1)
        end
end
mean(abs(u_pred(num,start:last)-u_cl(num,start:last)))
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