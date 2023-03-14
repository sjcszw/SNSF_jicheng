%% 16052022
% 1. System id: check model with training and validation data
% 2. Design UIO (test linear system)

%% Save u1 D~=0 u2 D=0 p3 UIO
% start = T_test_start-5; last = T_tot-1;
% u2_p3 = u_pred(:,start:last); 
% u2_r3 = u_cl(:,start:last); % r1: real u
% save('./u2_p3_1805.mat','u2_p3','u2_r3');
% start = T_test_start-5; last = T_tot;
% u1_p3 = u_pred(:,start:last); 
% u1_r3 = u_cl(:,start:last); % r1: real u
% save('./u1_p3_1805.mat','u1_p3','u1_r3');
% start = 7*96+1; last = 12*96;
% ur_UIO = u_cl(num,start:last);
% u_UIO = u_pred(num,start:last);
% save('./UIO.mat','ur_UIO','u_UIO')
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
        w_cl = zeros(nw,T_tot)-1;        
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

if ifExp
    pred_step = 4;
    opt.Focus = 'prediction'; 
    opt.N4Horizon = 1*3600/T_samp;
    order = 2;
    [m1, X0] = ssest(exp_iden, order, 'Ts', T_samp, opt);
    order = 3;
    [m2, X0] = ssest(exp_iden, order, 'Ts', T_samp, opt);
    order = 8;
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
%% ====== 2. Design UIO ======
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


%% --- 2.2 Check conditions, compute N L E
rng(123)
if ~ifD
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
    %
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
else
    E = zeros(nx,ny);
    L = B/D; %*inv(D)
    N = A - L*C;
end
%% --- 2.3 estimate input
u_pred = zeros(nu,T_tot);
z_pred = zeros(nx,T_tot);
x_pred = zeros(nx,T_tot);
if ifExp
    x_pred(:,1) = z_pred(:,1) + E*y_cl(:,1); 
    tem = 2;%T_test_start;
    for t = tem:T_tot
                if ifTemp
                    z_pred(:,t) = N*z_pred(:,t-1) + L*y_cl(:,t-1) + (eye(nx)-E*C)*Bw*w_cl(:,t-1);
                else
                    z_pred(:,t) = N*z_pred(:,t-1) + L*y_cl(:,t-1);
                end
                x_pred(:,t) = z_pred(:,t) + E*y_cl(:,t);     
    end
    for t = tem:T_tot
            if ifTemp || ifU
                u_pred(:,t-1) = B\(x_pred(:,t)-A*x_pred(:,t-1) - Bw*w_cl(:,t-1));
            else
                u_pred(:,t-1) = B\(x_pred(:,t)-A*x_pred(:,t-1)); %inv(B)*
    %             x_pred(:,t)-A*x_pred(:,t-1) - B*u_pred(t-1)
            end
    end    
else
    for t = T_test_start-1:T_tot
                if ifU
                    z_pred(:,t) = N*z_pred(:,t-1) + L*y_cl(:,t-1) + (eye(nx)-E*C)*Bw*w_cl(:,t-1);
                else
                    z_pred(:,t) = N*z_pred(:,t-1) + L*y_cl(:,t-1);
                end
                x_pred(:,t) = z_pred(:,t) + E*y_cl(:,t);  
    end
    %Comupte the input
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


%%
% for i = start:last
%     if mod(i,96)<=28 || mod(i,96)>=93  % not work time
%         u_pred(1,i) = 0;
%     end
%     u_pred(1,i) = max(u_pred(1,i),0);       
%     if u_pred(1,i)>70
%         u_pred(1,i) = mean(u_pred(1,i-1:i));
%     end    
% end
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
    if u_pred(1,i)>80 || u_pred(1,i)<0
        u_pred(1,i) = mean(u_pred(1,i-1:i));
    end    
end
mean(abs(u_pred(1,start:last)-u_cl(1,start:last)))

% mean(abs(u_pred_basic(num,start:last)-u_cl(num,start:last)))
% CO2 m2
%   No de: 12.2048 0night, 9.840 % delete neg: 8.54  smooth: 7.70
%   De, right, 9.936(8.5) 0night, 8.39 % delete neg: 6.42 smooth: 5.85
%       m=2 9.0(7.5)7.4; m=7 12.0(9.4)6.3
% CO2+temp
%   No de: K0 m2 8.4665 0night, 11.9 % delete neg: 5.17 smooth: 4.83
%            K0  m1 14.0 0night, 11.5 % delete neg: 5.23 smooth: 4.79
%   De: 
%           deCo2, m1 15.484 0night,11.822 delete neg: 5.357
%           deAll, m1 8.986 0night, 8.52 delete neg: 5.84 smooth: 5.32
1
%%
u_pred = zeros(nu,T_tot);
z_pred = zeros(nx,T_tot);
x_pred = zeros(nx,T_tot);
x_pred(:,1) = z_pred(:,1) + E*y_cl(:,1); 
for i = 0:14
    z_pred(:,i*96+1)=0; % 1705 morning has this, wrong; but no influence?
    for t = i*96+1:(i+1)*96
        if t+2>T_tot
            break;
        end
        if t<=i*96+27 || t>=(i+1)*96-4 %!!!!!! should be<=28 because u(t) >=93
            if ifTemp
                x_pred(:,t+1) = A*x_pred(:,t) + Bw*w_cl(:,t);
            else
                x_pred(:,t+1) = A*x_pred(:,t) ;
            end
            z_pred(:,t+1) = x_pred(:,t+1) - E*y_cl(:,t+1);
        else
            if ifTemp
                z_pred(:,t+1) = N*z_pred(:,t) + L*y_cl(:,t) + (eye(nx)-E*C)*Bw*w_cl(:,t-1);
            else
                z_pred(:,t+1) = N*z_pred(:,t) + L*y_cl(:,t);
            end
            x_pred(:,t+1) = z_pred(:,t+1) + E*y_cl(:,t+1);
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
start = 7*96+1; last = 12*96;
mean(abs(u_pred(1,start:last)-u_cl(num,start:last)))
for i = start:last
    u_pred(1,i) = max(u_pred(1,i),0);
    if u_pred(1,i)>80
        u_pred(1,i) = mean(u_pred(1,i-1:i));
    end
end
mean(abs(u_pred(1,start:last)-u_cl(num,start:last)))
figure()
hold on; grid on;
plot(u_cl(1,start:last),'b--');
plot(u_pred(1,start:last),'r','LineWidth',1.5);


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