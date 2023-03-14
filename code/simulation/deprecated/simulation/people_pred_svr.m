% polydome_hankel_prediction_mystruct
clc
clear
rng(125)
%%
load('raw_06-12-2021_19-12-2021.mat')

% y_cl = mean([exp.sensor_temp_1,exp.sensor_temp_2,exp.sensor_temp_3,exp.sensor_temp_4],2)';
% y_cl = [mean([exp.sensor_temp_1,exp.sensor_temp_2,exp.sensor_temp_3,exp.sensor_temp_4],2)';
%         mean([exp.co2_1,exp.co2_2,exp.co2_3,exp.co2_4],2)'/100];
y_cl = mean([exp.co2_1,exp.co2_2,exp.co2_3,exp.co2_4],2)'/100;
% detrend
%     y_cl_env = mean([exp.co2_1,exp.co2_2,exp.co2_3,exp.co2_4]')/100;
%     y_cl_env = [y_cl_env(:,96*7+1:96*7+32), y_cl_env(:,96*6+33:96*7)];
%     y_cl_env = kron(ones(ny,14),y_cl_env);
%     y_cl_env = y_cl_env(1:size(y_cl,2));
%     y_cl = y_cl - y_cl_env;
temp_cl = mean([exp.sensor_temp_1,exp.sensor_temp_2,exp.sensor_temp_3,exp.sensor_temp_4],2)';

u_cl = exp.people';
w_cl = [exp.power-2.35,exp.weather_temp,exp.weather_rad/100]';

T = 96*4+50; % T >= (nu + 1)*(N_ini + N_pred + nx) - 1
T_tot = size(y_cl,2); % total steps for simulation
N_init = 4;
N_pred = 0;
%%
[x_train, y_train] = select_feature_co2(u_cl(:,1:96*4+50+N_pred),y_cl(:,1:96*4+50+N_pred),N_init,N_pred,0);
[x_test, y_test] = select_feature_co2(u_cl(:,7*96+7*4+1-N_init:end),y_cl(:,7*96+7*4+1-N_init:end),N_init,N_pred,0);
% [x_test, y_test] = select_feature_co2(u_cl(:,7*96+1-N_init:end),y_cl(:,7*96+1-N_init:end),N_init,N_pred,0);

% [x_train, y_train] = select_feature_temp(u_cl(:,1:96*4+50+N_pred),temp_cl(:,1:96*4+50+N_pred),w_cl(:,1:96*4+50+N_pred),N_init,N_pred,0);
% [x_test, y_test] = select_feature_temp(u_cl(:,7*96+7*4+1-N_init:end),temp_cl(:,7*96+7*4+1-N_init:end),w_cl(:,7*96+7*4+1-N_init:end),N_init,N_pred,0);
% [x_train, y_train] = select_feature_co2temp(u_cl(:,1:96*4+50+N_pred),y_cl(:,1:96*4+50+N_pred),temp_cl(:,1:96*4+50+N_pred),w_cl(:,1:96*4+50+N_pred),N_init,N_pred,0);
% [x_test, y_test] = select_feature_co2temp(u_cl(:,7*96+7*4+1-N_init:end),y_cl(:,7*96+7*4+1-N_init:end),temp_cl(:,7*96+7*4+1-N_init:end),w_cl(:,7*96+7*4+1-N_init:end),N_init,N_pred,0);
%%
% 1:96*4+50;
% i*96+7*4+1:N_pred:(i+1)*96
%%
Mdl = fitrsvm(x_train',y_train','KernelFunction','rbf','KernelScale','auto',... % rbf linear
    'Standardize',false);
y_pred = predict(Mdl,x_train');
figure()
plot(y_train,'b.');
hold on;
plot(y_pred,'r','LineWidth',1.5);
xlabel('time');
ylabel('y');
legend('Data','GPR predictions');
hold off
%%
% gprMdl2 = fitrgp(x_train',y_train','KernelFunction','squaredexponential',...
%     'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
%     struct('AcquisitionFunctionName','expected-improvement-plus'));
%%
y_test_pred = predict(Mdl,x_test');
figure()
start = 1; last = size(y_test,2)-96*2;
plot(y_test(start:last),'b--');
hold on;
plot(y_test_pred(start:last),'r','LineWidth',1.5);
xlabel('time');
ylabel('y');
legend('Data','GPR predictions');
hold off
% Only CO2:N_ini=4,N_pred=0,rbf (3.5594 use std, 3.5476), linear 4.3614
%       N_ini=4,N_pred=4,rbf 3.6880, linear 4.5820
%       N_ini=10,N_pred=4,rbf 4.2968, linear 4.6685
%       N_ini=10,N_pred=0,rbf 3.7250, linear 4.5310
% [CO2; diff]:  N_ini=4,N_pred=0,rbf 3.5672, linear 4.3500
%       N_ini=4,N_pred=4,rbf 3.8514, linear 4.5051
%       N_ini=10,N_pred=4,rbf 4.4158, linear 4.6063
%       N_ini=10,N_pred=0,rbf 3.8391, linear 4.5309
% [CO2; temperature;power;wea]
%       N_ini=4,N_pred=0,rbf 13.2418, linear 8.5056
%       N_ini=10,N_pred=4,rbf 12.4865, linear 9.5291
%       N_ini=10,N_pred=0,rbf 13.5532, linear 8.6242
% [temperature;power;wea]
%       N_ini=4,N_pred=0,rbf 11.4212, linear 5.1673
%       N_ini=10,N_pred=4,rbf 9.02, linear 6.5893
%       N_ini=10,N_pred=0,rbf 9.3247 linear 5.9937
mean(abs(y_test_pred(start:last)'-y_test(start:last)))
% No weekends
% Only CO2:N_ini=4,N_pred=0,rbf 4.8168 , linear  5.7890
%%
start = 1; last = 5*96;
u22 = y_test(start:last);
u2 = y_test_pred(start:last);
save('./huatu/svr.mat','u22','u2')
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

%%
num = 1;
figure()
hold on; grid on;
start = 7*96+1; last = size(u_pred,2);
plot( u_cl(num,start:last), 'Color',[0 0.4470 0.7410],'LineWidth',1,'Marker','s')
plot( u_pred(num,start:last), 'Color',[0.8500 0.3250 0.0980],'LineWidth',1,'Marker','o')
% plot( u_pred_basic(num,start:last), 'Color',[0.4660 0.6740 0.188],'LineWidth',1)
% scatter( [(start:N_pred:last)-start+1],u_pred(num,start:N_pred:last), 100,'red')
legend('Real people','Estimated people','FontSize',18) % 'Estimation point',

title('Estimate occpancy number in Polydome','fontsize',20, 'interpreter', 'latex')
mean(abs(u_pred(num,start:last)-u_cl(num,start:last)))

%% 
function [x,y] = select_feature_co2(people,co2, N_init,N_pred, diff)
    % At time t, use current CO2, N_init previous and N_pred future steps of CO2 level as x(t), 
    % use people(t) as y(t)
    len = size(co2,2);
    i = N_init + 1;
    while(i+N_pred<=len)
        x(:,i-N_init) = reshape(co2(:,i-N_init:i+N_pred),[],1);
        y(:,i-N_init) = people(:,i-N_init);
        i = i + 1;
    end
    % add the co2 difference?
    if diff == 1
        for i = 1:size(x,2)
            tem(:,i) = x(2:end,i) - x(1:end-1,i);
        end
        x = [x;tem];
    end
    
end
%%
function [x,y] = select_feature_temp(people,temp,other_input,N_init,N_pred, diff)
    % At time t, use current CO2, N_init previous and N_pred future steps of CO2 level as x(t), 
    % use people(t) as y(t)
    len = size(temp,2);
    i = N_init + 1;
    while(i+N_pred<=len)
         tem2 = reshape(temp(:,i-N_init:i+N_pred),[],1);
         tem3 = reshape(other_input(:,i-N_init:i+N_pred),[],1);
         x(:,i-N_init) = [tem2;tem3];
        y(:,i-N_init) = people(:,i-N_init);
        i = i + 1;
    end
    % add the co2 difference?
    if diff == 1
        for i = 1:size(x,2)
            tem(:,i) = x(2:end,i) - x(1:end-1,i);
        end
        x = [x;tem];
    end
    
end
%%
function [x,y] = select_feature_co2temp(people,co2,temp,other_input,N_init,N_pred, diff)
    % At time t, use current CO2, N_init previous and N_pred future steps of CO2 level as x(t), 
    % use people(t) as y(t)
    len = size(co2,2);
    i = N_init + 1;
    while(i+N_pred<=len)
         tem1 = reshape(co2(:,i-N_init:i+N_pred),[],1);
         tem2 = reshape(temp(:,i-N_init:i+N_pred),[],1);
         tem3 = reshape(other_input(:,i-N_init:i+N_pred),[],1);
         x(:,i-N_init) = [tem1;tem2;tem3];
        y(:,i-N_init) = people(:,i-N_init);
        i = i + 1;
    end
    % add the co2 difference?
    if diff == 1
        for i = 1:size(x,2)
            tem(:,i) = x(2:end,i) - x(1:end-1,i);
        end
        x = [x;tem];
    end
    
end