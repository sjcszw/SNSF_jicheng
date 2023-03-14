%%
clc
clear; 
% close all;
% addpath(genpath('./utilities'))

%%
% exp3 = exp1;
exp3.time = exp1.power.time;
exp3.sensor_temp_1 = exp1.sensor_temp.value(:,1);
exp3.sensor_temp_2 = exp1.sensor_temp.value(:,2);
exp3.sensor_temp_3 = exp1.sensor_temp.value(:,3);
exp3.sensor_temp_4 = exp1.sensor_temp.value(:,4);
exp3.power = exp1.power.value;
exp3.supply_temp = exp1.supply_temp.value;
exp3.return_temp = exp1.return_temp.value;
exp3.supply_flow = exp1.supply_flow.value;
exp3.weather_temp = exp1.air_temp.value;
exp3.weather_rad = exp1.solar_GHI.value;
exp3.co2 = exp1.co2.value;

exp3.mode = exp1.mode.value;
exp3.setpoint_cool = exp1.setpoint_winter.value;
exp3.setpoint_heat = exp1.setpoint_summer.value;

%%%%%%%%%%
exp4.time = exp2.power.time;
exp4.sensor_temp_1 = exp2.sensor_temp.value(:,1);
exp4.sensor_temp_2 = exp2.sensor_temp.value(:,2);
exp4.sensor_temp_3 = exp2.sensor_temp.value(:,3);
exp4.sensor_temp_4 = exp2.sensor_temp.value(:,4);
exp4.power = exp2.power.value;
exp4.supply_temp = exp2.supply_temp.value;
exp4.return_temp = exp2.return_temp.value;
exp4.supply_flow = exp2.supply_flow.value;
exp4.weather_temp = exp2.air_temp.value;
exp4.weather_rad = exp2.solar_GHI.value;
exp4.co2 = exp2.co2.value;

exp4.mode = exp2.mode.value;
exp4.setpoint_cool = exp2.setpoint_winter.value;
exp4.setpoint_heat = exp2.setpoint_summer.value;
%%
exp3.time = [exp3.time; exp1.power.time];
exp3.sensor_temp_1 = [exp3.sensor_temp_1; exp1.sensor_temp.value(:,1)];
exp3.sensor_temp_2 = [exp3.sensor_temp_2; exp1.sensor_temp.value(:,2)];
exp3.sensor_temp_3 = [exp3.sensor_temp_3; exp1.sensor_temp.value(:,3)];
exp3.sensor_temp_4 = [exp3.sensor_temp_4; exp1.sensor_temp.value(:,4)];
exp3.power = [exp3.power; exp1.power.value];
exp3.supply_temp = [exp3.supply_temp; exp1.supply_temp.value];
exp3.return_temp = [exp3.return_temp ; exp1.return_temp.value];
exp3.supply_flow = [exp3.supply_flow; exp1.supply_flow.value];
exp3.weather_temp = [exp3.weather_temp; exp1.air_temp.value];
exp3.weather_rad = [exp3.weather_rad; exp1.solar_GHI.value];
exp3.co2 = [exp3.co2; exp1.co2.value];

exp3.mode = [exp3.mode; exp1.mode.value];
exp3.setpoint_cool = [exp3.setpoint_cool; exp1.setpoint_winter.value];
exp3.setpoint_heat = [exp3.setpoint_heat; exp1.setpoint_summer.value];

%%%%%%
exp4.time = [exp4.time; exp2.power.time];
exp4.sensor_temp_1 = [exp4.sensor_temp_1; exp2.sensor_temp.value(:,1)];
exp4.sensor_temp_2 = [exp4.sensor_temp_2; exp2.sensor_temp.value(:,2)];
exp4.sensor_temp_3 = [exp4.sensor_temp_3; exp2.sensor_temp.value(:,3)];
exp4.sensor_temp_4 = [exp4.sensor_temp_4; exp2.sensor_temp.value(:,4)];
exp4.power = [exp4.power; exp2.power.value];
exp4.supply_temp = [exp4.supply_temp; exp2.supply_temp.value];
exp4.return_temp = [exp4.return_temp ; exp2.return_temp.value];
exp4.supply_flow = [exp4.supply_flow; exp2.supply_flow.value];
exp4.weather_temp = [exp4.weather_temp; exp2.air_temp.value];
exp4.weather_rad = [exp4.weather_rad; exp2.solar_GHI.value];
exp4.co2 = [exp4.co2; exp2.co2.value];

exp4.mode = [exp4.mode; exp2.mode.value];
exp4.setpoint_cool = [exp4.setpoint_cool; exp2.setpoint_winter.value];
exp4.setpoint_heat = [exp4.setpoint_heat; exp2.setpoint_summer.value];
%% change power of mode: cooling, -power
for i = 1:96*6*5
    if exp3.mode(i) ==  0
        exp3.power(i) = -exp3.power(i);
    end
    if exp4.mode(i) ==  0
        exp4.power(i) = -exp4.power(i);
    end    
end
%% check time
for i = 0:4
    datestr(exp3.time(i*96*6+1))
end
%%
figure(); hold on;
plot(exp3.power-exp4.power);
% plot(exp3.power);
plot(10*exp3.mode);
%%
figure(); hold on;
% plot(exp3.co2(:,1)-exp4.co2(:,1));
plot(exp3.co2(:,1));
plot(exp4.co2(:,1));
plot(10*exp3.mode);
%%
figure(); hold on;
plot(exp3.sensor_temp_1(:,1),'g');
plot(exp3.sensor_temp_2(:,1),'r');
plot(exp3.sensor_temp_3(:,1),'b');
plot(10*exp3.mode);
%% save .mat
exp = exp3;
exp.time_str = datestr(exp.time);
save('./raw_2021-11-12_2021-11-23.mat','exp')
% save('./data_07/raw_2021_1112_2021-1118.mat','exp')
%%
exp = exp4;
exp.time_str = datestr(exp.time);
save('./clear_2021-11-12_2021-11-23.mat','exp')
%%
exp = exp3;
% exp = exp4;
ddd = table;
ddd.time_str = datestr(exp.time);
ddd.time = exp.time;
ddd.sensor_temp_1 = exp.sensor_temp_1;
ddd.sensor_temp_2 = exp.sensor_temp_2;
ddd.sensor_temp_3 = exp.sensor_temp_3;
ddd.sensor_temp_4 = exp.sensor_temp_4;
ddd.power = exp.power;
ddd.supply_temp = exp.supply_temp;
ddd.return_temp = exp.return_temp;
ddd.supply_flow = exp.supply_flow;
ddd.weather_temp = exp.weather_temp;
ddd.weather_rad = exp.weather_rad;

ddd.mode = exp.mode;
ddd.setpoint_cool = exp.setpoint_cool;
ddd.setpoint_heat = exp.setpoint_heat;
writetable(ddd, './database/raw_2021-07-15_2021-08-15.csv', 'Delimiter',',');






