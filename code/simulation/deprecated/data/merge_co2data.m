%% local time zone: because we want to combine people, make sure the time slot is correct
clc
clear; 
% close all;
% addpath(genpath('./utilities'))
%%
% save('people_1206-1219.mat','people')
% people = reshape(people,[],1);
% save('people_1206-1219_vector.mat','people')
%%
% exp3 = exp1;
exp3.time = exp1.co2.time;
exp3.co2_1 = exp1.co2.value(:,1);
exp3.co2_2 = exp1.co2.value(:,2);
exp3.co2_3 = exp1.co2.value(:,3);
exp3.co2_4 = exp1.co2.value(:,4);
exp3.supply_flow = exp1.supply_flow.value;

exp3.sensor_temp_1 = exp1.sensor_temp.value(:,1);
exp3.sensor_temp_2 = exp1.sensor_temp.value(:,2);
exp3.sensor_temp_3 = exp1.sensor_temp.value(:,3);
exp3.sensor_temp_4 = exp1.sensor_temp.value(:,4);
exp3.power = exp1.power.value;
exp3.weather_temp = exp1.air_temp.value;
exp3.weather_rad = exp1.solar_GHI.value;

%%%%%%%%%%
exp4.time = exp2.co2.time;
exp4.co2_1 = exp2.co2.value(:,1);
exp4.co2_2 = exp2.co2.value(:,2);
exp4.co2_3 = exp2.co2.value(:,3);
exp4.co2_4 = exp2.co2.value(:,4);
exp4.supply_flow = exp2.supply_flow.value;

exp4.sensor_temp_1 = exp2.sensor_temp.value(:,1);
exp4.sensor_temp_2 = exp2.sensor_temp.value(:,2);
exp4.sensor_temp_3 = exp2.sensor_temp.value(:,3);
exp4.sensor_temp_4 = exp2.sensor_temp.value(:,4);
exp4.power = exp2.power.value;
exp4.weather_temp = exp2.air_temp.value;
exp4.weather_rad = exp2.solar_GHI.value;

%%
exp3.time = [exp3.time; exp1.power.time];
exp3.co2_1 = [exp3.co2_1; exp1.co2.value(:,1)];
exp3.co2_2 = [exp3.co2_2; exp1.co2.value(:,2)];
exp3.co2_3 = [exp3.co2_3; exp1.co2.value(:,3)];
exp3.co2_4 = [exp3.co2_4; exp1.co2.value(:,4)];
exp3.supply_flow = [exp3.supply_flow; exp1.supply_flow.value];

exp3.sensor_temp_1 = [exp3.sensor_temp_1; exp1.sensor_temp.value(:,1)];
exp3.sensor_temp_2 = [exp3.sensor_temp_2; exp1.sensor_temp.value(:,2)];
exp3.sensor_temp_3 = [exp3.sensor_temp_3; exp1.sensor_temp.value(:,3)];
exp3.sensor_temp_4 = [exp3.sensor_temp_4; exp1.sensor_temp.value(:,4)];
exp3.power = [exp3.power; exp1.power.value];
exp3.weather_temp = [exp3.weather_temp; exp1.air_temp.value];
exp3.weather_rad = [exp3.weather_rad; exp1.solar_GHI.value];

%%%%%%
exp4.time = [exp4.time; exp2.co2.time];
exp4.co2_1 = [exp4.co2_1; exp2.co2.value(:,1)];
exp4.co2_2 = [exp4.co2_2; exp2.co2.value(:,2)];
exp4.co2_3 = [exp4.co2_3; exp2.co2.value(:,3)];
exp4.co2_4 = [exp4.co2_4; exp2.co2.value(:,4)];
exp4.supply_flow = [exp4.supply_flow; exp2.supply_flow.value];

exp4.sensor_temp_1 = [exp4.sensor_temp_1; exp2.sensor_temp.value(:,1)];
exp4.sensor_temp_2 = [exp4.sensor_temp_2; exp2.sensor_temp.value(:,2)];
exp4.sensor_temp_3 = [exp4.sensor_temp_3; exp2.sensor_temp.value(:,3)];
exp4.sensor_temp_4 = [exp4.sensor_temp_4; exp2.sensor_temp.value(:,4)];
exp4.power = [exp4.power; exp2.power.value];
exp4.weather_temp = [exp4.weather_temp; exp2.air_temp.value];
exp4.weather_rad = [exp4.weather_rad; exp2.solar_GHI.value];
%%
% change the time zone to Zurich 
exp3.time = exp3.time + 1/24;
exp4.time = exp4.time + 1/24;
%% check time
for i = 0:10
    datestr(exp3.time(i*96*6+1))
end
%%
figure(); hold on;
% plot(exp3.power-exp4.power);
plot(exp3.power,'r');
plot(exp4.power,'g');
% plot(exp3.co2_4,'r');
% plot(exp4.co2_4,'g');
% plot(people*10)
%% save .mat
exp = exp3;
exp.people = people;
exp.time_str = datestr(exp.time);
save('./merged/raw_06-12-2021_19-12-2021.mat','exp')
% save('./data_07/raw_2021_1112_2021-1118.mat','exp')
%%
exp = exp4;
exp.people = people;
exp.time_str = datestr(exp.time);
save('./merged/clear_06-12-2021_19-12-2021.mat','exp')
%%
exp = exp3;
% exp = exp4;
ddd = table;
ddd.time_str = datestr(exp.time);
ddd.time = exp.time;
ddd.co2_1 = exp.co2_1;
ddd.co2_2 = exp.co2_2;
ddd.co2_3 = exp.co2_3;
ddd.co2_4 = exp.co2_4;
ddd.supply_flow = exp.supply_flow;

writetable(ddd, './merged/raw_06-12-2021_14-12-2021.csv', 'Delimiter',',');
% writetable(ddd, './merged/clear_06-12-2021_14-12-2021.csv', 'Delimiter',',');