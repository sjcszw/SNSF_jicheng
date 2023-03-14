
y1 = record_ada.y_cl;
u1 = record_ada.u_cl;
y11 = record_L2.y_cl;
u11 = record_L2.u_cl;

y2 = record_ada.y_cl;
u2 = record_ada.u_cl;
y21 = record_L2.y_cl;
u21 = record_L2.u_cl;

y3 = record_ada.y_cl;
u3 = record_ada.u_cl;
y31 = record_L2.y_cl;
u31 = record_L2.u_cl;

y4 = record_ada.y_cl;
u4 = record_ada.u_cl;
y41 = record_L2.y_cl;
u41 = record_L2.u_cl;

y5 = record_ada.y_cl;
u5 = record_ada.u_cl;
y51 = record_L2.y_cl;
u51 = record_L2.u_cl;

%%
record_ada.u_cl = [ u1,u2];
record_ada.y_cl = [ y1,y2];
record_L2.u_cl = [ u11,u21];
record_L2.y_cl = [ y11,y21];
record_ada.u_cl = [ u1,u2,u3,u4,u5];
record_ada.y_cl = [ y1,y2,y3,y4,y5];
record_L2.u_cl = [ u11,u21,u31,u41,u51];
record_L2.y_cl = [ y11,y21,y31,y41,y51];
save('simall_144_Q10R0.01_N6.mat','record_ada','record_L2','Q' ,'R', 'r','T')

%%
y_cl = record_ada.y_cl{1};
u_cl = record_ada.u_cl{1};
for num_tem = 2:10
    y_cl = y_cl + record_ada.y_cl{num_tem};
    u_cl = u_cl + record_ada.u_cl{num_tem};
end
y_cl = y_cl/num_tem;
u_cl = u_cl/num_tem;
%% ======read
num_mon = 40;
y_cl(1,:) = record_ada.y_cl{1}(145:216);
u_cl(1,:) = record_ada.u_cl{1}(145:216);
for num_tem = 2:num_mon
    y_cl(num_tem,:) = record_ada.y_cl{num_tem}(145:216);
    u_cl(num_tem,:) = record_ada.u_cl{num_tem}(145:216);
end
y_std = std(y_cl);
y_ada = mean(y_cl);
u_std = std(u_cl);
u_ada = mean(u_cl);
y1_ada = y_ada+y_std;
y2_ada = y_ada-y_std;
u1_ada = u_ada+u_std;
u2_ada = u_ada-u_std;
y_cl_ada = y_cl; u_cl_ada = u_cl;
%
y_cl(1,:) = record_L2.y_cl{1}(145:216);
u_cl(1,:) = record_L2.u_cl{1}(145:216);
for num_tem = 2:num_mon
    y_cl(num_tem,:) = record_L2.y_cl{num_tem}(145:216);
    u_cl(num_tem,:) = record_L2.u_cl{num_tem}(145:216);
end
y_std = std(y_cl);
y_L2 = mean(y_cl);
u_std = std(u_cl);
u_L2 = mean(u_cl);
y1_L2 = y_L2+y_std;
y2_L2 = y_L2-y_std;
u1_L2 = u_L2+u_std;
u2_L2 = u_L2-u_std;
y_cl_L2 = y_cl; u_cl_L2 = u_cl;
%%   ===== plot

%
a = datenum(2021, 5, 9, 0, 0, 0):(1200/86400):datenum(2021, 5, 9, 23, 40, 0);
figure
hold on


patch([a fliplr(a)], [y1_mpc fliplr(y2_mpc)], [0.8 1 0.8])
patch([a fliplr(a)], [y1_ada fliplr(y2_ada)], [0.9 1 1])

plot(a,y_ada,'-sb');
plot(a,y_mpc,'-sg');
plot(a,[20*ones(1,24),22*ones(1,30),20*ones(1,18)],'r');
tem = [30*ones(1,24), 26*ones(1,30), 30*ones(1,18)];
plot(a,tem,'r')
tem = [19*ones(1,24), 21*ones(1,30), 19*ones(1,18)];
plot(a,tem,'r')



% h = legend('y: closed loop', 'reference');
% set(h,'fontsize',10)
datetick('x','mm/dd/yy HH:MM');
%%
figure
hold on
figure
hold on


patch([a fliplr(a)], [u1_L2 fliplr(u2_L2)], [0.8 1 0.8])
patch([a fliplr(a)], [u1_ada fliplr(u2_ada)], [0.9 1 1])

plot(a,u_ada,'-sb');
plot(a,u_L2,'-sg');
%% ====== cost
display('ada');
mean(u_ada)
mean(y_ada)
cost = 0;
for i = 1:num_mon
    tem = compute_costQ( y_cl_ada(i, :),u_cl_ada(i, :),Q ,R, r(:,145:216), 72);
    cost = cost + tem(1);
end
cost= cost/72/num_mon
bb = reshape(u_cl_ada,1,[]);
stdu = std(bb)
bb = reshape(y_cl_ada,1,[]);
stdy = std(bb)

display('L2');
mean(u_L2)
mean(y_L2)
cost = 0;
for i = 1:num_mon
    tem = compute_costQ( y_cl_L2(i, :),u_cl_L2(i, :),Q ,R, r(:,145:216), 72);
    cost = cost + tem(1);
end
cost= cost/72/num_mon

%% ========= std

for i = 1:num_mon
    std_L2(1,i) = std(u_cl_L2(i, :)-u_L2);
    std_L2(2,i) = std(y_cl_L2(i, :)-y_L2);
    std_ada(1,i) = std(u_cl_ada(i, :)-u_ada);
    std_ada(2,i) = std(y_cl_ada(i, :)-y_ada);
end
disp('L2')
u = mean(std_L2(1,i))
y = mean(std_L2(2,i))
disp('ada')
u = mean(std_ada(1,i))
y = mean(std_ada(2,i))

% for i = 1:num_mon
%     std_L2(1,i) = std(u_cl_L2(i, :));
%     std_L2(2,i) = std(y_cl_L2(i, :));
%     std_ada(1,i) = std(u_cl_ada(i, :));
%     std_ada(2,i) = std(y_cl_ada(i, :));
% end
% disp('L2')
% u = mean(std_L2(1,i))
% y = mean(std_L2(2,i))
% disp('ada')
% u = mean(std_ada(1,i))
% y = mean(std_ada(2,i))
u_cl_L2d = u_cl_L2(:,2:end)-u_L2(:,2:end) - u_cl_L2(:,1:end-1)+u_L2(:,1:end-1);
y_cl_L2d = y_cl_L2(:,2:end)-y_L2(:,2:end) - y_cl_L2(:,1:end-1)+y_L2(:,1:end-1);
u_cl_adad = u_cl_ada(:,2:end)-u_ada(:,2:end) -u_cl_ada(:,1:end-1)+u_ada(:,1:end-1);
y_cl_adad = y_cl_ada(:,2:end)-y_ada(:,2:end) -y_cl_ada(:,1:end-1)+y_ada(:,1:end-1);
for i = 1:num_mon
    std_L2d(1,i) = std(u_cl_L2d(i, :));
    std_L2d(2,i) = std(y_cl_L2d(i, :));
    std_adad(1,i) = std(u_cl_adad(i, :));
    std_adad(2,i) = std(y_cl_adad(i, :));
end
disp('delta')
disp('L2')
u = mean(std_L2d(1,i))
y = mean(std_L2d(2,i))
disp('ada')
u = mean(std_adad(1,i))
y = mean(std_adad(2,i))
%% ========write
date_list = datestr(a, 'yyyy-mm-dd HH:MM');
ddd = table;
% ddd.tt = vec(cache.cl.time(1:end-1));
ddd.tt = date_list;
ddd.t = vec([1:72]);

ddd.max = vec([30*ones(1,24),26*ones(1,30),30*ones(1,18)]);
ddd.min = vec([19*ones(1,24),21*ones(1,30),19*ones(1,18)]);
ddd.ref = vec([20*ones(1,24),22*ones(1,30),20*ones(1,18)]);

ddd.y_ada = vec(y_ada);
ddd.u_ada = vec(u_ada);
ddd.y_L2 = vec(y_L2);
ddd.u_L2 = vec(u_L2);
ddd.y1_ada = vec(y1_ada);
ddd.y2_ada = vec(y2_ada);
ddd.u1_ada = vec(u1_ada);
ddd.u2_ada = vec(u2_ada);

ddd.y1_L2 = vec(y1_L2);
ddd.y2_L2 = vec(y2_L2);
ddd.u1_L2 = vec(u1_L2);
ddd.u2_L2 = vec(u2_L2);
%%
writetable(ddd, 'building_sim_regulation.dat', 'Delimiter',',');

%%
ddd = buildingsimminimizeu;
ddd.y_mpc = vec(y_mpc);
ddd.u_mpc = vec(u_mpc);
ddd.y1_mpc = vec(y1_mpc);
ddd.y2_mpc = vec(y2_mpc);
ddd.u1_mpc = vec(u1_mpc);
ddd.u2_mpc = vec(u2_mpc);
writetable(ddd, 'building_sim_minimizeu.dat', 'Delimiter',',');
%% ======read
num_mon = 30;
y_cl(1,:) = record_mpc.y_cl{1}(145:216);
u_cl(1,:) = record_mpc.u_cl{1}(145:216);
for num_tem = 2:num_mon
    y_cl(num_tem,:) = record_mpc.y_cl{num_tem}(145:216);
    u_cl(num_tem,:) = record_mpc.u_cl{num_tem}(145:216);
end
y_std = std(y_cl);
y_mpc = mean(y_cl);
u_std = std(u_cl);
u_mpc = mean(u_cl);
y1_mpc = y_mpc+y_std;
y2_mpc = y_mpc-y_std;
u1_mpc = u_mpc+u_std;
u2_mpc = u_mpc-u_std;
y_cl_mpc = y_cl; u_cl_mpc = u_cl;
%%
for i = 1:num_mon
    std_mpc(1,i) = std(u_cl_mpc(i, :)-u_mpc);
    std_mpc(2,i) = std(y_cl_mpc(i, :)-y_mpc);
end
disp('mpc')
u = mean(std_mpc(1,i))
y = mean(std_mpc(2,i))

u_cl_mpcd = u_cl_mpc(:,2:end)-u_mpc(:,2:end) - u_cl_mpc(:,1:end-1)+u_mpc(:,1:end-1);
y_cl_mpcd = y_cl_mpc(:,2:end)-y_mpc(:,2:end) - y_cl_mpc(:,1:end-1)+y_mpc(:,1:end-1);
for i = 1:num_mon
    std_mpcd(1,i) = std(u_cl_mpcd(i, :));
    std_mpcd(2,i) = std(y_cl_mpcd(i, :));
end
disp('delta')
disp('mpc')
u = mean(std_mpcd(1,i))
y = mean(std_mpcd(2,i))

disp('cost')
mean(mean(u_cl_mpc))

%%
ddd =table;
ddd.u_ada_reg = mean(u_cl(73:144));
ddd.y_ada_reg = mean(y_cl(73:144));
cost = compute_costQ( y_cl(:, 73:144),u_cl(:, 73:144),Q ,R, r(:,73:144), 72);
ddd.cost_ada_reg = mean(cost(1));
%%
ddd.u_ada_min = mean(u_cl(73:144));
ddd.y_ada_min = mean(y_cl(73:144));
cost = compute_costQ( y_cl(:, 73:144),u_cl(:, 73:144),Q ,R, r(:,73:144), 72);
ddd.cost_ada_min = mean(cost(1));
%%
ddd.u_L2_reg = mean(u_cl(73:144));
ddd.y_L2_reg = mean(y_cl(73:144));
cost = compute_costQ( y_cl(:, 73:144),u_cl(:, 73:144),Q ,R, r(:,73:144), 72);
ddd.cost_L2_reg = mean(cost(1));
%%
ddd.u_L2_min = mean(u_cl(73:144));
ddd.y_L2_min = mean(y_cl(73:144));
cost = compute_costQ( y_cl(:, 73:144),u_cl(:, 73:144),Q ,R, r(:,73:144), 72);
ddd.cost_L2_min = mean(cost(1));
%%
mean(u_cl(73:144))
mean(y_cl(73:144))
cost = compute_costQ( y_cl(:, 73:144),u_cl(:, 73:144),Q ,R, r(:,73:144), 72);
cost(1)/72
%% ============== mean std write
writetable(buildingsimmean, 'building_sim_mean.dat', 'Delimiter',';');

%% ================= disturbances
w1 = record_ada.w_real;
w2 = record_ada.w_real;
w3 = record_ada.w_real;
w4 = record_ada.w_real;
w5 = record_ada.w_real;
w_real = [ w1,w2,w3,w4,w5];
save('disturbance_mu','w_real')
%%
figure()
hold on
grid on

b = plot(mean(w1(:,145:216)), '--m');
bb  = plot(mean(w2(:,145:216)), '--b');
xlabel('Time t', 'interpreter', 'latex','fontsize',20);
ylabel('Disturbance ', 'interpreter', 'latex','fontsize',20);
bbb  = plot( mean(w3(:,145:216)), '--g');

h = legend([b,bb,  bbb],   'w1_{real}', 'w2_{real}', 'w3_{real}','Location','northeast');
set(h,'fontsize',24, 'interpreter', 'tex')

%% ============ read
num_mon = 50;
w1(1,:) = w_real{1}(1,:);
w2(1,:) = w_real{1}(2,:);
w3(1,:) = w_real{1}(3,:);
for num_tem = 2:num_mon
    w1(num_tem,:) = w_real{num_tem}(1,:);
    w2(num_tem,:) = w_real{num_tem}(2,:);
    w3(num_tem,:) = w_real{num_tem}(3,:);    
end

%% ========= write
a = datenum(2021, 5, 9, 0, 0, 0):(1200/86400):datenum(2021, 5, 9, 23, 40, 0);
date_list = datestr(a, 'yyyy-mm-dd HH:MM');
ddd =table;
% ddd.tt = vec(cache.cl.time(1:end-1));
ddd.tt = date_list;
ddd.t = vec([1:72]);

ddd.w1 = vec(mean(w1(:,145:216)));
ddd.w2 = vec(mean(w2(:,145:216)));
ddd.w3 = vec(mean(w3(:,145:216)));

%%
writetable(ddd, 'building_disturbance.dat', 'Delimiter',',');