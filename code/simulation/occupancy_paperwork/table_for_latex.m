%% L-CSS version
u11-u33
u11-u44
u11-ur_UIO
u11-ur_DUIO
%% 8am-7pm
for i = 0:4
    u1(i*96+1:i*96+32) = 0;
    u1(i*96+77:i*96+96) = 0;   
    u3(i*96+1:i*96+32) = 0;
    u3(i*96+77:i*96+96) = 0;
    u4(i*96+1:i*96+32) = 0;
    u4(i*96+77:i*96+96) = 0;   
    u_UIO(i*96+1:i*96+32) = 0;
    u_UIO(i*96+77:i*96+96) = 0;
    u_DUIO(i*96+1:i*96+32) = 0;
    u_DUIO(i*96+77:i*96+96) = 0;  
end

%%
figure()
hold on; grid on;
plot( u11, 'k','LineWidth',1)
plot( u1, 'Color',[0 0.4470 0.7410],'LineWidth',3.5)
plot( u3, 'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
plot( u4, 'Color',[0.4660 0.6740 0.1880],'LineWidth',1.5)
plot( u_UIO, 'LineWidth',1)
plot( u_DUIO, 'y','LineWidth',1)
a = sum(abs(u11-u1))/220;
b = sum(abs(u11'-u3))/220;
c = sum(abs(u11'-u4))/220;
d = sum(abs(u11-u_UIO))/220;
e = sum(abs(u11-u_DUIO))/220;
legend({'real occupancy',sprintf('Proposed methods, MAE: %0.3f',a), sprintf('Liear, MAE: %0.3f',b),...
    sprintf('GPR, MAE: %0.3f',c),sprintf('UIO, MAE: %0.3f',d),sprintf('DD-UIO, MAE: %0.3f',e)},'FontSize',18) %,sprintf('SVR, MAE: %0.3f',d)

xlabel('Step: 1 step = 15 minutes', 'interpreter', 'latex','fontsize',20);
ylabel('Occupancy', 'interpreter', 'latex','fontsize',20);
title('Occupancy estimation','fontsize',20, 'interpreter', 'latex')

%%
ddd = table;
% ddd.tt = vec(cache.cl.time(1:end-1));
ddd.t = vec([1:480]);

ddd.u = vec(u11(:));
ddd.u_Luen = vec(u1(:));
ddd.u_LR = vec(u3(:));
ddd.u_GR = vec(u4(:));
ddd.u_UIO = vec(u_UIO(:));
ddd.u_DUIO = vec(u_DUIO(:));
writetable(ddd, 'exp_LCSS.dat', 'Delimiter',',');
%% First submit version
for i = 0:4
    u2(i*96+1:i*96+28) = 0;
    u2(i*96+93:i*96+96) = 0;
    u3(i*96+1:i*96+28) = 0;
    u3(i*96+93:i*96+96) = 0;
    u4(i*96+1:i*96+28) = 0;
    u4(i*96+93:i*96+96) = 0;   
end

%% 8am-7pm
for i = 0:4
    u1(i*96+1:i*96+32) = 0;
    u1(i*96+77:i*96+96) = 0;   
    u2(i*96+1:i*96+32) = 0;
    u2(i*96+77:i*96+96) = 0;
    u3(i*96+1:i*96+32) = 0;
    u3(i*96+77:i*96+96) = 0;
    u4(i*96+1:i*96+32) = 0;
    u4(i*96+77:i*96+96) = 0;   
end
%%
mean(abs(u11 - u1))
% mean(abs(u11 - u2))
mean(abs(u11'- u3))
mean(abs(u11'-u4))
%%
figure()
hold on; grid on;
plot( u11, 'k','LineWidth',1)
plot( u1, 'Color',[0 0.4470 0.7410],'LineWidth',3.5)
plot( u3, 'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
plot( u4, 'Color',[0.4660 0.6740 0.1880],'LineWidth',1.5)
plot( u2, 'y','LineWidth',1)
a = sum(abs(u11-u1))/220
b = sum(abs(u11'-u3))/220
c = sum(abs(u11'-u4))/220
d = sum(abs(u11-u2))/220
legend({'real occupancy',sprintf('Proposed methods, MAE: %0.3f',a), sprintf('Liear, MAE: %0.3f',b),...
    sprintf('GPR, MAE: %0.3f',c),sprintf('SVR, MAE: %0.3f',d)},'FontSize',18) %,sprintf('SVR, MAE: %0.3f',d)

xlabel('Step: 1 step = 15 minutes', 'interpreter', 'latex','fontsize',20);
ylabel('Occupancy', 'interpreter', 'latex','fontsize',20);
title('Occupancy estimation','fontsize',20, 'interpreter', 'latex')

%%
ddd = table;
% ddd.tt = vec(cache.cl.time(1:end-1));

ddd.t = vec([1:480]);

ddd.u = vec(u11(:));
ddd.u_Luen = vec(u1(:));
ddd.u_LR = vec(u3(:));
ddd.u_GR = vec(u4(:));
writetable(ddd, 'exp.dat', 'Delimiter',',');

