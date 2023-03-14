%%
for i = 0:4
    u2(i*96+1:i*96+28) = 0;
    u2(i*96+93:i*96+96) = 0;
    u3(i*96+1:i*96+28) = 0;
    u3(i*96+93:i*96+96) = 0;
    u4(i*96+1:i*96+28) = 0;
    u4(i*96+93:i*96+96) = 0;   
end

%% 7am-7pm
for i = 0:4
    u1(i*96+1:i*96+28) = 0;
    u1(i*96+77:i*96+96) = 0;   
    u2(i*96+1:i*96+28) = 0;
    u2(i*96+77:i*96+96) = 0;
    u3(i*96+1:i*96+28) = 0;
    u3(i*96+77:i*96+96) = 0;
    u4(i*96+1:i*96+28) = 0;
    u4(i*96+77:i*96+96) = 0;   
end
%%
mean(abs(u11-u1))
mean(abs(u11'-u2))
mean(abs(u11'-u3))
mean(abs(u11'-u4))
%%
figure()
hold on; grid on;
plot( u11, 'k','LineWidth',1)
plot( u1, 'Color',[0 0.4470 0.7410],'LineWidth',3.5)
plot( u3, 'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
plot( u4, 'Color',[0.4660 0.6740 0.1880],'LineWidth',1.5)
plot( u2, 'y','LineWidth',1)
a = mean(abs(u11-u1))
b = mean(abs(u11'-u3))
c = mean(abs(u11'-u4))
d = mean(abs(u11'-u2))
legend({'real occupancy',sprintf('Proposed methods, MAE: %0.3f',a), sprintf('Liear, MAE: %0.3f',b),...
    sprintf('GPR, MAE: %0.3f',c),sprintf('SVR, MAE: %0.3f',d)},'FontSize',18) %,sprintf('SVR, MAE: %0.3f',d)

xlabel('Step: 1 step = 15 minutes', 'interpreter', 'latex','fontsize',20);
ylabel('Occupancy', 'interpreter', 'latex','fontsize',20);
title('Occupancy estimation','fontsize',20, 'interpreter', 'latex')