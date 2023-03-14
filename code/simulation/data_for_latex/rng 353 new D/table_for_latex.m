 % u1: D~=0; u2: D = 0 
 % p1: LMI u
 % p2: Luen u
 % p3: UIO
 % p4: data-driven UIO
 % r L real u
% check

%  % u1: D~=0
u1_r1-u1_r2
u1_r1-u1_r3
u1_r1-u1_r4

 % u1: D=0
% u2_r1-u2_r2
% u2_r1-u2_r3
% u2_r1-u2_r4

%% difference
du1_LMI = u1_r1 - u1_p1;
du1_Luen = u1_r2 - u1_p2;
du2_LMI = u2_r1 - u2_p1;
du2_Luen = u2_r2 - u2_p2;

du1_UIO = u1_r3 - u1_p3;
du1_DUIO = u1_r4 - u1_p4;
du2_UIO = u2_r3 - u2_p3;
du2_DUIO = u2_r4 - u2_p4;


%% L-CSS version
ddd = table;
% ddd.tt = vec(cache.cl.time(1:end-1));

ddd.t = vec([-5:49]);
ddd.u1_r11 = vec(u1_r1(1,:));
ddd.u1_r12 = vec(u1_r1(2,:));
ddd.u1_p11 = vec(u1_p1(1,:));
ddd.u1_p12 = vec(u1_p1(2,:));
writetable(ddd, 'sim_D1_u1_LCSS.dat', 'Delimiter',',');
%%
ddd = table;
ddd.t = vec([0:49]);
ddd.du1_LMI1 = vec(du1_LMI(1,6:end));
ddd.du1_LMI2 = vec(du1_LMI(2,6:end));
ddd.du1_Luen1 = vec(du1_Luen(1,6:end));
ddd.du1_Luen2 = vec(du1_Luen(2,6:end));
ddd.du1_UIO1 = vec(du1_UIO(1,6:end));
ddd.du1_UIO2 = vec(du1_UIO(2,6:end));
ddd.du1_DUIO1 = vec(du1_DUIO(1,6:end));
ddd.du1_DUIO2 = vec(du1_DUIO(2,6:end));
writetable(ddd, 'sim_D1_LCSS.dat', 'Delimiter',',');
%%
ddd = table;
ddd.t = vec([0:48]);
ddd.du2_LMI1 = vec(du2_LMI(1,6:end));
ddd.du2_LMI2 = vec(du2_LMI(2,6:end));
ddd.du2_Luen1 = vec(du2_Luen(1,6:end));
ddd.du2_Luen2 = vec(du2_Luen(2,6:end));
ddd.du2_UIO1 = vec(du2_UIO(1,6:end));
ddd.du2_UIO2 = vec(du2_UIO(2,6:end));
ddd.du2_DUIO1 = vec(du2_DUIO(1,6:end));
ddd.du2_DUIO2 = vec(du2_DUIO(2,6:end));
writetable(ddd, 'sim_D0_LCSS.dat', 'Delimiter',',');

%% Plot it here
% D~=0
num = 2;
figure; hold on;

plot(du1_LMI(num,6:end))
plot(du1_Luen(num,6:end))
plot(du1_UIO(num,6:end))
plot(du1_DUIO(num,6:end))
% plot(u1_p1(num,:))
% plot(u1_p2(num,:))
% plot(u1_p3(num,:))
% plot(u1_p4(num,:))
%%
% D=0
num = 1;
figure; hold on;
plot(du2_LMI(num,6:end))
plot(du2_Luen(num,6:end))
plot(du2_UIO(num,6:end))
plot(du2_DUIO(num,6:end))

% plot(u2_p1(num,:))
% plot(u2_p2(num,:))
% plot(u2_p3(num,:))
% plot(u2_p4(num,:))
%% First submit version
ddd = table;
% ddd.tt = vec(cache.cl.time(1:end-1));

ddd.t = vec([-5:49]);
ddd.u1_r11 = vec(u1_r1(1,:));
ddd.u1_r12 = vec(u1_r1(2,:));
ddd.u1_p11 = vec(u1_p1(1,:));
ddd.u1_p12 = vec(u1_p1(2,:));
writetable(ddd, 'sim_D1_u1.dat', 'Delimiter',',');
%%
ddd = table;
ddd.t = vec([0:49]);
ddd.du1_LMI1 = vec(du1_LMI(1,6:end));
ddd.du1_LMI2 = vec(du1_LMI(2,6:end));
ddd.du1_Luen1 = vec(du1_Luen(1,6:end));
ddd.du1_Luen2 = vec(du1_Luen(2,6:end));
writetable(ddd, 'sim_D1.dat', 'Delimiter',',');
%%
ddd = table;
ddd.t = vec([0:48]);
ddd.du2_LMI1 = vec(du2_LMI(1,6:end));
ddd.du2_LMI2 = vec(du2_LMI(2,6:end));
ddd.du2_Luen1 = vec(du2_Luen(1,6:end));
ddd.du2_Luen2 = vec(du2_Luen(2,6:end));
writetable(ddd, 'sim_D0.dat', 'Delimiter',',');
%%
