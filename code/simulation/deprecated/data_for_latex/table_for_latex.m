 % u1: D~=0; u2: D = 0 
 % p1: LMI u
 % p2: Luen u
 % r L real u
% check
% u1_r1-u1_r2
% u2_r1-u2_r2
% difference
du1_LMI = u1_r1 - u1_p1;
du1_Luen = u1_r2 - u1_p2;
du2_LMI = u2_r1 - u2_p1;
du2_Luen = u2_r2 - u2_p2;

%%
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
