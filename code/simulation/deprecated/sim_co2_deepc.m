% polydome_hankel_prediction_mystruct
clc
clear

%%
N_ini = 10;
N_pred = 12;
% 136 - 1810
% start = 1; last = 700;
start = 401; last = 1111;
T = last-start+1; % -700-1+1; % T >= (nu + 1)*(N_ini + N_pred + nx) - 1
[exp_elec_iden, exp_energy_iden] = create2data(exp,start,last);
% start = 701; last = 1111;
start = 1; last = 400;
[exp_elec_val, exp_energy_val] = create2data(exp,start,last);
% % % detrend
% exp_elec_iden.y = exp_elec_iden.y - mean(exp_elec_iden.y,2);
% exp_energy_iden.y = exp_energy_iden.y - mean(exp_energy_iden.y,2);
% exp_elec_val.y = exp_elec_val.y - mean(exp_elec_val.y,2);
% exp_energy_val.y = exp_energy_val.y - mean(exp_energy_val.y,2);
ny = 1; nu = 1; nw = 3; nx = 3;
time_val = exp.time(start:last);
co2_val = exp.co2(start:last,1); people_val = exp.people(start:last,1);
%%
% time = exp.time + 1/24;
% h=figure;
% hold on
% set(h,'Units','normalized','Position',[0 0 1 .5]); 
% yyaxis left
% plot(exp_elec_iden.u,'b','LineWidth',1)
% yyaxis right
% ylim([5 30])
% %     plot(tExp,exp.setpoint_winter.value/10.0,'r','LineWidth',1)
% plot(exp_elec_iden.w(1,:),'r','LineWidth',1)
% plot(exp_elec_iden.y(1,:),'g','LineWidth',1)
% % plot(time,exp.supply_temp/10.0,'y','LineWidth',1.5)
% title('Meteorological data','FontSize',18)
% legend({'power','air temperature','room temperature','supply temperature'},'FontSize',18)
% ylabel('[°C]','FontSize',18)

%%

% Define the weights for g vector: L2 norm or adaptive term
% diag_vector = 0.02:- 0.018/(T-N_ini-N_pred):0.002;
% diag_vector = 2:- 0.1/(T-N_ini-N_pred):1.9;
diag_vector = 0.02:- 0.001/(T-N_ini-N_pred):0.019;
Qg = diag(diag_vector);

% build the Henkal matrix 
disp('Computing Hankel...')
N_total = N_ini + N_pred;
Hankel_col = T - N_total + 1;
Hankel_U_elec = compute_Hankel_matrix(exp_elec_iden.u,nu,N_total, Hankel_col);
Hankel_W_elec = compute_Hankel_matrix(exp_elec_iden.w,nw,N_total, Hankel_col);
Hankel_U_energy = compute_Hankel_matrix(exp_energy_iden.u,nu,N_total, Hankel_col);
Hankel_W_energy = compute_Hankel_matrix(exp_energy_iden.w,nw,N_total, Hankel_col);

N_total = N_ini + N_pred + 1;
Hankel_col = T - N_total + 2;
Hankel_Y_elec = compute_Hankel_matrix(exp_elec_iden.y,ny,N_total, Hankel_col);
Hankel_Y_energy = compute_Hankel_matrix(exp_energy_iden.y,ny,N_total, Hankel_col);

%%
N_total = N_ini + N_pred + nx;
Hankel_col = T - N_total +1 ;

% Need: u is persistently exciting of order is N_ini + N_pred + nx
Hankel_UW_check = compute_Hankel_matrix([exp_elec_iden.u; exp_elec_iden.w], nu+nw, N_total, Hankel_col);


% Check if the past u trajectory is persistently exciting of order N_ini + N_pred + nx
if rank(Hankel_UW_check)== (nu+nw)*(N_ini + N_pred + nx) 
    disp('Hankel rank is ok')
else
    disp('Exciting of order of Hu is samller than N_ini + N_pred + nx')
end
%% prediction
% y_pred = Hg_pred*[y_ini;u_ini;u_pred;w_ini;w_pred];
Hg_pred_elec = Hg_problem(Hankel_U_elec,Hankel_W_elec,Hankel_Y_elec,ny,N_ini,Qg);
Hg_pred_energy = Hg_problem(Hankel_U_energy,Hankel_W_energy,Hankel_Y_energy,ny,N_ini,Qg);

%
num = size(exp_elec_val.u,2);
t = 49;
while(t<=num-1)
    if t+N_pred<=num
        uip_elec = reshape(exp_elec_val.u(:,t-N_ini:t+N_pred-1), [(N_ini+N_pred)*nu,1]);
        wip_elec = reshape(exp_elec_val.w(:,t-N_ini:t+N_pred-1), [(N_ini+N_pred)*nw,1]);
        yi_elec = reshape(exp_elec_val.y(:,t-N_ini:t), [(N_ini+1)*ny,1]); 
        y_pred_elec(t+1:t+N_pred) = Hg_pred_elec*[yi_elec;uip_elec;wip_elec];

        uip_energy = reshape(exp_energy_val.u(:,t-N_ini:t+N_pred-1), [(N_ini+N_pred)*nu,1]);
        wip_energy = reshape(exp_energy_val.w(:,t-N_ini:t+N_pred-1), [(N_ini+N_pred)*nw,1]);
        yi_energy = reshape(exp_energy_val.y(:,t-N_ini:t), [(N_ini+1)*ny,1]); 
        y_pred_energy(t+1:t+N_pred) = Hg_pred_energy*[yi_energy;uip_energy;wip_energy];

        t = t+N_pred;
    else
        break;
    end
end
%% plot
figure()
hold on; grid on;
start = 50; last = size(y_pred_elec,2);
plot( exp_elec_val.y(:,start:last), 'Color',[0 0.4470 0.7410],'LineWidth',1)
plot( y_pred_elec(:,start:last), 'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
plot( y_pred_energy(:,start:last), 'Color',[0.4660 0.6740 0.1880],'LineWidth',1)
scatter( [(start:N_pred:last)-start+1],y_pred_energy(:,start:N_pred:last), 100,'red')
'elec'
mae_elec = mean(abs(exp_elec_val.y(:,start:last)-y_pred_elec(:,start:last)))
fit = 1 - goodnessOfFit(y_pred_elec(:,start:last)', exp_elec_val.y(:,start:last)', 'NRMSE')
'energy'
mae_energy = mean(abs(exp_energy_val.y(:,start:last)-y_pred_energy(:,start:last)))
fit = 1 - goodnessOfFit(y_pred_energy(:,start:last)', exp_energy_val.y(:,start:last)', 'NRMSE')
legend({'real temperature',sprintf('prediction with u: electricity, MAE: %0.3f °C',mae_elec),...
    sprintf('prediction with u: energy, MAE: %0.3f °C',mae_energy),'prediction point'},'FontSize',18)

xlabel('Step: 1 step = 15 minutes', 'interpreter', 'latex','fontsize',20);
ylabel('Temperature [$^{\circ}$C] ', 'interpreter', 'latex','fontsize',20);
title('Prediction by bi-level deepc: $N_{pred}$ = 12 = 3 hours, T = 5 days','fontsize',20, 'interpreter', 'latex')
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
%% compute the prediction matrix
function Hg_pred = Hg_problem(Hankel_U,Hankel_W,Hankel_Y,ny,N_ini,Qg)
    % y_pred = Hg_pred*[y_ini;u_ini;u_pred;w_ini;w_pred];
%             Hu_init = Hankel_U(1:nu*N_ini, :);
%             Hw_init = Hankel_W(1:nw*N_ini, :);
            Hy_init = Hankel_Y(1:ny*(N_ini+1), :);
%             Hu_pred = Hankel_U(nu*N_ini+1:end, :);
%             Hw_pred = Hankel_W(nw*N_ini+1:end, :);
            Hy_pred = Hankel_Y(ny*(N_ini+1)+1:end, :);           
%             H = [Hu_init;Hu_pred;Hw_init;Hw_pred];
            H = [Hankel_U;Hankel_W];
            
            % ====== Compute g_norm from KKT condition
            temp = inv(Qg);
            % desired inversion
            temp = temp - temp*Hy_init'*((eye((N_ini+1)*ny)... %/(2*w_s)...
                +Hy_init*temp*Hy_init')\Hy_init)*temp;
            temp_inv = temp*H'/(H*temp*H');
         
            Hg_pred = Hy_pred*[(temp-temp*H'*((H*temp*H')\H)*temp)*Hy_init', temp_inv];
%             y_pred = Hg_pred*[obj.y_ini;obj.u_ini;obj.u_pred;obj.w_ini;obj.w_norm];
%             Hg_init = Hy_init*[(temp-temp*H'*((H*temp*H')\H)*temp)*Hy_init', temp_inv];
%             slack = Hg_init*[obj.y_ini;obj.u_ini;obj.u_pred;obj.w_ini;obj.w_norm]-obj.y_ini;
end