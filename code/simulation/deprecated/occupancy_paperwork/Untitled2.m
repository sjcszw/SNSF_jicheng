%% Luenberger observer and Kalman filter for state estimation
ifExp=1;

if ifExp
    % Kalman filter
    W = model.K*model.K'*model.NoiseVariance; 
    V = NoiseVariance;
    M = model.K;
    x_est(:,1) = zeros(nx,1);
    P_post = W;    
    for i = 1:T_tot-1
 		% Time Update
        x_pri	= A * x_est(:,i) + B * u_cl(:, i);
        P_pri	= A * P_post * A' + W;        
        if ifIno
            % Measurements Update
            k_post	= (P_pri * C' + M) * (C * P_pri * C' + V + C*M + M'*C' )^-1;
            x_est(:, i+1)	= x_pri + k_post * (y - C * x_pri);
            P_post	= (eye(nx) - k_post*C) * P_pri - k_post*M;           
        else
            % Measurements Update
            k_post	= P_pri * C'* (C * P_pri * C' + V)^-1;
            x_est(:, i+1)	= x_pri + k_post * (y - C * x_pri);
            P_post	= (eye(nx) - k_post*C) * P_pri;
        end
%         L_kal = A*qvar_kal*C'*(C*qvar_kal*C' + V)^(-1);
%         x_est(:,i+1) = A*x_est(:,i) + B*u_cl(:,i) + L_kal*(y_cl(:,i)-C*x_est(:,i));
%         qvar_kal = qvar_kal - qvar_kal*C'*(C*qvar_kal*C' + V)^(-1)*C*qvar_kal;
%         qvar_kal = A*qvar_kal*A' + W;
    end   
    figure()
    hold on
    plot(C*x_est(:,50:T_tot) - y_cl(:,50:T_tot), 'r')
    mean(abs(C*x_est(:,50:T_tot) - y_cl(:,50:T_tot)))

else
    Q = eye(nx); 
    R = 1;
    [K,~,~] = dlqr(A',C',Q,R);
    K = -K;                             
    eig(A'+ C'*K)
    % Luenberger observer
    clear x_elec_est y_pred_elec
    x_est(:,1) = zeros(nx,1);

    for i = 1:T_tot-1
        x_est(:,i+1) = A*x_est(:,i) + B*u_cl(:,i) - K'*(y_cl(:,i)-C*x_est(:,i));
    end
    figure()
    hold on
    plot(C*x_est(:,50:T_tot) - y_cl(:,50:T_tot), 'r')
    mean(abs(C*x_est(:,50:T_tot) - y_cl(:,50:T_tot)))
    
end


