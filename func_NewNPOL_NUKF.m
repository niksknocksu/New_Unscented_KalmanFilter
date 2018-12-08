function func_NewNPOL_NUKF(sheet,MC)
% Use excel file NPOL_new.xlsx to generate observer data and target data
file_in = 'NPOL_new.xlsx';
%sheet = 12 ;
n = 4;  % Dimension on state vector
A = xlsread(file_in,sheet,'');
T =1;
%kmax = numel(A.data.Sheet1(:,1));
x_obs(1,:) = A(:,1)';
x_obs(2,:) = A(:,2)';
x_obs(3,:) = A(:,3)';
x_obs(4,:) = A(:,4)';

x_tar(1,:) = A(:,5)';
x_tar(2,:) = A(:,6)';
x_tar(3,:) = A(:,7)';
x_tar(4,:) = A(:,8)';

F =[1 0 T 0;...
    0 1 0 T;...
    0 0 1 0;...
    0 0 0 1]; % State Transition matrix
q=[T^3/3 0 T^2/2 0;...
   0 T^3/3 0 T^2/2;...
   T^2/2 0 T 0;...
   0 T^2/2 0 T];
Q=((1.944E-6)/216000)*q; %Process noise covariance

kmax = numel(x_obs(1,:));
sigma_theta = 2*pi/180;
R = sigma_theta^2;
%MC = 1000;
track_loss = 0;
 for mm = 1:MC
    fprintf('Run %d/%d: ' ,mm,MC);
    % Relative state vector
    x_rel(:,1) = x_tar(:,1)-x_obs(:,1);
    
   for k = 2:kmax
    x_rel(:,k) = (F*x_rel(:,k-1)) + mvnrnd([0,0,0,0],Q)' - [x_obs(1,k) - x_obs(1,k-1) - (T*x_obs(3,k-1));...
                                                x_obs(2,k) - x_obs(2,k-1) - (T*x_obs(4,k-1));...
                                                x_obs(3,k) - x_obs(3,k-1);...
                                                x_obs(4,k) - x_obs(4,k-1)];
   end
   
   x_tar = x_rel + x_obs;
   init_r = norm([x_rel(1,1),x_rel(2,1)]);
   tar_speed = norm([x_tar(3,1),x_tar(4,1)]);
   
    % Measurements
    z = zeros(1,kmax);
    for i = 1:kmax
        z(i) = atan2(x_rel(1,i),x_rel(2,i)) + mvnrnd(0,R);
    end
    
    % Initialization 
    x_rel_obs = zeros(4,kmax);
    sigma_r = 1;
    init_r = init_r + sigma_r*randn(1);
    init_sp = tar_speed;
    sigma_sp = 1*0.0308667/60; %kms/sec
    init_sp = init_sp + sigma_sp*randn(1);
    init_co = z(1);
    sigma_c = pi/sqrt(12);
    x_rel_obs(:,1) = [init_r*sin(z(1));...
                      init_r*cos(z(1));...
                      (init_sp*sin(init_co))-x_obs(3,1);...
                      (init_sp*cos(init_co))-x_obs(4,1);];
    P_11 = (init_r^2)*R*(cos(z(1))^2) + (sigma_r^2)*(sin(z(1))^2);
    P_22 = (init_r^2)*R*(sin(z(1))^2) + (sigma_r^2)*(cos(z(1))^2);
    P_12 = ((sigma_r^2) - ((init_r^2)*R))*sin(z(1))*cos(z(1));
    P_21 = P_12;
    P_33 = ((init_sp^2)*(sigma_c^2)*(cos(init_co)^2)) + (sigma_sp^2)*(sin(init_co)^2);
    P_44 = ((init_sp^2)*(sigma_c^2)*(sin(init_co)^2)) + (sigma_sp^2)*(cos(init_co)^2);
    P_34 = ((sigma_sp^2) - (init_sp^2)*(sigma_c^2))*sin(init_co)*cos(init_co);
    P_43 = P_34;
    
    
    P_xx = [P_11,P_12,0,0;...
            P_21,P_22,0,0;...
            0,0,P_33,P_34;...
            0,0,P_43,P_44;];
     
   error_pos(mm,1) = ((x_tar(1,1) - (x_rel_obs(1,1) +  x_obs(1,1)))^2) + ((x_tar(2,1) - (x_rel_obs(2,1) +  x_obs(2,1)))^2);
   
   error_vel(mm,1) = ((x_tar(3,1) - (x_rel_obs(3,1) +  x_obs(3,1)))^2) + ((x_tar(4,1) - (x_rel_obs(4,1) +  x_obs(4,1)))^2);
   
   %**********************************************************************
                % NUKF Algortithm starts here :
   %**********************************************************************
   m = 0.5;
   for t = 2:kmax
        
    x_rel_obs(:,t) = F*x_rel_obs(:,t-1)  - [x_obs(1,t) - x_obs(1,t-1) - (T*x_obs(3,t-1));...
                                            x_obs(2,t) - x_obs(2,t-1) - (T*x_obs(4,t-1));...
                                            x_obs(3,t) - x_obs(3,t-1);...
                                            x_obs(4,t) - x_obs(4,t-1)];
   P_xx = (F*P_xx*F') + Q;
   
   % Create alpha
    
    for i = 1:n
            alpha(i) = (abs(dot(x_rel_obs(:,t), P_xx(:,i))))/(norm(x_rel_obs(:,t))*norm(P_xx(:,i)));
    end
    
    beta_n = ((0.25*max(m.*alpha)) - (0.5*sum(alpha))) + 1;
    shi = sum(alpha) + beta_n;
    
    % Create Wieghts
         w = zeros(1,4*n+1);
         w(1) = 1 - (sum(alpha)*0.5/shi);
        for i = 2:4*n+1
            if i <= n+1
                w(i) = m*alpha(i-1)*0.25/shi;
            elseif i>n+1 && i<=2*n+1
                w(i) = m*alpha(i-(n+1))*0.25/shi;
            elseif i>2*n+1 && i<=3*n+1
                w(i) = (1-m)*alpha(i-(1+2*n))*0.25/shi;
            elseif i>3*n+1 && i<=4*n+1
                w(i) = (1-m)*alpha(i-(3*n+1))*0.25/shi;
            end
        end
        
     % Generate Sigma Points
      
        U = chol(P_xx , 'lower');
        x_p = zeros(n,4*n+1);
        x_p(:,1) = x_rel_obs(:,t);
        
        for i = 2:4*n+1
            if i<=n+1
                x_p(:,i) = x_rel_obs(:,t) + (sqrt(shi/(m*alpha(i-1)))*U(:,i-1));
            elseif i>n+1 && i<=2*n+1
                x_p(:,i) = x_rel_obs(:,t) - (sqrt(shi/(m*alpha(i-(n+1))))*U(:,i-(n+1)));
            elseif i>2*n+1 && i<=3*n+1
                x_p(:,i) = x_rel_obs(:,t) + (sqrt(shi/((1-m)*alpha(i-(2*n+1))))*U(:,i-(2*n+1)));
            elseif i>3*n+1 && i<=4*n+1
                x_p(:,i) = x_rel_obs(:,t) - (sqrt(shi/((1-m)*alpha(i-(3*n+1))))*U(:,i-(3*n+1)));
            end
        end
        
        %Pass sigma points through measurement functions and find gamma_m
      
        gamma_m = zeros(1,4*n+1);
        
        for i = 1:4*n+1
            gamma_m(:,i) = atan2(x_p(1,i),x_p(2,i));
        end
        
        % Predicted measurement
        z_pred = 0;
        for i = 1:4*n+1
        z_pred = z_pred + (w(i)*gamma_m(:,i));
        end
        
         %Predicted Measurement covariance
        P_zz = 0;
        for i = 1:4*n+1
            P_zz = P_zz + (w(i).*((gamma_m(:,i) - z_pred)*(gamma_m(:,i) - z_pred)'));
        end
        P_zz = P_zz + R;
        
        % Measurement-state cross covariance
        
        P_xz = zeros(4,1);
        for i = 1:4*n+1
            P_xz = P_xz + (w(i).*((x_p(:,i) - x_rel_obs(:,t))*((gamma_m(:,i) - z_pred)')));
        end
        
        K = P_xz/(P_zz); % Kalman Gain
        x_rel_obs(:,t) = x_rel_obs(:,t) + (K*(z(t)-z_pred)); % Posterior estimare
        P_xx = P_xx - (K*P_zz*K'); % Posterior state covariance
        
        error_pos(mm,t) = ((x_tar(1,t) - (x_rel_obs(1,t) +  x_obs(1,t)))^2) + ((x_tar(2,t) - (x_rel_obs(2,t) +  x_obs(2,t)))^2);
        error_vel(mm,t) = ((x_tar(3,t) - (x_rel_obs(3,t) +  x_obs(3,t)))^2) + ((x_tar(4,t) - (x_rel_obs(4,t) +  x_obs(4,t)))^2);
        
   end
   if sqrt(error_pos(mm,kmax)) > 1
        fprintf('filter diverged, erros > threshold \n')
        track_loss = track_loss +1;
        error_pos(mm,:) = zeros(1,kmax);
        error_vel(mm,:) = zeros(1,kmax);
    else 
        fprintf('passed \n')
    end
   
 end
 
 x_fil_org = x_rel_obs + x_obs;

error_pos = error_pos(any(error_pos,2),:);
error_vel = error_vel(any(error_vel,2),:);
RMSE_pos = ((1/(MC-track_loss)).*sum(error_pos,1)).^0.5;
RMSE_vel = ((1/(MC-track_loss)).*sum(error_vel,1)).^0.5;
 
fprintf('No of times track diverged is %d of %d MC runs\n', track_loss,MC);
fprintf('Track loss for %d Monte carlo runs:  %.4f %% \n', MC, track_loss*100/MC);

file_out = 'RMSE_results.xlsx';

xlswrite(file_out,{'RMSE Position in kilometers',},sheet,'A1');
xlswrite(file_out,RMSE_pos',sheet,'A2');

xlswrite(file_out,{'RMSE Velocity in kmps'},sheet,'B1');
xlswrite(file_out,RMSE_vel',sheet,'B2');


xlswrite(file_out,{'No. of MC runs'},sheet,'C1');
xlswrite(file_out,MC,sheet,'C2');


xlswrite(file_out,{'Track loss'},sheet, 'D1');
xlswrite(file_out,track_loss,sheet, 'D2');

xlswrite(file_out,{'Track loss %'},sheet, 'E1');
xlswrite(file_out,(track_loss*100/MC),sheet, 'E2');

figure
plot(x_obs(1,:),x_obs(2,:),'b'); % Observer Plot X vs Y
hold on
plot(x_tar(1,:),x_tar(2,:),'r'); % Target Plot X vs Y
hold on
plot(x_fil_org(1,:),x_fil_org(2,:),'g','linewidth',1.5) % Est Target Plot X vs Y
hold on
plot(x_obs(1,1),x_obs(2,1),'r*');
hold on
plot(x_fil_org(1,1),x_fil_org(2,1),'r*');
hold on
plot(x_tar(1,1),x_tar(2,1),'r*');
xlabel('Distance in kms','fontsize',14);
ylabel('Distance in kms','fontsize',14);
%legend('Observer','Target','estimated Target','initial observer position','initial estimated position');
title('Observer and Target ','fontsize',14);
   