clear;clc;close all;
format long

%% Data Loading
kf_example02a

% Modifications for the exam
Qk     =    6.00000000000000;                         % for all k
Rk     =    0.05000000000000;                         % for all k

%% Initialization
Nsamp = 1000;
kmax = size(thist,1);
[nx,nv] = size(Gammak);
nz = size(Hk,1);
x_hat_vec = zeros(nx,(kmax+1));
x_hat_vec(:,1) =  xhat0;
P_xx_vec = zeros(nx,nx,(kmax+1));
P_xx_vec(:,:,1) = P0;

x1_hat_vec = zeros(Nsamp,(kmax+1));
x2_hat_vec = zeros(Nsamp,(kmax+1));

x1_vec = zeros(Nsamp,(kmax+1));
x2_vec = zeros(Nsamp,(kmax+1));
W_vec = zeros(nx,nz,kmax);

P11_KF = zeros(1,(kmax+1));
P22_KF = zeros(1,(kmax+1));
P12_KF = zeros(1,(kmax+1));

P11_KF(1) = P_xx_vec(1,1,1);
P22_KF(1) = P_xx_vec(2,2,1);
P12_KF(1) = P_xx_vec(1,2,1);

%% Propogation of Cov and KF Gain Computation

for k = 1:kmax
    
    % Propagation state-error cov
    P_bar_xx_kp1 = Fk*P_xx_vec(:,:,k)*Fk'+Gammak*Qk*Gammak';
    
    % Gain Computation state-error cov
    S_kp1 = Hk*P_bar_xx_kp1*Hk'+Rk;
    W_vec(:,:,k) = P_bar_xx_kp1*Hk'/S_kp1;
    temp_joseph = eye(nx) - W_vec(:,:,k)*Hk;
    
    % Update state-error cov
%     P_xx_vec(:,:,k+1) = P_bar_xx_kp1 - W_vec(:,:,k)*S_kp1*W_vec(:,:,k)';
    P_xx_vec(:,:,k+1) = temp_joseph*P_bar_xx_kp1*temp_joseph'...
                        + W_vec(:,:,k)*Rk*W_vec(:,:,k)';
    
    P11_KF(k+1) = P_xx_vec(1,1,k+1);
    P22_KF(k+1) = P_xx_vec(2,2,k+1);
    P12_KF(k+1) = P_xx_vec(1,2,k+1);
    
end

%% Truth-model Monte-Carlo and State estimation

for n = 1:Nsamp
    
    [xtruehist,zhist] = kf_truthmodel(Fk,Gammak,Hk,Qk,Rk,xhat0,...
                                           P0,kmax);
for k = 1:kmax
    
    % Propagation of state and state-error cov
    x_bar_kp1 = Fk*x_hat_vec(:,k);
    
    % Measurement update of state and state-error cov
    z_kp1 = zhist(k);
    W_kp1 = W_vec(:,:,k);
    Innov_kp1 = z_kp1 - Hk*x_bar_kp1;
    x_hat_vec(:,k+1) = x_bar_kp1 + W_kp1*Innov_kp1;
  
end
    x1_hat_vec(n,:) = x_hat_vec(1,:);
    x2_hat_vec(n,:) = x_hat_vec(2,:);
    
    x1_vec(n,:) = xtruehist(:,1)';
    x2_vec(n,:) = xtruehist(:,2)';
    
end

x1_tilde_vec = x1_vec-x1_hat_vec;
x2_tilde_vec = x2_vec-x2_hat_vec;

P11_samp = mean(x1_tilde_vec.^2);
P22_samp = mean(x2_tilde_vec.^2);
P12_samp = mean(x1_tilde_vec.*x2_tilde_vec);

%% output

disp(['Sample Mean of state estimation error at k=10, for a number of samples:' num2str(Nsamp)]);
disp([mean(x1_tilde_vec(:,11));mean(x2_tilde_vec(:,11))]);
disp(['Sample Mean of state estimation error at k=35, for a number of samples:' num2str(Nsamp)]);
disp([mean(x1_tilde_vec(:,36));mean(x2_tilde_vec(:,36))]);

Psamp_10 = [P11_samp(11) P12_samp(11);...
            P12_samp(11) P22_samp(11)];
        
Psamp_35 = [P11_samp(36) P12_samp(36);...
            P12_samp(36) P22_samp(36)];

disp(['KF-Predicted state estimation error Covariance at k=10, for a number of samples:' num2str(Nsamp)]);
disp(P_xx_vec(:,:,11));
disp(['Sample state estimation error Covariance at k=10, for a number of samples:' num2str(Nsamp)]);
disp(Psamp_10);


disp(['KF-Predicted state estimation error Covariance at k=35, for a number of samples:' num2str(Nsamp)]);
disp(P_xx_vec(:,:,36));
disp(['Sample state estimation error Covariance at k=35, for a number of samples:' num2str(Nsamp)]);
disp(Psamp_35);



%% Plots

t_vec = [0,thist'];

figure;
subplot(2,1,1)
ps = plot(t_vec,x1_vec);hold on;grid on;

    ax = gca;
    ax.LineWidth = 1;
    ax.GridColor = [0 0 0];
    ax.MinorGridColor = 'k';
    ylabel('x_1')
    legend([ps(1)],'Monte-Carlo')
subplot(2,1,2)
plot(t_vec,x2_vec);hold on;grid on;

    ax = gca;
    ax.LineWidth = 1;
    ax.GridColor = [0 0 0];
    ax.MinorGridColor = 'k';
    xlabel('time, t (s)')
    ylabel('x_2')

    figure;
subplot(2,1,1)
ps = plot(t_vec,x1_hat_vec);hold on;grid on;

    ax = gca;
    ax.LineWidth = 1;
    ax.GridColor = [0 0 0];
    ax.MinorGridColor = 'k';
    ylabel('x_1 hat ')
    legend([ps(1)],'Monte-Carlo')
subplot(2,1,2)
plot(t_vec,x2_hat_vec);hold on;grid on;

    ax = gca;
    ax.LineWidth = 1;
    ax.GridColor = [0 0 0];
    ax.MinorGridColor = 'k';
    xlabel('time, t (s)')
    ylabel('x_2 hat')

    figure;
subplot(2,1,1)
ps = plot(t_vec,x1_tilde_vec);hold on;grid on;
pm = plot(t_vec,mean(x1_tilde_vec),'k','LineWidth',1.5);

    ax = gca;
    ax.LineWidth = 1;
    ax.GridColor = [0 0 0];
    ax.MinorGridColor = 'k';
    ylabel('x_1 tilde')
    legend([ps(1) pm],'Monte-Carlo','sample Mean')
subplot(2,1,2)
plot(t_vec,x2_tilde_vec);hold on;grid on;
pm = plot(t_vec,mean(x2_tilde_vec),'k','LineWidth',1.5);
    ax = gca;
    ax.LineWidth = 1;
    ax.GridColor = [0 0 0];
    ax.MinorGridColor = 'k';
    xlabel('time, t (s)')
    ylabel('x_2 tilde')
    
%% Sample covariances

    figure;
subplot(3,1,1)
ps = plot(t_vec,P11_samp);hold on;grid on;
pm = plot(t_vec,P11_KF,'k','LineWidth',1.5);

    ax = gca;
    ax.LineWidth = 1;
    ax.GridColor = [0 0 0];
    ax.MinorGridColor = 'k';
    ylabel('P11')
    legend([ps(1) pm],'sample Cov','KF Predicted')
subplot(3,1,2)
plot(t_vec,P22_samp);hold on;grid on;
pm = plot(t_vec,P22_KF,'k','LineWidth',1.5);
    ax = gca;
    ax.LineWidth = 1;
    ax.GridColor = [0 0 0];
    ax.MinorGridColor = 'k';
    xlabel('time, t (s)')
    ylabel('P22')


subplot(3,1,3)
plot(t_vec,P12_samp);hold on;grid on;
pm = plot(t_vec,P12_KF,'k','LineWidth',1.5);
    ax = gca;
    ax.LineWidth = 1;
    ax.GridColor = [0 0 0];
    ax.MinorGridColor = 'k';
    xlabel('time, t (s)')
    ylabel('P12')






