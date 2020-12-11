clear;clc;close all;
format longg

%% Data Loading

% kf_example03a
kf_example03b

%% Initialization

kmax = size(thist,1);
[nx,nv] = size(Gammak);
nz = size(Hk,1);
x_hat_vec = zeros(nx,(kmax+1));
x_hat_vec(:,1) =  xhat0;
P_xx_vec = zeros(nx,nx,(kmax+1));
P_xx_vec(:,:,1) = P0;
W_vec = zeros(nx,nz,kmax);

P11_KF = zeros(1,(kmax+1));
P22_KF = zeros(1,(kmax+1));
P12_KF = zeros(1,(kmax+1));
P33_KF = zeros(1,(kmax+1));
P31_KF = zeros(1,(kmax+1));
P32_KF = zeros(1,(kmax+1));

P11_KF(1) = P_xx_vec(1,1,1);
P22_KF(1) = P_xx_vec(2,2,1);
P12_KF(1) = P_xx_vec(1,2,1);
P33_KF(1) = P_xx_vec(3,3,1);
P31_KF(1) = P_xx_vec(3,1,1);
P32_KF(1) = P_xx_vec(3,2,1);

P11_sirf = zeros(1,(kmax+1));
P22_sirf = zeros(1,(kmax+1));
P12_sirf = zeros(1,(kmax+1));
P33_sirf = zeros(1,(kmax+1));
P31_sirf = zeros(1,(kmax+1));
P32_sirf = zeros(1,(kmax+1));

P11_sirf(1) = P_xx_vec(1,1,1);
P22_sirf(1) = P_xx_vec(2,2,1);
P12_sirf(1) = P_xx_vec(1,2,1);
P33_sirf(1) = P_xx_vec(3,3,1);
P31_sirf(1) = P_xx_vec(3,1,1);
P32_sirf(1) = P_xx_vec(3,2,1);


%% KF State estimation and Cov Prop

for k = 1:kmax
    
    % Propagation of state and state-error cov
    x_bar_kp1 = Fk*x_hat_vec(:,k);
    P_bar_xx_kp1 = Fk*P_xx_vec(:,:,k)*Fk'+Gammak*Qk*Gammak';
    
    % Measurement update of state and state-error cov
    S_kp1 = Hk*P_bar_xx_kp1*Hk'+Rk;
    W_vec(:,:,k) = P_bar_xx_kp1*Hk'/S_kp1;
    z_kp1 = zhist(k);
    W_kp1 = W_vec(:,:,k);
    Innov_kp1 = z_kp1 - Hk*x_bar_kp1;
    temp_joseph = eye(nx) - W_kp1*Hk;
    
    x_hat_vec(:,k+1) = x_bar_kp1 + W_kp1*Innov_kp1;
%     P_xx_vec(:,:,k+1) = P_bar_xx_kp1 - W_vec(:,:,k)*S_kp1*W_vec(:,:,k)';
    P_xx_vec(:,:,k+1) = temp_joseph*P_bar_xx_kp1*temp_joseph'...
        + W_kp1*Rk*W_kp1';
    
    P11_KF(k+1) = P_xx_vec(1,1,k+1);
    P22_KF(k+1) = P_xx_vec(2,2,k+1);
    P12_KF(k+1) = P_xx_vec(1,2,k+1);
    P33_KF(k+1) = P_xx_vec(3,3,k+1);
    P31_KF(k+1) = P_xx_vec(3,1,k+1);
    P32_KF(k+1) = P_xx_vec(3,2,k+1);
  
end

%% SIRF Computations

x_hat_vec_srif = zeros(nx,kmax+1);
x_hat_vec_srif(:,1) =  xhat0;
P_xx_vec_srif  = zeros(nx,nx,(kmax+1));
P_xx_vec_srif(:,:,1) = P0;

z_x_vec = zeros(nx,kmax+1);
R_xx_vec = zeros(nx,nx,kmax+1);
I_scr = eye(nx)/P0;
R_xx_vec(:,:,1) = chol(I_scr);
z_x_vec(:,1) = R_xx_vec(:,:,1)*xhat0;

R_vv = (eye(nv)/chol(Qk))'; 

F_inv = eye(nx)/Fk;
R_a = chol(Rk);
R_a_invTr = (eye(nz)/R_a)';

H_a = R_a_invTr*Hk;

for k = 1:kmax
    
    % Compute R_xx_bar_kp1, R_vv_bar_k, R_vx_bar_kp1 and T_a_k. Dynamic
    % Propagation
    R_xx_mult_F_inv = R_xx_vec(:,:,k)*F_inv;
    Block_Mat_for_Ta = [R_vv zeros(nv,nx) ;...
           -R_xx_mult_F_inv*Gammak R_xx_mult_F_inv];
    [T_a_Tr,R_block] = qr(Block_Mat_for_Ta);
    T_a_k = T_a_Tr';
    
    R_vv_bar_k = R_block(1:nv,1:nv);
    R_vx_bar_kp1 = R_block(1:nv,nv+1:end);
    R_xx_bar_kp1 = R_block(nv+1:end,nv+1:end);
    
    % Compute z_bar_vc. Dynamic Propagation
    z_bar_vec = T_a_k*[zeros(nv,1);...
                    z_x_vec(:,k)];
    z_x_bar_kp1 = z_bar_vec(nv+1:end);
    
    % Compute R_xx_kp1 and T_b_kp1. Measurment update
    Block_Mat_for_Tb = [R_xx_bar_kp1 ;...
                         H_a];
    [T_b_Tr,R_block2] = qr(Block_Mat_for_Tb);
    T_b_kp1 = T_b_Tr';
    R_xx_vec(:,:,k+1) = R_block2(1:nx,:);
    
    % Compute z_kp1_vec. Measurment update
    z_vec_kp1 = T_b_kp1*[z_x_bar_kp1 ;...
                    R_a_invTr*zhist(k)];
    z_x_vec(:,k+1) = z_vec_kp1(1:nx);

    R_xx_inv_kp1 = eye(nx)/R_xx_vec(:,:,k+1);
    x_hat_vec_srif(:,k+1) = R_xx_inv_kp1*z_x_vec(:,k+1);
    P_xx_vec_srif(:,:,k+1) = R_xx_inv_kp1*R_xx_inv_kp1';
    
    P11_sirf(k+1) = P_xx_vec_srif(1,1,k+1);
    P22_sirf(k+1) = P_xx_vec_srif(2,2,k+1);
    P12_sirf(k+1) = P_xx_vec_srif(1,2,k+1);
    P33_sirf(k+1) = P_xx_vec_srif(3,3,k+1);
    P31_sirf(k+1) = P_xx_vec_srif(3,1,k+1);
    P32_sirf(k+1) = P_xx_vec_srif(3,2,k+1);
    
end

%% Output

disp('Error in Final Estimated State:')
disp(x_hat_vec(:,end)-x_hat_vec_srif(:,end));

disp('Error in Final State-estimation Error Covariance:')
disp(P_xx_vec(:,:,end)-P_xx_vec_srif(:,:,end));

disp('Matrix Error Metric:')
MEM = abs(P_xx_vec(:,:,end)-P_xx_vec_srif(:,:,end))./(abs(P_xx_vec_srif(:,:,end))+eps^6);
disp(MEM)


