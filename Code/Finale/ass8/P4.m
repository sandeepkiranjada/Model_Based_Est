%% Assignment 8 Problem 3

clc; clear; %close all;

%% Initialization

load('measdata_pfexample');

x_hat_0 = xhat0;
P_0 = P0;

nx = 1;
nz = 1;
N = 101;


x_hat = zeros(nx,N);
P = zeros(nx,N);

x_hat(1) = x_hat_0;
P(1) = P_0;

Gamma = 1



%% EKF

for kk = 1:N-1
    
    x_bar_kp1 = 2*atan(x_hat(kk)) + 0.5*cos(pi*kk/3);
    F_k = 2  /(1 + x_hat(kk)^2 );
    
    P_bar_kp1 = F_k^2*P(kk) + Q;
    
    z_bar_kp1 = x_bar_kp1 + x_bar_kp1^2 + x_bar_kp1^3;
    
    H_kp1 = 1 + 2*x_bar_kp1 + 3*x_bar_kp1^2;
    
    nu_kp1 = zkhist(kk) - z_bar_kp1;
    
    s_kp1 = H_kp1^2*P_bar_kp1+R;
    w_kp1 = P_bar_kp1*H_kp1/s_kp1;
    
    x_hat(kk+1) = x_bar_kp1 + w_kp1*nu_kp1;
    P(kk+1) = P_bar_kp1 - w_kp1^2*s_kp1;
    
end

figure(1);
plot([0:100],x_hat,'-b');grid on;hold on;
ylabel('x')
xlabel('sample index, k')
% ylim([-5 5])
    