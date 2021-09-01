%% Assignment 8 Problem 7

clc; clear; close all;

%% Initialization

load('measdata_pfexample02');

x_hat_0 = xhat0;
P_0 = P0;


nx = 1;
nz = 1;
N = 101;
Ns = 1e5;

x_hat = zeros(nx,N);
P = zeros(nx,N);

x_hat(1) = x_hat_0;
P(1) = P_0;

x_math = zeros(Ns,N+1);

x_math_0 = x_hat_0 + randn(Ns,1)*sqrt(P_0);

x_math_k = x_math_0;

z_hat_kp1 = zeros(Ns,N-1);
W_k = ones(Ns,1)/Ns;

%% Particle Filter

for kk = 1:N-1
    
    x_math_kp1 = 2*atan(x_math_k) + 0.5*cos(pi*kk/3) + randn(Ns,1)*sqrt(Q);
    
    z_hat_kp1(:,kk) = x_math_kp1 + x_math_kp1.^2 + x_math_kp1.^3;
    
    expo_temp = (zkhist(kk,1)*ones(Ns,1) - z_hat_kp1(:,kk)).^2/R/2;
    
    log_wk_m_expo = log(W_k) - expo_temp;
    
    [~,imax] = max(log_wk_m_expo);
    
    W_til_kp1 = exp(log_wk_m_expo - log_wk_m_expo(imax));
    
    W_kp1 = W_til_kp1./sum(W_til_kp1);
    
    x_hat(kk+1) = x_math_kp1'*W_kp1;
    
    
    del_x_kp1 = x_math_kp1-x_hat(kk+1);
    
    P(kk+1) = (del_x_kp1.^2)'*W_kp1;
    
    c = zeros(1,Ns+1);
    for ii = 2:Ns+1
        c(ii) = c(ii-1) + W_kp1(ii-1);
    end
    
    % Resamp, a version first approach
    
    [eta,indv] = sort(rand(1,Ns));
    
    nn = 1;
    x_temp = zeros(1,Ns);
    
    for ii = 1:Ns
        
        flag = 1;
        while flag && nn<=Ns
            
            if eta(nn) >= c(ii) && eta(nn) < c(ii+1)
                
                x_temp(indv(nn))=x_math_kp1(ii);
                nn = nn +1;
                flag = 1;
            else
                
                flag = 0;
            end
        end
        
    end
    
    x_math(:,kk) = (x_math_k);

    
    % reset
    x_math_k = x_temp';
    W_k = ones(Ns,1)/Ns;
    kk
end

x_math(:,N) = (x_math_k);

%% Smoothing

nu_sum = zeros(Ns,1);
for ii = 1:Ns
    nu_sum(ii,1) = sum(zkhist'-z_hat_kp1(ii,:));
end
expo_sum  = -0.5 * nu_sum.^2/R;

    [~,imax] = max(expo_sum);
    
    W_smooth = exp(expo_sum - expo_sum(imax));
    
    W_smooth = W_smooth./sum(W_smooth);
    
    x_star = zeros(nx,N);
    P_star = zeros(nx,N);

for kk = 1:N
    
    x_star(kk) = x_math(:,kk)'*W_smooth;
    
    
    del_x_kp1 = x_math(:,kk)-x_hat(kk);
    
    P_star(kk+1) = (del_x_kp1.^2)'*W_smooth;

end

%% Results

figure(1);
plot([0:100],x_hat,'--r');grid on;hold on;
plot([0:100],x_star,'--b');grid on;hold on;
ylabel('x')
xlabel('sample index, k')
legend('PF','Fixed-Int PS')



disp('Final State Estimate')
disp(x_hat(end))

disp('Final State Covariance')
disp(P(end))

disp('Final smoothed-state Estimate')
disp(x_hat(end))

disp('Final smoothed-State Covariance')
disp(P(end))
