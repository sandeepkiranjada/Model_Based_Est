%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_                                   Problem Set #6 Problem #2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
ProblemSet = 6;
ProblemNum = 2;
flag_save_figs = 0;
format longg

%-----------------------------------------------------------------------------------------------------
%  Problem (a) (kf_example03a.m)
%-----------------------------------------------------------------------------------------------------
% Load problem matrices and data for the linear Kalman filter problem of this assignment.
kf_example03a
part = 'a';
solvePS06_N02(Fk,Gammak,Hk,Qk,Rk,xhat0,P0,zhist,thist,part)

%-----------------------------------------------------------------------------------------------------
%  Problem (b) (kf_example03b.m)
%-----------------------------------------------------------------------------------------------------
% Load problem matrices and data for the linear Kalman filter problem of this assignment.
kf_example03b
part = 'b';
solvePS06_N02(Fk,Gammak,Hk,Qk,Rk,xhat0,P0,zhist,thist,part)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_                                           Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function solvePS06_N02(Fk,Gammak,Hk,Qk,Rk,xhat0,P0,zhist,thist,part)

z = [0 ; zhist];
t = [0 ; thist];
N = length(z);
n_x = length(xhat0);
n_v = length(Qk);

%
% Estimation using standard Kalman Filter
%
[xhat_kf,P_kf] = kf(xhat0,P0,Fk,Gammak,Qk,Hk,Rk,z,n_x,N);
%
% Cell-wise computation of (diag(P(k)))^(1/2)
sigma_kf = cellfun(@(x) diag(x)'.^0.5,P_kf,'un',0);
sigma_kf = cell2mat(sigma_kf);

%
% Estimation using Square-Root Information Filter
%
[xhat_sqif,P_sqif] = sqif(xhat0,P0,Fk,Gammak,Qk,Hk,Rk,z,n_x,n_v,N);
%
% Cell-wise computation of (diag(P(k)))^(1/2)
sigma_sqif= cellfun(@(x) diag(x)'.^0.5,P_sqif,'un',0);
sigma_sqif = cell2mat(sigma_sqif);

fprintf('\n-------------------------------------\n');
xhat_kf_N = xhat_kf(:,end)
xhat_sqif_N = xhat_sqif(:,end)

xhat_kf_N - xhat_sqif_N 

P_kf_N = P_kf{end}
P_sqif_N = P_sqif{end}

P_kf_N - P_sqif_N

diff_kf_sqif_xhat = xhat_kf(:,end)-xhat_sqif(:,end)
eig_P_kf = eig(P_kf{end})
eig_P_sqif = eig(P_sqif{end})
diff_eig_P = eig(P_kf{end}) - flip(eig(P_sqif{end}))
fprintf('\n-------------------------------------\n');

%
% Plotting
%
% Open figure
picsize = [1 1 30 20];
f1 = figure();
f1.Units = 'centimeters';
f1.Position = picsize;
%
% Plot KF results
subplot(3,1,1); plot(t,xhat_kf(1,:),'-^k','MarkerSize',8,'LineWidth',1); hold on; grid on; ylabel('x_1');
subplot(3,1,2); plot(t,xhat_kf(2,:),'->k','MarkerSize',8,'LineWidth',1); hold on; grid on; ylabel('x_2');
subplot(3,1,3); plot(t,xhat_kf(3,:),'-vk','MarkerSize',8,'LineWidth',1); hold on; grid on; ylabel('x_3');
%
% Plot SQIF results
subplot(3,1,1); plot(t,xhat_sqif(1,:),'--^b','MarkerFaceColor','b','MarkerSize',2,'LineWidth',1); hold on; grid on;
subplot(3,1,2); plot(t,xhat_sqif(2,:),'-->b','MarkerFaceColor','b','MarkerSize',2,'LineWidth',1); hold on; grid on;
subplot(3,1,3); plot(t,xhat_sqif(3,:),'--vb','MarkerFaceColor','b','MarkerSize',2,'LineWidth',1); hold on; grid on;
subplot(3,1,1); title(['Problem ' part ' State Estimates']);
legend('KF','SQIF');
%
% Open figure
picsize = [1 1 30 20];
f1 = figure();
f1.Units = 'centimeters';
f1.Position = picsize;
%
% Plot KF results
subplot(3,1,1); plot(t,sigma_kf(:,1),'-^k','MarkerSize',8,'LineWidth',1); hold on; grid on; ylabel('\sigma_{x_1}');
subplot(3,1,2); plot(t,sigma_kf(:,2),'-^k','MarkerSize',8,'LineWidth',1); hold on; grid on; ylabel('\sigma_{x_2}');
subplot(3,1,3); plot(t,sigma_kf(:,3),'-^k','MarkerSize',8,'LineWidth',1); hold on; grid on; ylabel('\sigma_{x_3}');
%
% Plot SQIF results
subplot(3,1,1); plot(t,sigma_sqif(:,1),'--^b','MarkerFaceColor','b','MarkerSize',2,'LineWidth',1); hold on; grid on;
subplot(3,1,2); plot(t,sigma_sqif(:,2),'--^b','MarkerFaceColor','b','MarkerSize',2,'LineWidth',1); hold on; grid on;
subplot(3,1,3); plot(t,sigma_sqif(:,3),'--^b','MarkerFaceColor','b','MarkerSize',2,'LineWidth',1); hold on; grid on;
subplot(3,1,1); title(['Problem ' part ' State Estimation Standard Deviations']);
legend('KF','SQIF');
%
end

function [xhat,P] = sqif(xhat0,P0,Fk,Gammak,Qk,Hk,Rk,z,n_x,n_v,N)
%
% Preallocate for speed
%
xhat = NaN*ones(n_x,N);
zx = NaN*ones(n_x,N);
Rxx = cell(N,1);

%
% Initialize the Kalman Filter
%
xhat(:,1) = xhat0;
P = cell(N,1);
%
P{1} = P0;
Rxx{1} = chol(inv(P0)); % sum(sum(abs(inv(P0) - P0\eye(3)))) is O(1e-17)
zx(:,1) = Rxx{1}*xhat0;
Rvv = inv(chol(Qk))'; 
Finv = inv(Fk);
Ra = chol(Rk);
RainvT = inv(Ra)'; % Which in both problems is just 1.
Ha = RainvT*Hk; % for all k

for k = 1:N-1
    %
    % Prediction Step/Dynamic Propagation
    %
    RxxFinv = Rxx{k}*Finv;
    MatA = [Rvv zeros(n_v,n_x) ; -RxxFinv*Gammak RxxFinv];
    [TaT,R] = qr(MatA);
    Tak = TaT';
    % Rbarvv{k} = R(1:n_v,1:n_v);
    % Rbarvx{k+1} = R(1:n_v,n_v+1:end);
    Rbarxxkp1 = R(n_v+1:end,n_v+1:end);
    %
    zbarvec = Tak*[zeros(n_v,1) ; zx(:,k)];
    % zbarv = zbarvec(1:n_v);
    zbarx = zbarvec(n_v+1:end);
    
    %
    % Measurement Update
    %
    MatB = [Rbarxxkp1 ; Ha];
    [TbT,R] = qr(MatB);
    Rxx{k+1} = R(1:n_x,:);
    %
    Tbkp1 = TbT';
    zvec = Tbkp1*[zbarx ; RainvT*z(k+1)];
    zx(:,k+1) = zvec(1:n_x);
    %
    Rxxinvkp1 = Rxx{k+1}\eye(n_x);
    xhat(:,k+1) = Rxxinvkp1*zx(:,k+1);
    P{k+1} = Rxxinvkp1*Rxxinvkp1';
end
end

function [xhat,P] = kf(xhat0,P0,Fk,Gammak,Qk,Hk,Rk,z,n_x,N)

%
% Preallocate for speed
%
xhat = NaN*ones(n_x,N);
xbar = NaN*ones(n_x,N);
%
P = cell(N,1);
Pbar = cell(N,1);
%
nu = NaN*ones(N,1);
W = cell(N,1);
S = NaN*ones(N,1);

%
% Initialize the Kalman Filter
%
xhat(:,1) = xhat0;
P{1} = P0;

for k = 1:N-1
    %
    % Prediction Step
    %
    xbar(:,k+1) = Fk*xhat(:,k);
    Pbar{k+1}   = Fk*P{k}*Fk' + Gammak*Qk*Gammak';
    
    %
    % Measurement Update
    %
    nu(k+1) = z(k+1) - Hk*xbar(:,k+1);
    S(k+1) = Hk*Pbar{k+1}*Hk' + Rk;
    W{k+1} = Pbar{k+1}*Hk'/S(k+1);
    dum = eye(n_x) - W{k+1}*Hk;
    %
    xhat(:,k+1) = xbar(:,k+1) + W{k+1}*nu(k+1);
    P{k+1} = dum*Pbar{k+1}*dum' + W{k+1}*Rk*W{k+1}';
    %
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_                                              Comments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%
% set(gca,'FontSize',12);
%
% post_save_figs
