%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_                                   Problem Set #6 Problem #2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
ProblemSet = 5;
ProblemNum = 6;
flag_save_figs = 0;
format longg

%-----------------------------------------------------------------------------------------------------
%  (kf_example02b.m)
%-----------------------------------------------------------------------------------------------------

%
% Load problem matrices and data for the linear Kalman filter problem of this assignment.
kf_example02b
%
% Define some useful dimensions
n_x = length(xhat0);        % Number of states
kmax = length(zhist);       % Number of k iterations for the truth model (= number of measurements)
N_mc = 50;                  % Number of Monte Carlo simulations
%
% Define different process noise covariances for the truth simulation and the Kalman filter
Qk_true = 40;               % for all k
Qk_kf = 0.004;              % for all k

%
% Preallocate for speed
%
epsilon_all = NaN*ones(kmax,N_mc);
x_all = cell(N_mc,1);
xhat_all = cell(N_mc,1);
P_all = cell(N_mc,1);
xtilde_all = cell(N_mc,1);

%
% Monte Carlo simulation
%
for idx_mc = 1:N_mc
    %
    % Simulate a perturbed xhat(0), and a noisy measurement time history
    %
    [xtruehist,zhist_noise] = kf_truthmodel(Fk,Gammak,Hk,Qk_true,Rk,xhat0,P0,kmax);
    % Output is the true estimate time history
    xhat0_noise = xtruehist(1,:)';
    zhist_noise = [0 ; zhist_noise];
    t = [0 ; thist];
    N = length(zhist_noise);
    %
    % Estimation using a standard Kalman Filter
    %
    [xhat,P] = kf(xhat0_noise,P0,Fk,Gammak,Qk_kf,Hk,Rk,zhist_noise,n_x,N);
    % Output is estimated state time history xhat which is an nx x N matrix, and state estimation 
    % covariance time history P which is an N x 1 cell array where each cell contains a 2x2 
    % covariance matrix.
    %
    % Compute for the error xtilde
    %
    xtilde = xhat' - xtruehist;
    %
    % Compute for the non-dimensional estimation error statistic, epsilon(k)
    %   
    xtildecell = mat2cell(xtilde,ones(length(xtilde),1),2); % Convert xtilde vector to cell array
    epsilon = cell2mat(cellfun(@(x,P) x*inv(P)*x',xtildecell,P,'un',0)); % Cell-wise computation of epsilon(k)
    %
    % Store data from this MC run
    %
    x_all{idx_mc} = xtruehist;
    xhat_all{idx_mc} = xhat;
    P_all{idx_mc} = P;
    xtilde_all{idx_mc} = xtilde;
    epsilon_all(:,idx_mc) = epsilon(2:end);
end
%
% Compute for the test statistic epsilon_bar(k) which is a kmax x N_mc vector.
%
epsilon_bar = (1/N_mc)*sum(epsilon_all,2);
%
% Compute for the upper and lower limits
%
alpha = 0.01;
r1 = (1/N_mc)*chi2inv(alpha/2,N_mc*n_x);
r2 = (1/N_mc)*chi2inv(1-alpha/2,N_mc*n_x);
%
% Check if epsilon_bar(k) values are within the limits 99% of the time
%
check = epsilon_bar > r1 & epsilon_bar < r2;
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_                                              Comments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For the case simulated where Qk is significantly different for the truth
% simulation, and for the Kalman filter, the filter fails the consistency
% evaluation.
%
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



