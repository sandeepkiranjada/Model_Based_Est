%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_                                   Problem Set #6 Problem #2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
ProblemSet = 6;
ProblemNum = 2;
flag_save_figs = 0;
format longg

%
% Load problem matrices and data for the linear Kalman filter problem of this assignment.
kf_example03a
z = [0 ; zhist];
t = [0 ; thist];
N = length(z);
n_x = length(xhat0);
n_v = length(Qk);

%
% Estimation using standard Kalman Filter
%
[xhat,P,xbar,Pbar] = kf(xhat0,P0,Fk,Gammak,Qk,Hk,Rk,z,n_x,N);

%
% Backward smoothing steps
%
xstar = NaN*ones(n_x,N);
Pstar = cell(N,1);
xstar(:,N) = xhat(:,N);
Pstar{N} = P{N};

for k = N-1:-1:1
    
    Ck = P{k}*Fk'*(Pbar{k+1}\eye(n_x));
    xstar(:,k) = xhat(:,k) + Ck*(xstar(:,k+1) - xbar(:,k+1));
    Pstar{k} = P{k} - Ck*(Pbar{k+1} - Pstar{k+1})*Ck';
    
end

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
plot(t,xhat(1,:),'-^k','MarkerSize',6); hold on; grid on;
plot(t,xhat(2,:),'->k','MarkerSize',6);
plot(t,xhat(3,:),'-vk','MarkerSize',6);
%
% Plot SQIF results
plot(t,xstar(1,:),'-^b','MarkerFaceColor','b','MarkerSize',2); hold on; grid on;
plot(t,xstar(2,:),'->b','MarkerFaceColor','b','MarkerSize',2);
plot(t,xstar(3,:),'-vb','MarkerFaceColor','b','MarkerSize',2);

xhat(:,11)
xstar(:,11)
P{11}
Pstar{11}
P{11}-Pstar{11}
eig(P{11})-eig(Pstar{11})

sqrt(diag(P{11}))
sqrt(diag(Pstar{11}))
sqrt(diag(P{11}))-sqrt(diag(Pstar{11}))

function [xhat,P,xbar,Pbar] = kf(xhat0,P0,Fk,Gammak,Qk,Hk,Rk,z,n_x,N)

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

%
% set(gca,'FontSize',12);
%
% post_save_figs

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_                                              Comments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The smoothed state time history estimate plots don't appear to be
% smoother at the full scale used in the plots of the time history.
% If we zoom in however, we see that the smoothed time history is
% "smoother" than the filter estimate time history.
%
% There is a 0.18113, 0.003733, 0.327480 1-sigma difference between 
% the filtered estimate xhat (10) and the smoothed estimate xstar(10), 
% where the variance corresponding to xhat(10) is higher than that 
% of xstar(10).
%
