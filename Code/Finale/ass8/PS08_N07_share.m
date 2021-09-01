%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_                                 Problem Set #8 Problem #7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear; 
close all;
format longg

%
%  Load input data
%
%    load measdata_pfexample; meas2ylims = [-250 300];
   load measdata_pfexample02; meas2ylims = [-150 150];

%
%  Solve the problem using PF (Assignment #8 Problem #3)
%
   Ns = 100;   % Number of particles
   Nt = 50;    % Threshold for effective number of particles
   [xhathist_pf,Phist_pf,sigmahist_pf] = pf(xhat0,P0,Q,R,zkhist,Ns,Nt);
   
%
%  Solve the problem using EKF (Assignment #8 Problem #4)
%
   [xhathist_ekf,Phist_ekf] = ekf(xhat0,P0,Q,R,zkhist);

%
%  Solve the problem using PF (Assignment #8 Problem #4)
%
   [xstarhist,Pstar] = fips(xhat0,P0,Q,R,zkhist,Ns);  
   
%
%  Plotting the state time histories
%
   N = size(zkhist,1);
   tvec = 0:N;
   plot(tvec,xhathist_ekf,':ob','LineWidth',1,'MarkerSize',4); hold on; grid on;
   plot(tvec,xhathist_pf,':or','LineWidth',1,'MarkerSize',5);
   plot(tvec,xstarhist,':oc','LineWidth',1,'MarkerSize',3);
   legend('EKF','PF','FIPS');
   xlim([0 N]);
   xlabel('Time Index');
   ylabel('x(t)');
    
%
%  Comparing with actual measurements
%
   figure()
   zhist_pf   =  xhathist_pf +  xhathist_pf.^2 +  xhathist_pf.^3;
   zhist_ekf  = xhathist_ekf + xhathist_ekf.^2 + xhathist_ekf.^3;
   zhist_fips =    xstarhist +    xstarhist.^2 +    xstarhist.^3; 
   N = size(zkhist,1);
   tvec = 0:N;
   plot(tvec,[NaN ; zkhist],':dk','LineWidth',1,'MarkerSize',6); hold on; grid on;
   plot(tvec,zhist_ekf,':ob','LineWidth',1,'MarkerSize',4);
   plot(tvec,zhist_pf,':or','LineWidth',1,'MarkerSize',5); 
   plot(tvec,zhist_fips,':oc','LineWidth',1,'MarkerSize',3);
   legend('Actual Measurements','EKF','PF','FIPS');
   xlim([0 N]); 
   xlabel('Time Index');
   ylabel('z(x)');

   ylim(meas2ylims);
