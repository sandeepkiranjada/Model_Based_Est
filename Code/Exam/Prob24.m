%% HW2 Prob 4

clc; clear; close all;

%% Data

z= [  -45.1800;
1.7900;
 -31.3800;
 26.7700;
 27.6400];

H = [   -4.9300, -1.3100, -1.5900;
13.2600, 9.7100, 30.7000;
-17.0800, -11.9100, -12.1300;
-24.0300, -2.9900, -26.9500;
-2.4000, -8.7000, 9.3900];

R = [5.9700, -0.9200, -1.1800, -7.0600, -1.7900;
-0.9200, 3.4500, 1.7100, -0.6000, -4.0500;
-1.1800, 1.7100, 1.1900, 0.5600, -1.6700;
-7.0600, -0.6000, 0.5600, 9.9200, 4.8500;
-1.7900, -4.0500, -1.6700, 4.8500, 6.8700];

%% Exam addition

z  =  z + 0.5;
R = R + eye(5)*1.5;

%% Estimation

[x_hat_vec,acc] = myWLS(R,H,z);

disp('Estimation of x, x_hat:')
disp(x_hat_vec)

disp('Norm of the violation of the first-order necessary conditions:')
disp(acc)


%% Function

function [x_hat_vec,acc] = myWLS(R,H,z)

Ra = chol(R);

Ha = Ra'\H;
za = Ra'\z;

[Qb,Rb] = qr(Ha);

zb  = Qb'*za;

N_non_zero = sum(Rb(:,end) ~= 0);

zb1 = zb(1:N_non_zero,:);
Rb1 = Rb(1:N_non_zero,:);

x_hat_vec = Rb1\zb1;


acc = norm(-(H'/R)*(z-H*x_hat_vec))/norm(-(H'/R)*z);

end 
