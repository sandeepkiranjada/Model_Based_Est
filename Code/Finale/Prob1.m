%% Assignment 7, Problem 1
clc;clear;close all;
format long

%% Inital State and Proc Noise, Etc..

N1 = 60;
N2 = 120;
t0 = 0;
tf = 3;
x0 = [-0.40; 0.85; -0.60; -1.65];
v0 =[-0.77; 1.30; 1.65];
idervflag = 1;

%% RK4 Calls

[xf1,dfprinted_dxk1,dfprinted_dvk1] = ...
             c2dnonlinear(x0,0,v0,t0,tf,N1,'fscript_ts01',idervflag);
[xf2,dfprinted_dxk2,dfprinted_dvk2] = ...
             c2dnonlinear(x0,0,v0,t0,tf,N2,'fscript_ts01',idervflag);


%% Truth Model
         
   A = ...
     [-0.43256481152822, -1.14647135068146,  0.32729236140865, ...
                    -0.58831654301419;...
      -1.66558437823810,  1.19091546564300,  0.17463914282092, ...
                     2.18318581819710;...
       0.12533230647483,  1.18916420165210, -0.18670857768144, ...
                    -0.13639588308660;...
       0.28767642035855, -0.03763327659332,  0.72579054829330, ...
                     0.11393131352081];
   D = ...
     [ 1.06676821135919,  0.29441081639264, -0.69177570170229;...
       0.05928146052361, -1.33618185793780,  0.85799667282826;...
      -0.09564840548367,  0.71432455181895,  1.25400142160253;...
      -0.83234946365002,  1.62356206444627, -1.59372957644748];

sysmodel_ct = ss(A,D,eye(4),zeros(4,3));
sysmodel_dt = c2d(sysmodel_ct,(3),'zoh');
[dfprinted_dxk_t,dfprinted_dvk_t] = ssdata(sysmodel_dt);
xft = dfprinted_dxk_t*x0 + dfprinted_dvk_t*v0;

%% Results

disp('Error in fprinted (Truth-RK4), NRK = 60');
disp((xft-xf1))

disp('Error in dfprinted_dxk (Truth-RK4), NRK = 60');
disp((dfprinted_dxk_t-dfprinted_dxk1))

disp('Error in dfprinted_dvk (Truth-RK4), NRK = 60');
disp((dfprinted_dvk_t-dfprinted_dvk1))


disp('Error in fprinted (Truth-RK4), NRK = 120');
disp((xft-xf2))

disp('Error in dfprinted_dxk (Truth-RK4), NRK = 120');
disp((dfprinted_dxk_t-dfprinted_dxk2))

disp('Error in dfprinted_dvk (Truth-RK4), NRK = 120');
disp((dfprinted_dvk_t-dfprinted_dvk2))










