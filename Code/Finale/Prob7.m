clc;clear;close all;
format longg

%% Inital State and Proc Noise, Etc..

cart_EKF_meas

NRK = 120;
idervflag = 1;

N = length(thist);

x_hat = zeros(5,N);

sig_rhoa = 0.002;
sig_rhob = 0.002;
del_t = thist(2)-thist(1);
qsteertil = 0.1;
qspeedtil = 21.25;
tausteer = 0.25;
l_a = -1;
l_b = 1;
sasb = sig_rhoa*sig_rhob;

c = zhist(1,1)^2 - zhist(1,2)^2;
d = c - l_a^2 + l_b^2;
x2_0 = d/2/(l_b-l_a);
x3_0 = sqrt(zhist(1,1)^2 - (l_a-x2_0)^2);

c = zhist(2,1)^2 - zhist(2,2)^2;
d = c - l_a^2 + l_b^2;
x2_1 = d/2/(l_b-l_a);
x3_1 = sqrt(zhist(2,1)^2 - (l_a-x2_1)^2);

vr_1 =  [ x2_1-x2_0 ; x3_1-x3_0 ]./del_t;

x1_0 = atan2(vr_1(2),vr_1(1));

x4_0 = 0;
x5_0 = 2.1;
x_0 = [x1_0 ; x2_0 ; x3_0 ; x4_0 ; x5_0];

P = zeros(5,5,N);
P1_11 = 2*sasb/(del_t^2)/(norm(vr_1)^2);
P1_12 = sasb;
P1_13 = sasb;
P1_14 = qsteertil*tausteer/2;
P1_15 = 2*sasb/(del_t^2);

P1 = diag([P1_11 P1_12 P1_13 P1_14 P1_15]);

P(:,:,1) = P1;

x_hat(:,1) = x_0 + randn(5,1) .* sqrt([P1_11 P1_12 P1_13 P1_14 P1_15]');

qsteer = qsteertil/del_t;
qspeed = qspeedtil/del_t;

Q = diag([qsteer qspeed]);

v = diag(sqrt([qsteertil qspeedtil]))*randn(2,N);

R = diag([sig_rhoa^2 sig_rhoa^2]);

eps_nu = zeros(1,N-1);


for k = 1:N-1

    [x_bar_kp1,Fk,Gammak] = ...
                 c2dnonlinear(x_hat(:,k),0,v(:,k),thist(k),thist(k+1),NRK,...
                 'fscript_cart',idervflag);
             
    P_bar_kp1 = Fk*P(:,:,k)*Fk' + Gammak*Q*Gammak';
    
    [z_bar_kp1,H_kp1] = h_cart(x_bar_kp1,idervflag);
    
    nu_kp1 = zhist(k+1,:)' - z_bar_kp1;
    
    S_kp1 = H_kp1*P_bar_kp1*H_kp1' + R;
    W_kp1 = P_bar_kp1*H_kp1'/S_kp1;
    
    x_hat(:,k+1) = x_bar_kp1 + W_kp1*nu_kp1;
    
    P(:,:,k+1) = P_bar_kp1 - W_kp1*S_kp1*W_kp1';
    
    eps_nu(k) = nu_kp1'/S_kp1*nu_kp1;
    
end

%%

figure;
plot(x_hat(2,:),x_hat(3,:),'b','Linewidth',1.5);hold on;grid on;
xlabel('x_2, m')
ylabel('x_3, m')
axis equal
xlim([-3 2])
ylim([1 6])
ax = gca;
ax.LineWidth = 1;
ax.GridColor = [0 0 0];
ax.MinorGridColor = 'k';

figure;
subplot(3,1,1)
plot(thist,x_hat(1,:),'b','Linewidth',1.5);hold on;grid on;
xlabel('time, s')
ylabel('x_1, m')
ax = gca;
ax.LineWidth = 1;
ax.GridColor = [0 0 0];
ax.MinorGridColor = 'k';

subplot(3,1,2)
plot(thist,x_hat(4,:),'b','Linewidth',1.5);hold on;grid on;
xlabel('time, s')
ylabel('x_4, rad')
ax = gca;
ax.LineWidth = 1;
ax.GridColor = [0 0 0];
ax.MinorGridColor = 'k';

subplot(3,1,3)
plot(thist,x_hat(5,:),'b','Linewidth',1.5);hold on;grid on;
xlabel('time, s')
ylabel('x_5, m/s')
ax = gca;
ax.LineWidth = 1;
ax.GridColor = [0 0 0];
ax.MinorGridColor = 'k';



figure;
plot(thist(2:end),eps_nu,'.b','Linewidth',1.5);hold on;grid on;
plot(thist(2:end),ones(1,N-1)*2,'k','Linewidth',1.5);hold on;grid on;
xlabel('time, s')
ylabel('\epsilon_\nu')
ax = gca;
ax.LineWidth = 1;
ax.GridColor = [0 0 0];
ax.MinorGridColor = 'k';

