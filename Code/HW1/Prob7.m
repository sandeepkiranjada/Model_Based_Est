clc;clear;close all;

A = [8 7;3 -14;-13 16];
R_and_0 = A.*0;

a1 = A(:,1);


r11 = norm(a1);

r11 = -1*sign(A(1,1))*r11;

R_and_0(1,1) = r11;
r1 = R_and_0(:,1);
H1 = householder(a1,r1);

A2 = H1*A;

a2 = A2((2:end),2);


r22 = norm(a2);
r22 = -1*sign(A(2,2))*r22;

H22 = householder(a2,[r22;0]);

H2 = eye(length(H1));

H2(2:end,2:end)= H22;

myR=H2*H1*A;
myQ=(H2*H1)';

[Q,R]=qr(A);

myQ-Q
myR-R

