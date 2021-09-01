function [xhathist,Phist] = ekf(xhat0,P0,Q,R,zkhist)

%
%  Get the problem dimensions and initialize the output arrays.
%
    N = size(zkhist,1);
    Np1 = N + 1;
    nz = length(zkhist(1,:));
    nx = length(xhat0);
    xhathist = zeros(Np1,nx);
    Phist = zeros(Np1,1);
%
%  Assign initial xhat and P to first elements of output arrays.
%
    xhathist(1) = xhat0;
    Phist(1) = P0;
    sigmahist(1) = sqrt(P0);
    zkhist = [0 ; zkhist ];

for k = 1:N-1

    xhatk = xhathist(k);
    Pk = Phist(k);

    %
    % Dynamic propagation
    %
        F = 2/(1+xhatk^2);
        xbarkp1 = 2*atan(xhatk) + 0.5*cos(pi*k/3);
        Pbarkp1 = F*Pk*F' + Q;  % Note that Gamma = 1

    %
    % Measurement update
    %
        zbarkp1 = xbarkp1 + xbarkp1^2 + xbarkp1^3;
        Hkp1 = 1 + 2*xbarkp1 + 3*xbarkp1^2;
        nukp1 = zkhist(k+1) - zbarkp1;
        Skp1 = Hkp1*Pbarkp1*Hkp1' + R;
        Wkp1 = Pbarkp1*Hkp1'/Skp1;
        dum = eye(nx) - Wkp1*Hkp1;
        xhatkp1 = xbarkp1 + Wkp1*nukp1;
        Pkp1 = dum*Pbarkp1*dum' + Wkp1*R*Wkp1';

    %
    % Store variables
    %
        xhathist(k+1) = xhatkp1;
        Phist(k+1) = Pkp1;

end

end