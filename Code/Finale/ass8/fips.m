function [xstarhist,Pstarhist] = fips(xhat0,P0,Q,R,zhist,Ns)

%
%  Get the problem dimensions and initialize the output arrays.
%
    N = size(zhist,1);
    Np1 = N + 1;
    nx = length(xhat0);
    nv = length(Q); 
    xhathist = zeros(Np1,nx);
    xstarhist = zeros(Np1,nx);
    Phist = zeros(nx,nx,Np1);
    Pstarhist = zeros(nx,nx,Np1);
    zhist = [0 ; zhist];
    
%
%  Assign initial xhat and P to first elements of output arrays.
%
    xhathist(1,:) = xhat0;
    Phist(:,:,1) = P0;

%
%  Define number of particles, preallocate arrays for w_i's, x_i's and v_i's
%
    w_i = zeros(Np1,Ns);
    x_i = zeros(nx,Ns,Np1);           % Chi underbar ^i in the lecture notes
    v_i = zeros(nv,Ns,Np1);           % script v underbar ^i in the lecture notes

%
%  Reshape zhist into a 3D matrix
%
    z3D = zeros(nv,Ns,Np1);
    for idx_i = 1:Ns
      z3D(1,idx_i,:) = zhist;
    end
    
%
%  Generate initial particles for x_i, initialize weights at 1/Ns, precompute Svj which uses a 
%  time-invariant Q, and precompute Rinv
%
    Sx0 = chol(Phist(:,:,1))';                                              % (37.10)
    x_i(:,:,1) = xhathist(1,:)' + Sx0*randn(nx,Ns);                         % (37.11)
    w_i(1,:) = (1/Ns)*ones(1,Ns);
    Svj = chol(Q)';                                                         % (37.12)
    Rinv = 1/R;    
    
%
%  Compute Ns simulated state time histories based on zero data. The sampling amounts to truth-model 
%  simulation of the system driven by simulated random process noise.
%

   for k = 1:N
      %
      %  Generate Ns independent process noise samples
      %
         v_i(:,:,k) = Svj*randn(nv,Ns);                                      % (37.13)
      
      %
      %  Propagate particles using x_i(k), v_i(k)
      %
         x_i(:,:,k+1) = 2*atan(x_i(:,:,k)) + 0.5*cos(pi*k/3) + v_i(:,:,k);   % (39.1), Gamma = 1
   end

%
%  Compute weights, one for each time history.
%
   zofx_i = x_i(:,:,:) + x_i(:,:,:).^2 + x_i(:,:,:).^3;
   delz = z3D - zofx_i;
   log_w_til_i = zeros(Ns,1);
   for idx_i = 1:Ns
      log_w_til_i(idx_i) = - 0.5*sum((delz(:,idx_i,:).^2)*Rinv);
   end
   [~,i_max] = max(log_w_til_i);                                           % (39.10)
   normlog_w_til_i = log_w_til_i - log_w_til_i(i_max);
   w_til_i = exp(normlog_w_til_i);                                         % (39.11)
   w_N = w_til_i/sum(w_til_i);                                             % (39.12)

% %
% %  Use the Ns randomly simulated points at the final sample, and the
% %  corresponding weights to compute Monte-Carlo approximations of the a
% %  posteriori mean and error covariance
% %
%    xhathist(N,:) = sum(w_N'.*x_i(:,:,N));
%    delx = x_i(:,:,N) - xhathist(N,:);
%    Phist(:,:,N) = sum(w_N'.*(delx.^2));
% 
% %
% %  Set x*(N) = xhat(N) and P*(N) = P(N) 
% %
%    xstarhist(N,:) = xhathist(N,:);
%    Pstarhist(:,:,N) = Phist(:,:,N);
   
%
%  Compute for the smoothed estimates for k = 1:N-1
%
   for k = 1:Np1
      xstarhist(k,:) = sum(w_N'.*x_i(:,:,k));
      delx = x_i(:,:,k) - xhathist(N,:);
      Phist(:,:,k) = sum(w_N'.*(delx.^2));
   end

end