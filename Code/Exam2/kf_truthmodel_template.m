function [xtruehist,zhist] = kf_truthmodel(F,Gamma,H,Q,R,xhat0,...
                                           P0,kmax)
%
%  Copyright (c) 2002 Mark L. Psiaki.  All rights reserved.  
% 
%  This function performs a truth-model Monte-Carlo simulation for
%  the discrete-time stochastic system model:
%
%     x(k+1) = F*x(k) + Gamma*v(k)
%     z(k)   = H*x(k) + w(k)
%
%  Where v(k) and w(k) are uncorrelated, zero-mean, white-noise
%  Gaussian random processes with covariances E[v(k)*v(k)'] = Q and
%  E[w(k)*w(k)'] = R.  The simulation starts from a true x(0)
%  that is drawn from the Gaussian distribution with mean xhat0
%  and covariance P0.  The simulation starts at time k = 0 and
%  lasts until time k = kmax.
%  
%
%  Inputs:
%
%    F          The (nx)x(nx) state transition matrix for the time-
%               invariant linear system.
%
%    Gamma      The (nx)x(nv) process noise gain matrix for the time-
%               invariant linear system.
%
%    H          The (nz)x(nx) output gain matrix for the time-
%               invariant linear system.
%
%    Q          The (nv)x(nv) symmetric positive definite process noise 
%               covariance matrix.
%
%    R          The (nz)x(nz) symmetric positive definite measurement 
%               noise covariance matrix.
%
%    xhat0      The (nx)x1 initial state estimate.
%
%    P0         The (nx)x(nx) symmetric positive definite initial state
%               estimation error covariance matrix.
%
%    kmax       The maximum discrete-time index of the simulation.
%
%  Outputs:
%
%    xtruehist  = [x(0)';x(1)';x(2)';...;x(kmax)'], the (kmax+1)x(nx) 
%               array that stores the truth time history for the state 
%               vector.  Note that Matlab does not allow zero indices
%               in its arrays.  Thus, x(0) = xtruehist(1,:)',
%               x(1) = xtruehist(2,:)', etc.
%
%    zhist      = [z(1)';z(2)';z(3)';...;z(kmax)'], the (kmax)x(nz)
%               array that stores the measurement time history.
%               z(1) = zhist(1,:)', z(2) = zhist(2,:)', etc.  Note
%               that the state vector xtruehist(j+1,:)' and 
%               the measurement vector zhist(j,:)' correspond
%               to the same time.
%

%
%  Get problem dimensions and set up the output arrays.
%
   [nx,nv] = size(Gamma);
   nz = size(H,1);
   xtruehist = zeros((kmax+1),nx);
   zhist = zeros(kmax,nz);
%
%  Calculate the appropriate matrix square root of P0 for use in
%  randomly generating an initial state.
%
   ????;
%
%  Generate the initial truth state with the aid of a random number 
%  generator and store it.
%
   x0 = ????;
   xtruehist(1,:) = x0';
%
%  Calculate the appropriate matrix square roots of Q and R for use in
%  randomly generating the process and measurement noise vectors.
%
   ????;
   ????;
%
%  Iterate through the kmax samples, propagating the state with
%  randomly generated process noise of the correct covariance and
%  corrupting the measurements with randomly generated measurement
%  noise of the correct covariance.
%
   xk = x0;
   for k = 1:kmax
      xkm1 = xk;
      vkm1 = ????;
      xk = F*xkm1 + Gamma*vkm1;
      wk = ????;
      zk = H*xk + wk;
      kp1 = k + 1;
      xtruehist(kp1,:) = xk';
      zhist(k,:) = zk';
   end