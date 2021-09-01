function [xhathist,Phist,sigmahist] = pf(xhat0,P0,Qk,Rk,zhist,Ns,Nt)

%
%  Get the problem dimensions and initialize the output arrays.
%
    N = size(zhist,1);
    Np1 = N + 1;
    nx = length(xhat0);
    nv = length(Qk); 
    xhathist = zeros(Np1,nx);
    Phist = zeros(nx,nx,Np1);
    sigmahist = zeros(Np1,nx);
    
%
%  Assign initial xhat and P to first elements of output arrays.
%
    xhathist(1,:) = xhat0;
    Phist(:,:,1) = P0;
    sigmahist(1,:) = sqrt(diag(P0))';
    zhist = [0 ; zhist ];
    
%
%  Preallocate arrays for w_i's, x_i's and v_i's
%
    w_i = zeros(Np1,Ns);
    x_i = zeros(nx,Ns,Np1);         % Chi underbar ^i in the lecture notes
    v_i = zeros(nv,Ns,Np1);         % script v underbar ^i in the lecture notes
    
%
%  Generate initial particles for x_i, initialize weights at 1/Ns, precompute Svj which uses a 
%  time-invariant Q, and precompute Rinv
%
    Sx0 = chol(Phist(:,:,1))';                                              % (37.10)
    x_i(:,:,1) = xhathist(1,:)' + Sx0*randn(nx,Ns);                         % (37.11)
    w_i(1,:) = (1/Ns)*ones(1,Ns);
    Svj = chol(Qk)';                                                        % (37.12)
    Rinv = 1/Rk;

%
%  Particle Filtering Algorithm with Traditional Resampling
%

for k = 1:N
    
    %  -------------------------------------------------------------------------------------------  %
    %                                     Dynamic Propagation
    %  -------------------------------------------------------------------------------------------  %

    %
    %  Generate Ns independent process noise samples
    %
        v_i(:,:,k) = Svj*randn(nv,Ns);                                      % (37.13)
        
    %
    %  Propagate particles using x_i(k), v_i(k)
    %
        x_i(:,:,k+1) = 2*atan(x_i(:,:,k)) + 0.5*cos(pi*k/3) + v_i(:,:,k);   % (39.1), Gamma = 1
        
    %  -------------------------------------------------------------------------------------------  %
    %                                     Measurement Update
    %  -------------------------------------------------------------------------------------------  %  
    
    %
    %  Compute for the normalized weights using natural log to avoid numerical underflow
    %
        zofx_i = x_i(:,:,k+1) + x_i(:,:,k+1).^2 + x_i(:,:,k+1).^3;
        delz = zhist(k+1) - zofx_i; 
        log_wtil_i_kp1 = log(w_i(k,:)) - 0.5*delz.^2*Rinv;                  % (39.9)
        [~,i_max] = max(log_wtil_i_kp1);                                    % (39.10)
        wtiltil_i = exp(log_wtil_i_kp1 - log_wtil_i_kp1(i_max));            % (39.11)
        w_i(k+1,:) = wtiltil_i/sum(wtiltil_i);                              % (39.12)
    
    %
    %  Compute for normalized weights using traditional exponential function
    %
    %     delz = zhist(k+1) - Hk*x_i(:,:,k+1);                                % [1xNs] vector of scalar delz's
    %     wtil_i = w_i(k,:).*exp(-0.5*delz.^2*Rinv);                          % (39.2)
    %     w_i(k+1,:) = wtil_i/sum(wtil_i);                                    % (39.3)
    %
    %  Note that this section was coded to check whether or not using the natural log is a better 
    %  scheme for calculating the normalized weights. In the few cases checked, the resulting weights 
    %  were equivalent. That being said, we're proceeding with the natural log scheme.
    
    %
    %  Compute for the a posteriori state estimate
    %
        xhathist(k+1,:) = sum( repmat(w_i(k+1,:),nx,1) .* x_i(:,:,k+1) , 2); % (39.4) 
        
    %
    %  Compute for the posteriori estimation error covariance matrix
    %
        delx = x_i(:,:,k+1) - xhathist(k+1,:)';
        Pdum = zeros(nx,nx);
        for idx_i = 1:Ns
            Pdum = Pdum + w_i(k+1,idx_i)*delx(:,idx_i)*delx(:,idx_i)';      % (39.5)
        end
        Phist(:,:,k+1) = Pdum;
        sigmahist(k+1,:) = sqrt(diag(Pdum))';
        
    %  -------------------------------------------------------------------------------------------  %
    %                                Traditional Resampling
    %  -------------------------------------------------------------------------------------------  %
    
    %
    %  Check the number of effective samples
    %
        Neff = 1/sum(w_i(k+1,:).^2);                                        % (40.8)
        
    %
    %  Resampling algorithm
    %    
    if Neff < Nt     
        %
        %  Compute for the constants
        %
            c_i = zeros(1,Ns+1);                                            % (38.18)
            c_i(Ns+1) = 1 + eps;                                            % (38.20)
            for idx_i = 2:Ns
                c_i(idx_i) = sum( w_i(k+1,1:idx_i-1) );                     % (38.19)
            end   
        %
        %  Initialize new particles
        %
            x_l = zeros(nx,Ns);   
        %
        %  Resampling x_i(k+1)                                              % Lec 38 Slides 14-15
        %
            eta_l = zeros(1,Ns);
            eta_l(1) = rand(1,1)/Ns;
            
            for idx_l = 1:Ns    
                %
                %  Compute eta_l = eta_1 + (l-1)/Ns
                %
                    eta_l(idx_l) = eta_l(1) + (idx_l - 1)/Ns;               % ( 
                %
                %  Reproduce high-weight x_i particles
                %
                    idx_i = 1;
                    while eta_l(idx_l) > c_i(idx_i+1)
                        idx_i = idx_i + 1;
                    end
                %
                %  Assign new x_l(k+1) = x_i(k+1)
                %
                    x_l(:,idx_l) = x_i(:,idx_i,k+1);
            end
            
            %
            %  Replace x_i(k+1) with new x_l(k+1)
            %
            x_i(:,:,k+1) = x_l;     
    end
    
    %
    %  Re-initialize weights
    %
        w_i(k+1,:) = (1/Ns)*ones(1,Ns);
    
end

end