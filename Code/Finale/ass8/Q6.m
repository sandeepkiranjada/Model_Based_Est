clear
clc
format long
load('measdata_pfexample02')

%% Memory allocations

Ns = 400;
k = size(zkhist,1);
xhathist = zeros(k+1,1);
xsmoothhist = zeros(k+1,1);
Phist = zeros(k+1,1);
Psmoothhist = zeros(k+1,1);
Chi= zeros(Ns,k+1);
Chi_new = zeros(Ns,k+1);
h = zeros(Ns,k+1);
wtil = zeros(Ns,k+1);
wtil_end = zeros(Ns,1);
log_wtil = zeros(Ns,k+1);
wtiltil = zeros(Ns,k+1); 
w = zeros(Ns,k+1);
v = zeros(Ns,k);
c = zeros(Ns+1,1);
measerr = zeros(Ns,1);
new_measerr = zeros(Ns,1);

Resampling_flag = 1;

%% Initial Values

xhathist(1,1) = xhat0;
Phist(1,1) = P0;
xsmoothhist(1,1) = xhat0;
Psmoothhist(1,1) = P0;

%% Particle Filter

% Sample from N(x0,P0)
Sx0           = chol(P0)';
Chi(:,1)      = xhat0 + Sx0*randn(Ns,1);
w(:,1)        = ones(Ns,1)*(1/Ns);
Svj  = chol(Q)';

for idx=1:k
    % Sample from N(0,Q)
    v(:,idx) = Svj*randn(Ns,1);
    Chi(:,idx+1) = 2*atan(Chi(:,idx)) + 0.5*cos(pi*idx/3) + v(:,idx);
    h(:,idx+1) = Chi(:,idx+1) + Chi(:,idx+1).^2 + Chi(:,idx+1).^3;
    log_wtil(:,idx+1)= log(w(:,idx))-0.5.*(ones.*zkhist(idx,1)-h(:,idx+1)).*inv(R).*(ones.*zkhist(idx,1)-h(:,idx+1));
    [~,imax] = max(log_wtil(:,idx+1));
    wtiltil(:,idx+1)=exp(log_wtil(:,idx+1)-log_wtil(imax,idx+1));
    w(:,idx+1) = wtiltil(:,idx+1)./sum(wtiltil(:,idx+1));
    xhathist(idx+1,1)  = sum(w(:,idx+1).*Chi(:,idx+1));
    Phist(idx+1,1) = sum(w(:,idx+1).*(Chi(:,idx+1)-xhathist(idx+1,1)).*(Chi(:,idx+1)-xhathist(idx+1,1)));
    
    if Resampling_flag
        % finding coefficient c(i) i=1...Ns+1
        c(1,1)=0; c(Ns+1,1)= 1+ 10^-10;
        for ii=2:Ns, c(ii,1)=sum(w(1:ii-1,idx+1)); end
        % Updating Chi and w
        eta(1) = randn(1,1)/Ns;
        for ll = 1:Ns
            eta(ll)=eta(1)+(ll-1)./Ns;
            for jj = 1:Ns
                if eta(ll)<c(jj+1)
                    break;
                end
            end
            Chi_new(ll,idx+1) = Chi(jj,idx+1);
        end
        
        % Resampling
        Chi(:,idx+1) = Chi_new(:,idx+1);
        w(:,idx+1) = ones(Ns,1)*(1/Ns);
    end
end

%% Fixed-Interval Particle smoothing

for idx1=1:Ns
    measerr(idx1,1) = (-0.5*sum((zkhist-h(idx1,2:end)').^2/R)') ;
end
[~,imax] = max(measerr(:,1));
new_measerr(:,1) = measerr(:,1) - measerr(imax,1);
wtil_end(:,1) = exp(new_measerr(:,1));
wtil_end(:,1) = wtil_end(:,1)./sum(wtil_end(:,1));
                
for idx2=1:k  
    xsmoothhist(idx2,1) = sum(wtil_end(:,1).*Chi(:,idx2));
    Psmoothhist(idx2,1) = sum(wtil_end(:,1).*(Chi(:,idx2)-xhathist(idx2,1)).*(Chi(:,idx2)-xhathist(idx2,1)));
end

xhist_end = sum(wtil_end(:,1).*Chi(:,end));
Phist_end = sum(wtil_end(:,1).*(Chi(:,end)-xhist_end).*(Chi(:,end)-xhist_end));

disp('a posterior computed estimation error covariance at final sample')
disp(Phist(end,1))

disp('a posterior computed estimation error covariance at final sample')
disp(Phist_end(end,1))


%% EKF

xhathist_EKF = ekf(zkhist,xhat0,P0,Q,R);

%% Plot

kk=0:1:100;
figure
plot(kk,xsmoothhist,'-b*',kk,xhathist,'-.ro',kk,xhathist_EKF,'-.kx');
ylabel('X')
xlabel('Sample Index - k')
legend('Particle-Smoothing','PF','EKF','Location','best')
ylim([-15 30])
grid on