%  kf_Q.m
%  
clear
clc
%  This Matlab script loads into your work space the problem matrices
%  and data for the linear Kalman filter example that was presented
%  in class.
%
   Fk     = [  0.81671934103521,  0.08791146849849;...
              -3.47061412053765,  0.70624978972000];     % for all k
   Gammak = [  0.00464254201630;...
               0.08791146849849];                        % for all k
   Hk     = [  2.00000000000000,  0.30000000000000];     % for all k
%
   Qk     =    4.00000000000000;                         % for all k
   Rk     =    0.01000000000000;                         % for all k
%
   xhat0   = [  0.20000000000000;...
               -2.50000000000000];
   P0      = [  0.25000000000000,  0.08000000000000;...
                0.08000000000000,  0.50000000000000];
   kmax    = 50;
   NumMC   = 50;
   xtrueMC = zeros(kmax+1,2,NumMC);
   zMC = zeros(kmax,NumMC);
%
%  Note that the following arrays define the sample times 
%  [t(1);t(2);t(3);...;t(50)] and the measurements 
%  [z(1);z(2);z(3);...;z(50)]
%
   thist = [1:kmax]'*0.1;
   for i=1:NumMC
       [xtruehist,zhist] = kf_truthmodel(Fk,Gammak,Hk,Qk,Rk,xhat0,P0,kmax);
       xtrueMC(:,:,i)=xtruehist;
       zMC(:,i)=zhist;
       
       for k=1:50
           if k==1
               xbar(:,1)=Fk*xhat0;
               Pbar(:,:,1)=Fk*P0*Fk'+Gammak*Qk*Gammak';
               Inno(:,1)=zhist(1,1)-Hk*xbar(:,1);
               S(:,1)=Hk*Pbar(:,:,1)*Hk'+Rk;
               W(:,1)=Pbar(:,:,1)*Hk'/S(:,1);
               xhat(:,1)=xbar(:,1)+W(:,1)*Inno(:,1);
               P(:,:,1)=Pbar(:,:,1)-W(:,1)*S(:,1)*W(:,1)';
               std11(1,1)=sqrt(P(1,1,1));
               std22(1,1)=sqrt(P(2,2,1));
           else
               xbar(:,k)=Fk*xhat(:,k-1);
               Pbar(:,:,k)=Fk*P(:,:,k-1)*Fk'+Gammak*Qk*Gammak';
               Inno(:,k)=zhist(k,1)-Hk*xbar(:,k);
               S(:,k)=Hk*Pbar(:,:,k)*Hk'+Rk;
               W(:,k)=Pbar(:,:,k)*Hk'/S(:,k);
               xhat(:,k)=xbar(:,k)+W(:,k)*Inno(:,k);
               P(:,:,k)=Pbar(:,:,k)-W(:,k)*S(:,k)*W(:,k)';
               std11(k,1)=sqrt(P(1,1,k));
               std22(k,1)=sqrt(P(2,2,k));
           end
       end
       xhatMC(:,:,i)=[xhat0,xhat]';
   end
% xhat=[xhat0,xhat];
xtildat=xtrueMC-xhatMC;

xtilad10_mean=mean(xtildat(10+1,:,:),3)
xtilad35_mean=mean(xtildat(35+1,:,:),3)
xtilad10_cov=cov(permute(xtildat(10+1,:,:),[3 2 1]))
xtilad35_cov=cov(permute(xtildat(35+1,:,:),[3 2 1]))