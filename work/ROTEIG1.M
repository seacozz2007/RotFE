function [r1,lambda1]=roteig(Rot,crit)
%[r1,lambda1]=roteig( RotorModel );
%
% Compute the eigevectors and eigenvalues
% of a Rotor system,
% Input:
%   RotorModel : a structure (or an object of type rotor)
%                 which was created by one of the two commands:
%						(1)  RotorModel=rotfe('model');
%						(2)  RotorModel=load('model.mtx','-mat');
%    Whithin the model the following fields are being used
%      M,K,D,G,W (i.e. RotorModel.M, RotorModel.K ...)
%
% By:		I. Bucher
% Date	15-8-98
% 
% Example:
%				
nM=length(Rot.M); 
       z=zeros(nM);
       
      A=[ Rot.M z ; z Rot.K]; B=[Rot.W*Rot.G+Rot.D Rot.K ; -Rot.K z ]; % see Lee page 156
      [r1 lambda]=eig(full(B),full(A));                   % 2N modes
      lambda1=-diag(lambda);   % a matter of definition
      %It seems a bit difficult to remove the redundent rigid body modes
      % which are created by the double-size (1st order) eignvalue problem being solved
      %therfore an accurate SVD based range space method is being used, IB 7-99
      
      [dd jj]=sort(abs(lambda1));
      lambda1=lambda1(jj); r1=r1(1:nM,jj);
      