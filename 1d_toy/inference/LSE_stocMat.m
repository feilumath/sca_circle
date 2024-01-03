function [lse_stoc,lse_all,Tmat_lse] = LSE_stocMat(local_p_all,Xt_all,stoCA_par,lse_stocON)
% Estimate the tranform stochastic matrix by LSE: only the first (K-1) row
% of the matrix as free parameter 
%  Abar = B otimes Amat; 
%  bbar (k) = sum  (1+ c^k - c^K) phi


M = length(local_p_all); 
K = stoCA_par.K; 
 
bbar = zeros( (K-1)*K,1);  

%% LSE estimator: each row, min (1- P(T)): phi= local_p_all{m}    || 1_{k} - T* phi ||^2 >> phi*phi' T = phi; phi KxN for each t
Amat = zeros(K,K);  % phi*phi'   with phi KxN for each t 
bvec = zeros(K,K);  % phi*1_[Xt1](k)             1_[Xt1] size 1xK   >>> lse is column stochastic: sum_j T(k,j) =1

lse_all = cell(1,M);
parfor m = 1:M
    Amat1traj = zeros(K,K); bvec1traj = zeros(K,K); 
          
    bbar1traj = zeros( (K-1)*K,1);  
    
    phi1traj  = local_p_all{m}; 
    Xt        = Xt_all{m};  % Xt size = [N,tN] 
    [K1,~,tN] = size(phi1traj); 
    if K1 ~= K; fprintf('\n local_p_all size miss match\n'); end
    for t=1:tN-1
        Xt1       = Xt(:,t+1); % size Nx1, with value in [1,K]; 
        N = length(Xt1); 
        indicator = zeros(N,K);   % N x K
        for n= 1:N
            indicator(n,Xt1(n)) = 1; 
        end 
        Tmat_lse = phi1traj(:,:,t);   % KxN
        Amat1traj = Amat1traj+ Tmat_lse*Tmat_lse'; 
        bvec1traj = bvec1traj+ Tmat_lse*indicator;  
        
        temp_bar  = Tmat_lse*(indicator(:,1:(K-1))+1- indicator(:,K)); 
        bbar1traj = bbar1traj + reshape(temp_bar,[],1); 
    end
    lse_all{m}.Amat = Amat1traj; 
    lse_all{m}.bvec = bvec1traj;
    Amat = Amat +  Amat1traj; 
    bvec = bvec + bvec1traj;
    
    lse_all{m}.bbar1traj = bbar1traj;
    bbar = bbar + bbar1traj;
end

Tmat_lse  = Amat\bvec; 
% lb = 0; ub = 1;  Tmat_lse =  lsqlin(Amat,bvec,[],[],[],[],lb,ub); 
% lse.cvec = Tmat_lse;
% lse.Amat = Amat; 
% lse.bvec = bvec; 
Tmat_lse = Tmat_lse'; 



%% LSE estimator: only the first (K-1) row as freee parameter, the last row by stochasticity
if lse_stocON ==1
    B    =  ones(K-1,K-1)+ eye(K-1);
    Abar = kron(B,Amat); lb = zeros(K*(K-1),1); ub = 1+lb;
    lse_stoc = lsqlin(Abar,bbar,[],[],[],[],lb,ub); %  Abar\bbar;
    lse_stoc = reshape(lse_stoc,[K,K-1]);
    lse_stoc = [lse_stoc, 1-sum(lse_stoc,2)];
    lse_stoc = lse_stoc';
else
    lse_stoc =[]; 
end
end

