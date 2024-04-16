function [local_p_all,Xt1_color] = all_local_density(Xt_all,stoCA)
% Assemble data to local population density phi(X_i,V_i)(t) \in \R^K. It is what needed for Tmat esitmation (LSE or MLE) 
% Size of phi(X_i,V_i)(t):  [K,N,tN] 
% local_p   size: K,N,tN   %% local_p(:,1,1) will be tranposed to  1 x K 


M = length(Xt_all); % Xt_all= cell(1,M): each cell is a trajectory of data: they can be of different length
K = stoCA.K; 
local_p_all = cell(1,M); 
Xt1_color   = cell(1,M);
% fprintf('\n Can use parallel to accelerate computation of local densities \n'); 
parfor m=1:M
    local_p_all{m} = local_density1Traj(Xt_all{m},K,stoCA); 
end


function local_p = local_density1Traj(Xt,K,stoCA)
% asseme local population density phi(X_i,V_i)(t) \in \R^K for each
% trajectory 
[N,tN]     = size(Xt); 
local_p    = zeros(K,N,tN);  %% local_p(:,1,1) will be tranposed to  1 x K  
edges      = stoCA.edges; 
% xt1color   = zeros(N,tN-1);

for t = 1:tN
    X0 = Xt(:,t);  
    for nn=1:N
        nhbr_ind       = stoCA.Ind_nhbr_k(nn);
        Xk_nhbr        = X0(nhbr_ind);
        nhbr_density   = histcounts(Xk_nhbr,edges);      % we can put more weights at k later
        nhbr_density   = nhbr_density/length(Xk_nhbr);
        local_p(:,nn,t) = nhbr_density;
    end
    % xt1color   = Xt(:,t+1); 
end
