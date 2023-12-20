
function [lkhd,lkhd_1step] = move1step_infer(X0,X1,stoCA)
% compuate the likelihood of moving from X0 to X1 with the StoCA.TMat 
% For future: we return the probablity P(X1|X0)
N = stoCA.N; 
K = stoCA.K;  % we can keep all P_k's  
edges = stoCA.edges; 

lkhd = 1;
Px0tox1 = zeros(N,1);
for k=1:N
    nhbr_ind = stoCA.Ind_nhbr_k(k); 
    Xk_nhbr  = X0(nhbr_ind);
    nhbr_density = histcounts(Xk_nhbr,edges);      % we can put more weights at k later
    nhbr_density = nhbr_density'/length(Xk_nhbr); 
    P_k          = stoCA.TMat*nhbr_density;   % a probability vector 1xK
    Px0tox1(k)   = P_k(X1(k)); 
    lkhd         = lkhd*Px0tox1(k);    
end    
    
end

