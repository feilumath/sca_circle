
function X1 = move1step(X,stoCA)
% StoCA move 1 step
N = stoCA.N; 
K = stoCA.K;  
edges = stoCA.edges; 
X1 = X;

for k=1:N
    nhbr_ind = stoCA.Ind_nhbr_k(k); 
    Xk_nhbr  = X(nhbr_ind);
    nhbr_density = histcounts(Xk_nhbr,edges);      % we can put more weights at k later
    nhbr_density = nhbr_density'/length(Xk_nhbr); 
    P_k      = stoCA.TMat*nhbr_density;   % a probability vector 1xK
    % mnvec   = 1+mnrnd(K-1,P_k,1);   % multi Normial vector X=(X_1,...,X_K), sum X_i = K; but P(X_k= k-1) !=p_k
    X1(k)   = gendist(P_k',1,1);   
end    
    
end

