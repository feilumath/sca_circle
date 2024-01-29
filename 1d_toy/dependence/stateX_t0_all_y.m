
function [Pmat,Tphi] = stateX_t0_all_y(X,K,N,stoCA,states_all)
% compute the probabilty from state X to all other states
KhatN = K^N; 
Tphi  = zeros(K,N);
Pmat  = zeros(1,KhatN); 
for n=1:N
    nhbr_ind     = stoCA.Ind_nhbr_k(n);
    Xk_nhbr      = X(nhbr_ind);
    nhbr_density = histcounts(Xk_nhbr,stoCA.edges);      % we can put more weights at k later
    nhbr_density = nhbr_density'/length(Xk_nhbr);
    P_n          = stoCA.TMat*nhbr_density;   % a probability vector 1xK for agent n 
    Tphi(:,n)    = P_n';   
end

for l= 1:KhatN
    targetState = states_all(l,:); 
    temp = 1; 
    state_to_ind = targetState;    % the state values are alphabeta; from alphabeta elements to index
    for n=1:N
        k_index = state_to_ind(n); 
        Pn = Tphi(k_index);        % Probability of X_n = targetsate_n
        temp = temp*Pn; 
    end
    Pmat(l) = temp; 
end

end