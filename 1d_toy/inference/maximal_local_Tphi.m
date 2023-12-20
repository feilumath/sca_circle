function [X1_mle_all,prob_X1mle_all] = maximal_local_Tphi(local_p_all,Tmat)
% find the most likely state (alphabet for all sites) from local densities: Tmat*local_p_all
M = length(local_p_all); 
X1_mle_all     = cell(1,M);
prob_X1mle_all = cell(1,M);
[K,N,tN] = size(local_p_all{1});
% X1_mle     = zeros(N,tNminus1); 
% prob_X1mle = X1_mle; 

Tmat = reshape(Tmat,[K,K]);
for m =1:M
    temp   = local_p_all{m}; 
    [~,N,tN] = size(temp); 
    for t = 1:tN
        temp(:,:,t) = Tmat*temp(:,:,t); 
    end
    [prob_X1mle,X1_mle] = max(temp);  
    X1_mle_all{m}     = squeeze(X1_mle);
    prob_X1mle_all{m} = squeeze(prob_X1mle); 
end

end


