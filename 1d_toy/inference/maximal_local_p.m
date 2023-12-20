function [X1_mle_all,prob_X1mle_all] = maximal_local_p(local_p_all)
% find the most likely state (alphabet for all sites) from local densities local_p_all
M = length(local_p_all); 
X1_mle_all     = cell(1,M);
prob_X1mle_all = cell(1,M);
% [K,N,tNminus1] = size(local_p_all{1});
% X1_mle     = zeros(N,tNminus1); 
% prob_X1mle = X1_mle; 
for m =1:M
    [prob_X1mle,X1_mle] = max(local_p_all{m}); 
    X1_mle_all{m}     = squeeze(X1_mle);
    prob_X1mle_all{m} = squeeze(prob_X1mle); 
end

end


