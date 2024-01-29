
function [Pmat,states_all]  = trans_prob_Mat_markov(stoCA,K,N)
% Compute the transition probability matrix of the Markov chain. 
% The matrix has a size K^N x K^N. Thus, we can only compute a small one. 

KhatN = K^N;

if KhatN >4*1024; fprintf('\n  Warning: K^N too large to be computed.\n'); return; end

Pmat       = zeros(KhatN,KhatN); 
states_all = generate_all_states(N, K);


states_all = states_all+1; % shift the alphabet from 0:K-1 to 1:K   ***** 

for kn =1:KhatN
    X = states_all(kn,:);  
    [Pmat_1row,~] = stateX_t0_all_y(X,K,N,stoCA,states_all); 
    Pmat(kn,:) = Pmat_1row; 
end

end
