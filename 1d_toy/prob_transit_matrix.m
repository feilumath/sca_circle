% compute the transition probability matrix of the Markov chain
close all; clear all; % clc

add_mypaths; 

%% setttings
N = 3; % Number of states
K = 2; % Number of possible values for each state
tN = 200; 
nhbrSize = 3; 
stoCA_par = settings_model(K,N,tN,nhbrSize);      % settings of the SCA model: N, K, graph, etc,  


%% compute transition probability matrix P, and invariant measure pi
[Pmat,states_all]  = trans_prob_Mat_markov(stoCA_par,K,N); 

% Display the states
disp(states_all);

figure; imagesc(Pmat); 

[V,D]  = eig(Pmat');              % [V,D] =eig(A) >> A*V= V*D = [D11*v1; D22*v2; ...]; 
                                  % eig is better than svd: 
                                  % [Usvd,S,Vsvd]= svd(Pmat);  >> singular value >1
inv_pi = V(:,1)'/sum(V(:,1));     % A'*v = D11*v >>> v'*A' = D11*v'; 
disp(inv_pi);

%% test dependence on T: adding perturbations to T, and test changes in Pmat and pi
 [diff_Tmat,diff_Pmat]= dependence_on_T(stoCA_par,K,N); 


%%% ============ 

 function [diff_Tmat,diff_Pmat]=dependence_on_T(stoCA,K,N)
% test dependence of P and inv_pi on T

n_test = 100; 
stoCA1 = stoCA; 
Tmat   = stoCA.TMat; 
Pmat   = trans_prob_Mat_markov(stoCA,K,N); 
[V,~]  = eig(Pmat');
inv_pi = V(:,1)'/sum(V(:,1));


diff_Tmat = zeros(1,n_test);
diff_Pmat = zeros(1,n_test);
diff_pi   = zeros(1,n_test);
for m=1:n_test
    T_perturb   = rand(size(Tmat));                                        % to make it stochastic 
    Tmat1 = Tmat+T_perturb;
    Tmat1 = Tmat1*diag(1./sum(Tmat1));
    stoCA1.TMat = Tmat1;
    Pmat1       = trans_prob_Mat_markov(stoCA1,K,N);
    [V1,~]  = eig(Pmat1'); 
    inv_pi1 = V1(:,1)'/sum(V1(:,1));   

    diff_Pmat(m)= norm(Pmat1-Pmat,"fro");
    diff_Tmat(m)= norm(T_perturb,"fro");
    diff_pi(m) = norm(inv_pi1-inv_pi);
end

[diff_Tmat1,ind]    = sort(diff_Tmat);
diff_Pmat1          = diff_Pmat(ind);
diff_pi1            = diff_pi(ind);
figure; 
scatter(diff_Tmat1,diff_Pmat1); hold on; 
scatter(diff_Tmat1,diff_pi1); 
legend('P change','pi change')

figure; 
semilogy(diff_Tmat1,diff_Pmat1./diff_Tmat1,'linewidth',1); hold on; 
semilogy(diff_Tmat1,diff_pi1./diff_Tmat1,'linewidth',2);
legend('|P-diff|/|T-dff|','|pi-diff|/|T-dff|')
xlabel('Frobenus norm T-diff');

figname = 'dependence_on_T'; 
set_positionFontsAll;
end


function [Pmat,states_all]  = trans_prob_Mat_markov(stoCA,K,N)
% Compute the transition probability matrix of the Markov chain. 
% The matrix has a size K^N x K^N. Thus, we can only compute a small one. 

KhatN = K^N;

if KhatN >1024; fprintf('K^N too large to be computed.\n'); return; end

Pmat       = zeros(KhatN,KhatN); 
states_all = generate_all_states(N, K);


states_all = states_all+1; % shift the alphabet from 0:K-1 to 1:K   ***** 

for kn =1:KhatN
    X = states_all(kn,:);  
    [Pmat_1row,~] = stateX_t0_all_y(X,K,N,stoCA,states_all); 
    Pmat(kn,:) = Pmat_1row; 
end

end


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



function states = generate_all_states(N, K)
    % Generate all K^N states for the cellular automata model.
    % Each state is represented as a row in the matrix 'states'.

    % Total number of states
    KhatN = K^N;
    
    % Initialize the matrix to store states
    states = zeros(KhatN, N);
    
    for i = 1:KhatN
        % Compute the ith state
        state = i-1;
        for j = 1:N
            % Convert the state to base K and store it
            states(i, N-j+1) = mod(state, K);
            state = floor(state / K);
        end
    end
end
