% settings of inference 

M = 1000;    % number of trajectories 
K = stoCA_par.K; 


% initialize the optimization
T_mat   = rand(K,K);  T_mat = T_mat*diag(1./sum(T_mat)); 
T_matIC = reshape(T_mat,[K^2,1]);

Aeq = zeros(K, K^2);
for i=1:K
    indA1        = (i-1)*K + (1:K);    
    Aeq(i,indA1) = 1; 
end
beq = ones(K,1); 
lb  = zeros(K^2,1); 
ub  = 1+lb; 


inferInfo  = struct(...
   'method', 'MLE',...     %  'LSE_MLE'    
   'M',       M, ...            % number of trajectories 
   'T_matIC', T_matIC,...
   'Aeq',     Aeq,...
   'beq',     beq,...
   'lb',      lb,...
   'ub',      ub,...
   'dataType', 'synthetic'...    
   ); 

