function [Tmat_est,fval,exitflag,output] = infer_MLE(Xt_all,stoCA,inferInfo)
% Estimate the tranform stochastic matrix by maximal likelihood
%{
to be improved:
  - provide derivative -- it can be computed directly. 
    Not sure if it will actually help. The gradient of the loss function is always nonzero, and there are constraints. 
  - random batch for large M 
%} 

M = inferInfo.M;  
K = stoCA.K; 
T_matIC  = inferInfo.T_matIC; 
%% optimization with respect to T_mat to minimize the likelihood: TBD
negLog_lkhd_Fn = @(Tmat) -1*likelihood_Mtraj(Xt_all,stoCA,Tmat,M);  % we can use random batch to make the optimization faster

fprintf('\n -Log likelihood at IC = %4.4f \n',negLog_lkhd_Fn(T_matIC) ) % tested good, 

fprintf('\n Optimization by fmincon with all data. To be improved by mini-batch/SGD \n')
A = []; B=[]; Aeq =inferInfo.Aeq; beq = inferInfo.beq; LB= inferInfo.lb; UB = inferInfo.ub; NONLCON=[]; 
[Tmat_est,fval,exitflag,output ]= fmincon(negLog_lkhd_Fn,T_matIC,A,B,Aeq,beq,LB,UB,NONLCON);

% display results: 
Tmat_est = reshape(Tmat_est,[K,K])   % Tmat_est = Tmat_est*diag(1./sum(Tmat_est))   % normalize to stochastic
Tmat_true= stoCA.TMat 


end


function log_lkhd = likelihood_Mtraj(Xt_all,stoCA,T_mat,M)
% compute the likelihood of M trajectories

T_mat = reshape(T_mat,[stoCA.K,stoCA.K]);
stoCA.TMat = T_mat;  % for iteration in estimation 
log_lkhd = 0; 
for m=1:M
    Xt= Xt_all{m}; 
    log_lkhd_1traj = likelihood_1tra(Xt,stoCA);   % compute the likelihood of 1-trajectory Xt
    log_lkhd       = log_lkhd+ log_lkhd_1traj;
end

end


function log_lkhd = likelihood_1tra(Xt,stoCA)
% compute the likelihood of 1-trajectory Xt
N    = stoCA.N;
tN   = stoCA.tN; 

log_lkhd = 0;   % the likelihood is very small, use log10
for t=1:tN
    X0         = Xt(:,t); X1 = Xt(:,t+1);
    lkhd_1step = move1step_infer(X0,X1,stoCA);   
    log_lkhd       = log_lkhd+ log10(lkhd_1step); 
end

end