function [Tmat_mle,fval,exitflag,output] = infer_LSE_MLE(local_p_all,Xt_all,stoCA,inferInfo)
% Estimate the tranform stochastic matrix by LSE and MLE 

M = length(local_p_all); 
K = stoCA.K; 

%% LSE estimator min (1- P(T)): phi= local_p_all{m}    || 1_{k} - T* phi ||^2 >> phi*phi' T = phi; phi KxN for each t
Amat = zeros(K,K);  % phi*phi'   with phi KxN for each t 
bvec = zeros(K,K);  % phi*1_[Xt1](k)             1_[Xt1] size 1xK   >>> lse is column stochastic: sum_j T(k,j) =1
lse_all = cell(1,M);
for m= 1:M
    Amat1traj = zeros(K,K); bvec1traj = zeros(K,K);
    phi1traj  = local_p_all{m}; 
    Xt        = Xt_all{m};  % Xt size = [N,tN] 
    [K1,~,tN] = size(phi1traj); 
    if K1 ~= K; fprintf('\n local_p_all size miss match\n'); end
    for t=1:tN-1
        Xt1       = Xt(:,t+1); % size Nx1, with value in [1,K]; 
        N = length(Xt1); 
        indicator = zeros(N,K); 
        for n= 1:N
            indicator(n,Xt1(n)) = 1; 
        end 
        Tmat_lse = phi1traj(:,:,t);   % KxN
        Amat1traj = Amat1traj+ Tmat_lse*Tmat_lse'; 
        bvec1traj = bvec1traj+ Tmat_lse*indicator;         
    end
    lse_all{m}.Amat = Amat1traj; 
    lse_all{m}.bvec = bvec1traj;
    Amat = Amat +  Amat1traj; 
    bvec = bvec + bvec1traj;
end
% lse.cvec = Amat\bvec; 
lb = 0; ub = 1; 
lse.cvec = lsqlin(Amat,bvec,[],[],[],[],lb,ub); 
lse.Amat = Amat; 
lse.bvec = bvec; 
Tmat_lse = lse.cvec';


%% MLE estimator 
%  T_matIC  = inferInfo.T_matIC; 
 T_matIC  = reshape(Tmat_lse',[K^2,1]);
% optimization with respect to T_mat to minimize the likelihood: TBD
negLog_lkhd_Fn = @(Tmat) -1*likelihood_Mtraj(local_p_all,Xt_all,stoCA,Tmat,M);  % we can use random batch to make the optimization faster

fprintf('\n -Log likelihood at IC = %4.4f \n',negLog_lkhd_Fn(T_matIC) ) % tested good, 

fprintf('\n Optimization by fmincon. \n')
options         = optimoptions('fmincon');
options.Display = 'iter';
A = []; B=[]; Aeq =inferInfo.Aeq; beq = inferInfo.beq; LB= inferInfo.lb; UB = inferInfo.ub; NONLCON=[]; 
[Tmat_mle,fval,exitflag,output ]= fmincon(negLog_lkhd_Fn,T_matIC,A,B,Aeq,beq,LB,UB,NONLCON,options);

% display results: 
Tmat_mle = reshape(Tmat_mle,[K,K])   % Tmat_est = Tmat_est*diag(1./sum(Tmat_est))   % normalize to stochastic
Tmat_true= stoCA.TMat 
Tmat_lse = lse.cvec' 


 %% plot the log-likelihood function when K==2; 
if K==2
    x1 = 0.8:0.01:1; x2= 0.03:0.01:.2;
    n1 = length(x1); n2 = length(x2);
    lkhd_val = zeros(n1,n2);
    for i= 1:n1
        for j=1:n2
            Tmat = [x1(i);x2(j); 1-x1(i);1-x2(j)];
            lkhd_val(i,j) =  negLog_lkhd_Fn(Tmat);
        end
    end
    
    [X,Y] = meshgrid(x2,x1);
    figure
    surf(X,Y,lkhd_val);xlabel('x2'); ylabel('x1'); title('Liklihood function');
    figure;
    imagesc(x1,x2,lkhd_val); xlabel('x1'); ylabel('x2'); colorbar; title('Liklihood values');
    print('fig_lkhd.pdf','-dpdf')
end
end


function log_lkhd = likelihood_Mtraj(local_p_all,Xt_all,stoCA,T_mat,M)
% compute the likelihood of M trajectories

if ~isequal( size(T_mat), [stoCA.K,stoCA.K] )     % if the two array are of different size 
   T_mat = reshape(T_mat,[stoCA.K,stoCA.K]);
end

stoCA.TMat = T_mat;  % for iteration in estimation 
log_lkhd = 0; 
for m=1:M
    Xt= Xt_all{m}; phi1traj  = local_p_all{m}; 
    log_lkhd_1traj = likelihood_1tra(Xt,phi1traj,T_mat);   % compute the likelihood of 1-trajectory Xt
    log_lkhd       = log_lkhd+ log_lkhd_1traj;
end

end


function log_lkhd = likelihood_1tra(Xt,phi1traj,T_mat)
% compute the likelihood of 1-trajectory Xt
[N,tN] = size(Xt);
log_lkhd = 0;   % the likelihood is very small, use log10
for t=1:tN-1
    X1 = Xt(:,t+1); % N x 1, with values in [1,K]
    temp       = T_mat*phi1traj(:,:,t);  % K x N 
    lkhd_1step = zeros(1,N);
    for n=1:N
        lkhd_1step(n) = temp(X1(n),n);
    end
    log_lkhd       = log_lkhd+ sum(log10(lkhd_1step)); 
end

end