%% test convergence of the estimator -- trajectory data


Mseq  = 10.^(0:4);  
nsimu = 10; 

n_Mseq     = length(Mseq);  
lse_seq    = zeros(K,K,n_Mseq,nsimu); 
errLSE_seq = zeros(n_Mseq,nsimu); 


settings_infer; 
inferInfo.M = 2*Mseq(end); 

tic
% generate/load data: all data
Xt_all = generateData(inferInfo,stoCA_par);   % Xt_all= cell(1,M): each cell is a trajectory of data 
estTime.generateData = toc; 


tic
% Assemble data to local population density phi(X_i,V_i)(t) \in \R^K. It is what needed for Tmat esitmation (LSE or MLE) 
local_p_all = all_local_density(Xt_all,stoCA_par); % cell(1,M): for each traj, (K,N,tN-1); 
estTime.Assemble_densities = toc; 


%% Estimate the parameter by LSE:
tic
lse_stocON = 0; % if 1: compute lse_stoc that estimates K-1 rows

for n= 1:nsimu
M_index = randperm(inferInfo.M);
for i=1:n_Mseq
    M = Mseq(i);
    data_indx = M_index(1:M);  
    [~,~,Tmat_lse] = LSE_stocMat(local_p_all(data_indx),Xt_all(data_indx),stoCA_par,lse_stocON);
    lse_seq(:,:,i,n) = Tmat_lse; 
    errLSE_seq(i,n)  = norm(Tmat_lse - stoCA_par.TMat,'fro');
end
errLSE_seq = errLSE_seq/norm(stoCA_par.TMat,'fro'); 
end

estTime.LSE = toc; 

%% present the results

figure
boxplot(errLSE_seq');
ax = gca; ax.YAxis.Scale ="log";
set(gca, 'XTickLabels',compose('%d',Mseq)); 
xlabel('Sample size M'); ylabel('Relative Error'); 
title(compose('Relative error in %d simulations',n)); 

