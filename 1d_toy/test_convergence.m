%% test convergence of the estimator -- trajectory data


Mseq  = 10.^(0:5);  
nsimu = 10; 



%% setup 
if ~exist([SAVE_DIR,'figures/'],'dir'), mkdir([SAVE_DIR,'figures/']); end    
figpath0 = [SAVE_DIR,'figures/']; 
estFilename =  [SAVE_DIR,'data_convTest.mat'];

n_Mseq     = length(Mseq);  
lse_seq    = zeros(K,K,n_Mseq,nsimu); 
errLSE_seq = zeros(n_Mseq,nsimu); 


stoCA_par   = settings_model(K,N,tN,nhbrSize);  
inferInfo.M = 2*Mseq(end); 

%%  Generate and assemble a large dataset: to randomly sample from it in n simulations 
tic
% generate/load data: all data
Xt_all = generateData(inferInfo,stoCA_par);   % Xt_all= cell(1,M): each cell is a trajectory of data 
estTime.generateData = toc; 


% Assemble data to local population density phi(X_i,V_i)(t) \in \R^K. It is what needed for Tmat esitmation (LSE or MLE) 
tic
local_p_all = all_local_density(Xt_all,stoCA_par); % cell(1,M): for each traj, (K,N,tN-1); 
estTime.Assemble_densities = toc; 


%% Part 1: LSE using trajectory information 
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
save(estFilename,"estTime","errLSE_seq","stoCA_par","Mseq","inferInfo");

%% present the results
figure
boxplot(errLSE_seq');
ax = gca; ax.YAxis.Scale ="log";
set(gca, 'XTickLabels',compose('%d',Mseq)); 
xlabel('Sample size M'); ylabel('Relative Error'); 
% title(compose('Relative error in %d simulations',n)); 
figname = [figpath0,'fig_LSE_traj'];
set_positionFontsAll;

%% Part 2: LSE without using trajectory information 
fprintf('\n Estimator from ensemble data without trajectory information\n');

lse_seq_EnsD    = zeros(K,K,n_Mseq,nsimu); 
errLSE_seq_EnsD = zeros(n_Mseq,nsimu); 
tic
for n= 1:nsimu
M_index = randperm(inferInfo.M);
for i=1:n_Mseq
    M = Mseq(i);
    data_indx = M_index(1:M);  
    Xm_all        = data_Xt2Xm(Xt_all(data_indx));         % data Xt_all = cell(1,tN); each time can has M(t) samples
    local_p_all_M = data_pt2pm(local_p_all(data_indx));    % data local_p_all_M = cell(1,tN); 
    [Tmat_lse,~] = infer_from_sitesPDF(Xm_all,local_p_all_M,K);  % 
    lse_seq_EnsD(:,:,i,n) = Tmat_lse; 
    errLSE_seq_EnsD(i,n)  = norm(Tmat_lse - stoCA_par.TMat,'fro');
end
errLSE_seq_EnsD = errLSE_seq_EnsD/norm(stoCA_par.TMat,'fro'); 
end
estTime.LSE_EnsD = toc;

save(estFilename,"estTime","errLSE_seq_EnsD","lse_seq_EnsD",'-append');

%% present the results
figure
boxplot(errLSE_seq_EnsD');
ax = gca; ax.YAxis.Scale ="log";
set(gca, 'XTickLabels',compose('%d',Mseq)); 
xlabel('Sample size M'); ylabel('Relative Error');  
figname = [figpath0,'fig_LSE_Ensemble'];
set_positionFontsAll;