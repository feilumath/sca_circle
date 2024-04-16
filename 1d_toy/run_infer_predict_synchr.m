% use estimator to predict synchr

% special examples
% 

 clear all; close all; clc; 
add_mypaths; 


% setttings
K   = 3;             % size of the alphabet set
N   = 8;   
nhbrSize = 3;       % number of sites and neighbor size in each side
tN  = 100;          % number of time steps

stoCA_par = settings_model(K,N,tN,nhbrSize);              % settings of the SCA model: N, K, graph, etc,  
  
if ~exist([SAVE_DIR,'figures/'],'dir'), mkdir([SAVE_DIR,'figures/']); end
figpath0 = [SAVE_DIR,'figures/'];
str_name = sprintf('N%i_K%i',N,K);

%% T periodic: permutation >>> synchronization 
Tmat = eye(K); 
Tmat = [Tmat(2:K,:); Tmat(1,:)];  
Tmat = diag(1./sum(Tmat,2))*Tmat; % row sum is 1. 
stoCA_par.TMat = Tmat; 


% demo of a trajectory
X0 = randi(K,N,1); % ones(N,1)*1; 
Xt_true = stoCA_model(stoCA_par,X0);
tN =  stoCA_par.tN; gap  = ceil(tN/100); 
tInd = 1:gap:K*6;
figure(100); 
subplot(131); imagesc(tInd,1:N, Xt_true(:,tInd)); xlabel('Time'); ylabel('Sites');
% colormap("parula"); % default


%% Estimator
% part I: multiple trajectory data
M  = 1e2; 
inferInfo = settings_infer(M,K); 

% generate/load data: the data size may be large, to make it efficient later 
Xt_all = generateData(inferInfo,stoCA_par);   % Xt_all= cell(1,M): each cell is a trajectory of data 

% Assemble data to local population density phi(X_i,V_i)(t) \in \R^K. It is what needed for Tmat esitmation (LSE or MLE) 
local_p_all = all_local_density(Xt_all,stoCA_par); % cell: for each traj, (K,N,tN-1); 

lse_stocON = 1; % if 1: compute lse_stoc that estimates K-1 rows
[lse_stoc,~,Tmat_lse] = LSE_stocMat(local_p_all,Xt_all,stoCA_par,lse_stocON);

% Inference-- part II:  ensemble data without trajectory information: use maximal likely local densities 
Xm_all        = data_Xt2Xm(Xt_all);         % data Xt_all = cell(1,tN); each time can has M(t) samples
local_p_all_M = data_pt2pm(local_p_all);    % data local_p_all_M = cell(1,tN); 

[Tmat_lse_Ens,~] = infer_from_sitesPDF(Xm_all,local_p_all_M,stoCA_par.K);  % 

Tmat_true = stoCA_par.TMat; 
norm(Tmat_lse - stoCA_par.TMat,'fro')
norm(Tmat_lse_Ens- stoCA_par.TMat,'fro')



%% predict synchronization
stoCA_par_LSE = stoCA_par; 

% TMat = TMat*diag(1./sum(TMat)); % Column sum is 1. 

stoCA_par_LSE.TMat = Tmat_lse; 
Xt_lse = stoCA_model(stoCA_par_LSE,X0); 

    
stoCA_par_lse_Ens = stoCA_par; 
stoCA_par_lse_Ens.TMat = Tmat_lse_Ens; 
Xt_lse_Ens= stoCA_model(stoCA_par_lse_Ens,X0); 

figure(100); clf
subplot(131); imagesc(Xt_true(:,tInd)); xlabel('Time'); title('True'); ylabel('Sites');
subplot(132);  imagesc(Xt_lse(:,tInd)); xlabel('Time'); title('LSE-trajectory') % ylabel('Sites');
subplot(133);  imagesc(Xt_lse_Ens(:,tInd)); xlabel('Time'); title('LSE-ensemble') % ylabel('Sites');


figname = [figpath0,'predict_synchronization',str_name];
set_positionFontsAll;

% for latex 
lse1 = Tmat_lse; 
fprintf('\n LSE-traj = %2.4f & %2.4f & %2.4f \\  %2.4f & %2.4f & %2.4f \\  %2.4f & %2.4f & %2.4f \n', ...
       lse1(1,1), lse1(1,2), lse1(1,3), lse1(2,1),lse1(2,2), lse1(2,3), lse1(3,1),lse1(3,2), lse1(3,3) );

lse1 = Tmat_lse_Ens; 
fprintf('\n LSE-ensemble = %2.4f & %2.4f & %2.4f \\  %2.4f & %2.4f & %2.4f \\  %2.4f & %2.4f & %2.4f \n', ...
       lse1(1,1), lse1(1,2), lse1(1,3), lse1(2,1),lse1(2,2), lse1(2,3), lse1(3,1),lse1(3,2), lse1(3,3) );




