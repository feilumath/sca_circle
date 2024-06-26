% script to run the code 
% Stochastic cellular Automata
% Copyright: Mauro Maggioni, Sichen Yang, Fei Lu
% Last edited: by Fei Lu, 2023-12-15 %2021-10-26

clear all; close all; clc; 
add_mypaths; 

testErgodicityON = 1; 

est_predict_synchr = 1;  

%% setttings
K   = 2;           % size of the alphabet set
N   = 6;   
nhbrSize = 3;  % number of sites and neighbor size
tN  = 100;          % number of time steps

stoCA_par = settings_model(K,N,tN,nhbrSize);              % settings of the SCA model: N, K, graph, etc,  
  
%% test the StoCA model 
Xt = stoCA_model(stoCA_par);
tN =  stoCA_par.tN; gap  = ceil(tN/100); 
tInd = 1:gap:tN;
figure; 
imagesc(Xt(:,tInd)); xlabel('Time'); ylabel('Sites');


%% test empirical distribution, ergodicity, and Consensus/clustering
if testErgodicityON
    test_ergodicity(stoCA_par);   % path average converges fast. Ergodic? Consensus?
end



%% Inference-- 
% % part I: multiple trajectory data
M  = 1e4; 
inferInfo = settings_infer(M,K); 

% generate/load data: the data size may be large, to make it efficient later 
Xt_all = generateData(inferInfo,stoCA_par);   % Xt_all= cell(1,M): each cell is a trajectory of data 

% Assemble data to local population density phi(X_i,V_i)(t) \in \R^K. It is what needed for Tmat esitmation (LSE or MLE) 
local_p_all = all_local_density(Xt_all,stoCA_par); % cell: for each traj, (K,N,tN-1); 

% Estimate the parameter by LSE, and LSE-stocMat: the two results are almost same
% They works well: high probability events happens the most; They are not the MLE. 
lse_stocON = 1; % if 1: compute lse_stoc that estimates K-1 rows
[lse_stoc,lse_all,Tmat_lse] = LSE_stocMat(local_p_all,Xt_all,stoCA_par,lse_stocON);
Tmat_lse                       % LSE without using T is column-stochastic 
lse_stoc                       % LSE of the only the first (K-1) rows
Tmat_true = stoCA_par.TMat 
norm(Tmat_lse- stoCA_par.TMat,'fro')

% Inference-- part II:  ensemble data without trajectory information: use maximal likely local densities 
% function infer_Ensemble_noTraj(local_p_all,Xt_all,stoCA_par,inferInfo)
% estimate Tmat when data are ensemble without trajectory information
% Idea: match the marginal density of each site. Need large ensemble size

Xm_all        = data_Xt2Xm(Xt_all);         % data Xt_all = cell(1,tN); each time can has M(t) samples
local_p_all_M = data_pt2pm(local_p_all);    % data local_p_all_M = cell(1,tN); 

[Tmat_lse,lse_struct] = infer_from_sitesPDF(Xm_all,local_p_all_M,stoCA_par.K);  % 
norm(Tmat_lse- stoCA_par.TMat,'fro')



%% test convergence for both estimators as Sample size increases
test_convergence;  


%% Inference-Predict synchronization
infer_predict_synchr; 

%% addition comments: 
%{
%% MLE: estimate the transform matrix Tmat by Maximal likelihood
% -- it re-estimates the local density for each optimization iteraction. Inefficient. Us the LSE
% [Tmat_est,fval,~,output] = infer_MLE(Xt_all,stoCA_par,inferInfo);   
% [Tmat_est,fval,~,output] = infer_LSE_MLE(local_p_all,Xt_all,stoCA_par,inferInfo);  
%}
