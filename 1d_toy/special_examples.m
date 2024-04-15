% special examples
% 

% clear all; close all; clc; 
add_mypaths; 

% Define a custom colormap (e.g., RGB triplets for red, green, and blue)
% Colors: Red, Green, Blue
myColormap2 = [255 0 0;   % Red, values between 0 and 1
              255,255,0; % Yellow, % 46,139,87;   % sea green 
              0,128,0]/255; % Green            % 30,144,255]/255;  % dodger blue

myColormap = [173,216,230;   % Red, values between 0 and 1
              255,255,0; % Yellow, % 46,139,87;   % sea green 
              0,206,209]/255; % Green            % 30,144,255]/255;  % dodger blue

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
% TMat = TMat*diag(1./sum(TMat)); % Column sum is 1. 
stoCA_par.TMat = Tmat; 

X0 = randi(K,N,1); % ones(N,1)*1; 
Xt = stoCA_model(stoCA_par,X0);
tN =  stoCA_par.tN; gap  = ceil(tN/100); 
tInd = 1:gap:K*6;
figure; 
imagesc(Xt(:,tInd)); xlabel('Time'); ylabel('Sites');
colormap("parula"); default
% colormap(myColormap);

figname = [figpath0,'Example_synchronization',str_name];
set_positionFontsAll;

eig(Tmat)


%% T from deterministic to stochastic
K=3; 
Tmat = eye(K); 
Tmat = [Tmat(2:K,:); ones(1,K)/K]; 
stoCA_par.TMat = Tmat; 

X0 = ones(N,1)*1; 
Xt = stoCA_model(stoCA_par,X0);
tN =  stoCA_par.tN; gap  = ceil(tN/100); 
tInd = 1:K*6;
figure; 
imagesc(Xt(:,tInd)); xlabel('Time'); ylabel('Sites');
str_name = sprintf('N%i_K%i',N,K);
figname = [figpath0,'Example_det_stoc',str_name];
set_positionFontsAll;



%% T periodic, clusters >>> synchronization 
K = 4; 
stoCA_par = settings_model(K,N,tN,nhbrSize);    

% Tmat = [0,1,0,0; 1,0, 0,0; 0, 0, 0, 1; 0, 0, 1,0]; 
% Tmat = [0,0,0.5,0.5; 0,0, 0.5,0.5; 0.5,0.5, 0, 0; 0.5,0.5, 0,0]; 

% Tmat = [0,0,0,1; 0,0, 1,0; 0,1, 0, 0; 1,0, 0,0]; 
Tmat = [0,1,0,0; 0,0, 0.5,0.5; 0, 0, 0, 1; 1, 0, 0,0]; 
% TMat = TMat*diag(1./sum(TMat)); % Column sum is 1. 
stoCA_par.TMat = Tmat; 

X0 =   randi(K,N,1);  % ones(N,1)*1; %
Xt = stoCA_model(stoCA_par,X0);
tN =  stoCA_par.tN; gap  = ceil(tN/100); 
tInd = 1:gap:tN; % K*10;
figure; 
imagesc(Xt(:,tInd)); xlabel('Time'); ylabel('Sites');
str_name = sprintf('N%i_K%i',N,K);
figname = [figpath0,'Example_synchro_clusters',str_name];
set_positionFontsAll;
eig(Tmat)