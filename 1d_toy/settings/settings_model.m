

% settings of the stochastic cellular automata model
K   = 2;           % size of the alphabet set
tN  = 100;          % number of time steps
alphabet = 1:K; 
N   = 3;   nhbrSize = 3;  % number of sites and neighbor size


P0       = 1+ sin(2*(1:K));   % Initial distribution
IC_distr = P0/sum(P0);       

%% transform matrix (to be learned from data)
TMat = rand(K,K);  TMat = TMat*diag(1./sum(TMat)); % Column sum is 1. 

% neighborhood of site k 
siteInd_ext = [N-nhbrSize+1:N,1:N,1:nhbrSize];                     % cyclic site indices --- on a circle
Ind_nhbr_k  = @(k)siteInd_ext( (-nhbrSize:nhbrSize)+k + nhbrSize); % index of neighborhood of site k


figure; 
imagesc(TMat); colorbar; 
title('Transformation matrix T KxK')

%% all parameters
stoCA_par = struct(...
    'alphabet',    alphabet,...   % alphabet set; 
    'K',           K,...          % size of the alphabet set
    'N',           N, ...       % number of sites
    'tN',          tN, ...        % number of time steps
    'nhbrSize',    nhbrSize,...   % neighbor size
    'siteInd_ext', siteInd_ext,...%
    'Ind_nhbr_k',  Ind_nhbr_k,... 
    'TMat',       TMat,...        % transform matrix 
    'edges',      [0,1:K]+.5,...  % edges for histcounts in nbhr
    'IC_distr',    IC_distr,...   % initial distribution
    'name',       'Toy1D'... 
);




%{
% test random number generator for each move 
X = 1:3;
P = 1+sin(X); P=P/sum(P);
N=1000; T=zeros(1,N); for i=1:N; T(i) = gendist(P,1,1); end
P_est = histcounts(T,3)/N; 

N=1000; T=zeros(N,length(P)); for i=1:N; T(i,:) = mnrnd(1,P,1); end
P_estMnd = sum(T,1)/N
%}



