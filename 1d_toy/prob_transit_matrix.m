% compute the transition probability matrix of the Markov chain
close all; clear all; % clc
add_mypaths; 

%% setttings
N = 3; % Number of states
K = 3; % Number of possible values for each state
tN = 200; 
nhbrSize = min(3,floor(N/2));  %% in both left and right;   neighbor size = 2*nhbrSize+1; 
stoCA_par = settings_model(K,N,tN,nhbrSize);      % settings of the SCA model: N, K, graph, etc,  


%% compute transition probability matrix P, and invariant measure pi
[Pmat,states_all]  = trans_prob_Mat_markov(stoCA_par,K,N); 
figure; imagesc(Pmat); 

[V,D]  = eig(Pmat');              % [V,D] =eig(A) >> A*V= V*D = [D11*v1; D22*v2; ...]; 
                                  % eig is better than svd: 
                                  % [Usvd,S,Vsvd]= svd(Pmat);  >> singular value >1
inv_pi = V(:,1)'/sum(V(:,1));     % A'*v = D11*v >>> v'*A' = D11*v'; 

% Display the states
if K^N<=64 
    disp(states_all); 
    disp(inv_pi); 
end

%% test dependence on T: adding perturbations to T, and test changes in Pmat and pi
nsimu  = 100; 
plotON = 1; 
randT  = 0;    % when randT =1, the eigenvector-based pi may have machine precision error.
[diff_Tmat,diff_Pmat,diff_pi,diff_l1]= dependence_on_T(stoCA_par,K,N,nsimu,plotON,randT); 
close all; 

%% a sequence of N: 
% Nseq = 3:9; K    = 2; 
 Nseq = 3:7; K    = 3; 

datafilename = [SAVE_DIR,'dependence_on_T_varyN_K',sprintf('%i_l1_randT%i.mat',K,randT)];
if ~exist(datafilename,'file')
    n_Nseq      = length(Nseq);
    diff_Tseq   = zeros(n_Nseq,nsimu);
    diff_Pseq   = zeros(n_Nseq,nsimu);
    diff_pi_seq = zeros(n_Nseq,nsimu);


    mean_ratio_PT_l1  = zeros(n_Nseq,1);
    mean_ratio_piT_l1 = zeros(n_Nseq,1);
    for nn= 1:n_Nseq
        N  = Nseq(nn);
        nhbrSize = min(3,floor(N/2));  %% in both left and right;   neighbor size = 2*nhbrSize+1;
        stoCA_par = settings_model(K,N,tN,nhbrSize);      % settings of the SCA model: N, K, graph, etc,
        [diff_Tseq(nn,:),diff_Pseq(nn,:),diff_pi_seq(nn,:),diff_l1]= dependence_on_T(stoCA_par,K,N,nsimu,0,randT);

        ratio_PT_l1       = diff_l1.diff_Pmat_l1./diff_l1.diff_Tmat_l1;
        ratio_piT_l1      = diff_l1.diff_pi_l1./diff_l1.diff_Tmat_l1;
        mean_ratio_PT_l1(nn)  = mean(ratio_PT_l1);
        mean_ratio_piT_l1(nn) = mean(ratio_piT_l1);
    end

    ratio_PT       = diff_Pseq./diff_Tseq;
    ratio_piT      = diff_pi_seq./diff_Tseq;
    mean_ratio_PT  = mean(ratio_PT,2);
    mean_ratio_piT = mean(ratio_piT,2);
    
    save(datafilename,'Nseq','mean_ratio_PT','mean_ratio_piT', "K","diff_pi_seq","diff_Tseq","diff_Pseq",'nsimu',...
        'mean_ratio_PT_l1',"mean_ratio_piT_l1");
end

load(datafilename); 

figure; 
subplot(1,2,1)
semilogy(Nseq,mean_ratio_PT,'--x','linewidth',2); hold on
semilogy(Nseq,mean_ratio_piT,'-.d','linewidth',2); 
semilogy(Nseq,mean_ratio_piT(1)*K.^(-( Nseq-Nseq(1)) /2),':','linewidth',2);
semilogy(Nseq,mean_ratio_PT_l1(1)*K.^(( Nseq-Nseq(1))),':','linewidth',2);
% legend('|P_1-P_2|_2/|T_1-T_2|_2','|\pi_1-\pi_2|_2/|T_1-T_2|_2','c_0 K^{-N/2}')
legend('|\Delta P|_2 /  |\Delta T |_2','|\Delta \pi| _2/ |\Delta T|_2','c_0 K^{-N/2}','');
legend('location','best')
xlabel(['N (with K=',sprintf('%i)',K)] );
ylabel('Mean in 100 simulations');

subplot(1,2,2)
semilogy(Nseq,mean_ratio_PT_l1,'--x','linewidth',2); hold on
semilogy(Nseq,mean_ratio_piT_l1,'-.d','linewidth',2); 
semilogy(Nseq,mean_ratio_piT_l1(1)*K.^(-( Nseq-Nseq(1)) /2),':','linewidth',2);
semilogy(Nseq,mean_ratio_PT_l1(1)*K.^(( Nseq-Nseq(1))),':','linewidth',2);
legend('|\Delta P|_1 /  |\Delta T |_1','|\Delta \pi| _1/ |\Delta T|_1','','c_0 K^{N}')
legend('location','best'); % legend('boxoff')
xlabel(['N (with K=',sprintf('%i)',K)] );
% ylabel('Mean in 100 simulations');


figname = [SAVE_DIR,'mean_ratio_varyN_K',sprintf('%i_l1_randT%i',K,randT)];
set_positionFontsAll;



%%% ============ 




