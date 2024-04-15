%% test convergence of the estimator -- trajectory data


Mseq  = 10.^(0:5);
nsimu = 10;


N =6; K=2; 
randomT = 1; % random T in nsimu simulations ---- longer time in generating data. 


%% setup
if ~exist([SAVE_DIR,'figures/'],'dir'), mkdir([SAVE_DIR,'figures/']); end
figpath0 = [SAVE_DIR,'figures/'];
str_name = sprintf('N%i_K%i_M1e%i_nsimu%i_randT%i',N,K,log10(Mseq(end)),nsimu,randomT);
estFilename =  [SAVE_DIR,'data_convTest_',str_name,'.mat'];

n_Mseq     = length(Mseq);


stoCA_par   = settings_model(K,N,tN,nhbrSize);
inferInfo.M = 5*Mseq(end);


if  ~exist(estFilename,'file')

    if randomT == 0
        %  Generate and assemble a large dataset: to randomly sample from it in n simulations
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
        lse_seq    = zeros(K,K,n_Mseq,nsimu);
        errLSE_seq = zeros(n_Mseq,nsimu);
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
    else      % random T
        lse_seq    = zeros(K,K,n_Mseq,nsimu);
        errLSE_seq = zeros(n_Mseq,nsimu);

        lse_seq_EnsD    = zeros(K,K,n_Mseq,nsimu);
        errLSE_seq_EnsD = zeros(n_Mseq,nsimu);

        inferInfo.M = 1*Mseq(end);
        tic 
   
        parfor n= 1:nsimu     % par for 
            TMat = rand(K,K);  TMat = TMat*diag(1./sum(TMat)); % Column sum is 1.
            stoCA_par1      = stoCA_par;
            stoCA_par1.TMat = TMat;
            % generate/load data: all data
            Xt_all      = generateData(inferInfo,stoCA_par1);        % Xt_all= cell(1,M): each cell is a trajectory of data
            local_p_all = all_local_density(Xt_all,stoCA_par1); % cell(1,M): for each traj, (K,N,tN-1);
            lse_stocON  = 0; % if 1: compute lse_stoc that estimates K-1 rows

           
            % LSE-traj
            fprintf('\n LSE-trajectory %i\n',n);
            tic 
            tempLSE = zero(n_Mseq,1); 
            % M_index = 1:inferInfo.M; % randperm(inferInfo.M);
            for i=1:n_Mseq
                M = Mseq(i);
                data_indx        = 1:M;
                [~,~,Tmat_lse]   = LSE_stocMat(local_p_all(data_indx),Xt_all(data_indx),stoCA_par1,lse_stocON);
                lse_seq(:,:,i,n) = Tmat_lse;
                tempLSE(i) = norm(Tmat_lse - stoCA_par1.TMat,'fro');
            end
            errLSE_seq(:,n)= tempLSE/norm(stoCA_par1.TMat,'fro'); 

            % Ensemble-LSE 
            fprintf('\n Ensemble-LSE \n');
            M_index = 1:inferInfo.M; % randperm(inferInfo.M);
            tempLSE_ens = zero(n_Mseq,1); 
            for i=1:n_Mseq
                M             = Mseq(i);
                data_indx     = M_index(1:M);
                Xm_all        = data_Xt2Xm(Xt_all(data_indx));         % data Xt_all = cell(1,tN); each time can has M(t) samples
                local_p_all_M = data_pt2pm(local_p_all(data_indx));    % data local_p_all_M = cell(1,tN);
                [Tmat_lse,~]  = infer_from_sitesPDF(Xm_all,local_p_all_M,K);  %
                lse_seq_EnsD(:,:,i,n) = Tmat_lse;
                tempLSE_ens(i)  = norm(Tmat_lse - stoCA_par1.TMat,'fro');
            end
            errLSE_seq_EnsD(:,n) = tempLSE_ens/norm(stoCA_par1.TMat,'fro');
        end      

       estTime.LSE_all = toc;
       save(estFilename,"estTime","errLSE_seq","stoCA_par","Mseq","inferInfo","errLSE_seq_EnsD","lse_seq_EnsD");
    end
end


%% present the results: 
% % with trajectory 
load(estFilename); 
figure
boxplot(errLSE_seq');
ax = gca; ax.YAxis.Scale ="log";
set(gca, 'XTickLabels',compose('%d',Mseq));
xlabel('Sample size M'); ylabel('Relative Error');
% title(compose('Relative error in %d simulations',n));
figname = [figpath0,'fig_LSE_traj_',str_name];
set_positionFontsAll;

% % without trajectory 
figure
boxplot(errLSE_seq_EnsD');
ax = gca; ax.YAxis.Scale ="log";
set(gca, 'XTickLabels',compose('%d',Mseq));
xlabel('Sample size M'); ylabel('Relative Error');
figname = [figpath0,'fig_LSE_Ensemble_',str_name];
set_positionFontsAll;

% two LSEs in one: 
figure; 
ylim_lb = min(errLSE_seq,[],'all');
ylim_ub = max(errLSE_seq_EnsD,[],'all');

subplot(121); boxplot(errLSE_seq');
ax = gca; ax.YAxis.Scale ="log"; ylim([ylim_lb,ylim_ub]); 
set(gca, 'XTickLabels',compose('%d',Mseq));
xlabel('Sample size M'); ylabel('Relative Error');
title('LSE-trajectory');

subplot(122); boxplot(errLSE_seq_EnsD');
ax = gca; ax.YAxis.Scale ="log";
set(gca, 'XTickLabels',compose('%d',Mseq));
xlabel('Sample size M'); ylabel('Relative Error'); ylim([ylim_lb,ylim_ub]); 
title('LSE-ensemble');
figname = [figpath0,'fig_LSE2_traj_ens',str_name];
set_positionFontsAll;
