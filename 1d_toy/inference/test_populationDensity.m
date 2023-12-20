

%% test the population density

    popuDistr = populationDensity_Mtraj(Xt_all,stoCA_par); % the true population density

    Tmat0     = inferInfo.T_matIC; Tmat0 = stoCA_par.TMat; 
    % 1. most likely states from local_p
    [X1_mle_all,~] = maximal_local_Tphi(local_p_all,Tmat0); 
    popuDistr1     = populationDensity_Mtraj(X1_mle_all,stoCA_par);

    aa = popuDistr{1}; bb= popuDistr1{1};
    figure; 
    subplot(211);  plot(aa'); title('True population density');
    subplot(212);  plot(bb'); title('True Tmat >>> MLE sates population density');
    
  % test if population density satisfies D(t+1) = Tmat* D(t); 
   aa2 = Tmat0*aa; 
fprintf('\n Population density D(t): D(t+1) - Tmat* D(t):\n');
disp(aa2(:,1:end-1)-aa(:,2:end));