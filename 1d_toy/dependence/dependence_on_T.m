
function [diff_Tmat1,diff_Pmat1,diff_pi1,diff_l1]=dependence_on_T(stoCA,K,N,nsimu,plotON)
% test dependence of P and inv_pi on T

% nsimu = 100; 
Tmat   = stoCA.TMat; 
Pmat   = trans_prob_Mat_markov(stoCA,K,N); 
[V,~]  = eig(Pmat');
inv_pi = V(:,1)'/sum(V(:,1));


diff_Tmat = zeros(1,nsimu);
diff_Pmat = zeros(1,nsimu);
diff_pi   = zeros(1,nsimu);


diff_Tmat_l1  = zeros(1,nsimu);
diff_Pmat_l1  = zeros(1,nsimu);
diff_pi_l1    =  zeros(1,nsimu);

parfor m=1:nsimu
    T_perturb   = rand(size(Tmat));                                        % to make it stochastic 
    Tmat1 = Tmat+T_perturb;
    Tmat1 = Tmat1*diag(1./sum(Tmat1));
    stoCA1 = stoCA; 
    stoCA1.TMat = Tmat1;
    Pmat1       = trans_prob_Mat_markov(stoCA1,K,N);
    [V1,~]  = eig(Pmat1'); 
    inv_pi1 = V1(:,1)'/sum(V1(:,1));   

    diff_Pmat(m)= norm(Pmat1-Pmat,"fro");
    diff_Tmat(m)= norm(T_perturb,"fro");
    diff_pi(m) = norm(inv_pi1-inv_pi);

    diff_Pmat_l1(m) = sum(abs(Pmat1-Pmat),"all");
    diff_Tmat_l1(m) = sum(abs(T_perturb),'all');
    diff_pi_l1(m)   = norm(inv_pi1-inv_pi,1); 
end

[diff_Tmat1,ind]    = sort(diff_Tmat);
diff_Pmat1          = diff_Pmat(ind);
diff_pi1            = diff_pi(ind);

[diff_Tmat1_l1,ind]    = sort(diff_Tmat_l1);
diff_Pmat1_l1          = diff_Pmat_l1(ind);
diff_pi_l1             = diff_pi_l1(ind);

diff_l1.diff_Tmat_l1 = diff_Tmat1_l1; 
diff_l1.diff_Pmat_l1 = diff_Pmat1_l1; 
diff_l1.diff_pi_l1    = diff_pi_l1; 

if plotON ==1
    figure;
    scatter(diff_Tmat1,diff_Pmat1); hold on;
    scatter(diff_Tmat1,diff_pi1);
    legend('P change','pi change')

    figure;
    semilogy(diff_Tmat1,diff_Pmat1./diff_Tmat1,'linewidth',1); hold on;
    semilogy(diff_Tmat1,diff_pi1./diff_Tmat1,'linewidth',2);
    legend('|P-diff|/|T-dff|','|pi-diff|/|T-dff|')
    xlabel('Frobenus norm T-diff');

   %  figname = 'dependence_on_T';
   %  set_positionFontsAll;
end
end
