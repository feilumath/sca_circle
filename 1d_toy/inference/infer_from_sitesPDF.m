
function [Tmat_lse,lse_struct] = infer_from_sitesPDF(Xm_all,local_p_all_M,K)
% % estimate Tmat when data are ensemble without trajectory information
% Idea: match the marginal density of each site. Need large ensemble size
% fprintf('\n Estimator from ensemble data without trajectory information\n');
pdf_sites     = siteDensity(Xm_all,K);               % size = KxNxtN
mean_phi      = mean_localDensity(local_p_all_M);    % size = KxNxtN    
[K,N,tN]      = size(pdf_sites);

%% LSE estimator: min || pdf_sites - Tmat*mean_phi ||^2  >> pdf_sites*phi' =   Tmat*phi*phi'
%           >> (phi*phi')'*Tmat' = phi*pdf_sites'
Amat = zeros(K,K);  % mean_phi(t)*mean_phi(t)'         with phi KxN for each t 
bvec = zeros(K,K);  % mean_phi(t)*pdf_sites(t)'          >>> lse is column stochastic: sum_j T(k,j) =1
lse_all = cell(1,tN);
for t = 1:tN-1
    Amat1t = zeros(K,K); bvec1t = zeros(K,K);
    phi    = mean_phi(:,:,t);    % size = KxN
    pdf1t  = pdf_sites(:,:,t+1); % size KxN, to get transpose to have value in [1,K];
    Amat1t = Amat1t+ phi*phi';
    bvec1t = bvec1t+ phi*pdf1t';
    lse_all{t}.Amat = Amat1t; 
    lse_all{t}.bvec = bvec1t;
    Amat = Amat +  Amat1t; 
    bvec = bvec + bvec1t;
end

lse_nonneg =zeros(K,K); 
for k=1:K
     lse_nonneg(:,k)  = lsqnonneg(Amat,bvec(:,k)); 
end

lse_struct.cvec = lse_nonneg; %  Amat\bvec;
lse_struct.Amat = Amat; 
lse_struct.bvec = bvec; 
Tmat_lse = lse_nonneg;  
 Tmat_lse = diag(1./sum(Tmat_lse,2))*Tmat_lse;
end
