function pdf_sites= siteDensity(Xm_all,K)
% compute the marginal densities of each site from ensemble of samples
% Need sample size reasonably large.
tN        = length(Xm_all);
[N,~]     = size(Xm_all{1});
pdf_sites = zeros(K,N,tN); 
for t=1:tN
     temp  = Xm_all{t}; 
     Mt = length(temp(1,:)); 
     % Mhalf = ceil(M/2);
     % random_M = randi([1,M],1,Mhalf);
     for n=1:N
          % pdf_n = histcounts(temp(n,random_M)); % use only half of M samples.   
          pdf_n = histcounts(temp(n,:),K); % use all samples
          pdf_sites(:,n,t) = pdf_n'; 
     end
     pdf_sites(:,:,t) = pdf_sites(:,:,t)/Mt;  %  
end


end