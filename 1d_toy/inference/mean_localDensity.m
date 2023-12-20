function mean_phi = mean_localDensity(local_p_all_M)
% compute the empirical of E[phi]: Tmat*E[phi] provides site densities
% Need sample size reasonably large.
tN          = length(local_p_all_M);
[K,N,M]     = size(local_p_all_M{1});
mean_phi    = zeros(K,N,tN); 
for t=1:tN
     temp            = local_p_all_M{t};
     Mhalf    = ceil(M/2);
     random_M = randi([1,M],1,Mhalf);
     mean_phi(:,:,t) = mean(temp(:,:,random_M),3); % E[phi]  mean of phi
end
end