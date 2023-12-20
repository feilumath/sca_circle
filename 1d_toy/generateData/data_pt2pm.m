function pm_all = data_pt2pm(local_p_all)
% change data from trajectories to ensembles at different times 
% use cell so that the trajectories may have different lengths; the
% ensembles may have different sample sizes
[K,N,tN] = size(local_p_all{1}); 
M        = length(local_p_all);
pm_all   = cell(1,tN);
arrayMt  = zeros(K,N,tN,M); 
for m=1:M
     arrayMt(:,:,:,m) = local_p_all{m};
end
for t=1:tN
     pm_all{t} = squeeze(arrayMt(:,:,t,:) );     
end
end
