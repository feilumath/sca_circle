function Xm_all = data_Xt2Xm(Xt_all)
% change data from trajectories to ensembles at different times 
% use cell so that the trajectories may have different lengths; the
% ensembles may have different sample sizes
[N,tN] = size(Xt_all{1}); 
M        = length(Xt_all);
Xm_all   = cell(1,tN);
arrayMt  = zeros(N,tN,M); 
for m=1:M
     arrayMt(:,:,m) = Xt_all{m};
end
for t=1:tN
     Xm_all{t} = squeeze(arrayMt(:,t,:) );     
end
end
