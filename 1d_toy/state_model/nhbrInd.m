
function Ind_nhbr_k = nhbrInd(k,N,nhbrSize)
% the index of sites in the neighbor of site k: cyclic indexing
siteInd_ext = [N-nhbrSize+1:N,1:N,N+1:N+nhbrSize-1];

Ind_nhbr_k = siteInd_ext( (-nhbrSize:nhbrSize)+k + nhbrSize); 
    
end
