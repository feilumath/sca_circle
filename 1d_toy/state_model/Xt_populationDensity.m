function Xt_distr = Xt_populationDensity(Xt,stoCA)
% compute population density
N    = stoCA.N;
tN   = length(Xt(1,:)); 

edges = stoCA.edges; 

Xt_distr = zeros(stoCA.K,tN); 
for t=1:tN
    Xt_distr(:,t) = histcounts(Xt(:,t),edges);
end

Xt_distr = Xt_distr/N; 