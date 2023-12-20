function popuD = populationDensity_Mtraj(Xt_all,stoCA_par)
% compute population density from M trajectories in cells. 

if iscell(Xt_all) 
    M = length(Xt_all); 
    popuD = cell(1,M); 
    for m=1:M
        popuD{m} = Xt_populationDensity(Xt_all{m},stoCA_par);
    end
else
    fprintf('\n Compute population density: treat Xt_all as 1 trajectory of size NxtN \n');
    popuD{1} = Xt_populationDensity(Xt_all,stoCA_par);
end

end