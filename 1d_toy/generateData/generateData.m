function Xt_all = generateData(inferInfo,stoCA_par)
% generate data consisting of many trajectories

M = inferInfo.M; 

Xt_all = cell(M,1);
parfor m= 1:M
    Xt = stoCA_model(stoCA_par); 
    Xt_all{m} = Xt; 
end