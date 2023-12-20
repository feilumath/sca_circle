
function test_ergodicity(stoCA_par)
% test ergodicity by running the system for a large time. 

% run a long path
stoCA_par.tN = 1e3; 
Xt = stoCA_model(stoCA_par);
tN = stoCA_par.tN; 
tInd = 1:ceil(tN/10):tN; 
Xt_distr = Xt_populationDensity(Xt(:,tInd),stoCA_par); 

% %% test if population density satisfies D(t+1) = Tmat* D(t); 
TT = stoCA_par.TMat; XX = Xt_distr(:,1:10); YY = TT*XX; 
fprintf('\n Population density D(t): D(t+1) - Tmat* D(t):\n');
disp(XX(:,2:end)-YY(:,1:end-1));

% plot path-average
N = stoCA_par.N; 
g = @(x)  x.^2; siteInd = 1:4; 
avgG_X  = zeros(length(siteInd),length(tInd)); 
for n =1:length(tInd)
    avgG_X(:,n) = mean(g(Xt(siteInd,1:tInd(n))),2);   
end
figure; plot(tInd, avgG_X(siteInd,:)','linewidth',1);  xlabel('Time'); ylabel('path average');
title('Path average'); 

% plot population density along a path
y = 1:stoCA_par.K; 
figure; 
subplot(121); imagesc(tInd,y, Xt_distr);  xlabel('Time'); ylabel('Alphabeta');   % seems not reaching a stationary distribution
title('Population density of all sites')

subplot(122);  plot(tInd, Xt_distr(1:2,:),'linewidth',1); xlabel('Time'); ylabel('Prob(k)');
legend('k=1','k=2'); 
title('Population density at each time'); 

end