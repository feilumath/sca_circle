function X = stoCA_model(stoCA,X0)
% the stochastic cellular automata model 
%  run it to generate the time series X, 
% which represents the index of elements in the alphabet set 

N    = stoCA.N;
tN   = stoCA.tN; 
X    = zeros(N,tN);

P0     = stoCA.IC_distr; 
X(:,1) = gendist(P0,N,1);
if exist('X0','var')
    X(:,1) = X0; 
end

for t= 1:tN
    temp = move1step(X(:,t),stoCA);
    X(:,t+1) = temp; 
end

end


