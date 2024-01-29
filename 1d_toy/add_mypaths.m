%% 
folder = fileparts(which(mfilename));
addpath(genpath(folder));


parentpath = pwd;

addpath([parentpath '/state_model/']);
addpath([parentpath '/settings/']); 
addpath([parentpath '/generateData/']);
addpath([parentpath '/inference/']);
addpath([parentpath '/plotsFn/']);
% addpath([parentpath '/output/']);
addpath([parentpath '/dependence/']);

%% save data to your local DIR:  
% making DIR to your root DIR
if ispc
    SAVE_DIR = [getenv('USERPROFILE'), '\DataAnalysis\StoCA\'];
else
    SAVE_DIR = [getenv('HOME'),'/DataAnalysis/StoCA/'];     
end
if ~exist(SAVE_DIR,'dir'), mkdir(SAVE_DIR); end    
addpath(SAVE_DIR);


