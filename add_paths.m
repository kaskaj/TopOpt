
% This script adds the absolute location of
% the shared functions to the path

path1 = cd;
cd ..           
path2 = cd;
if ( isunix )
    addpath(genpath([path2,'/path/']));
else
    addpath(genpath([path2,'\path\']));
end
cd(path1)
clear all

addpath('../PlotData')
addpath('../PrepareData')
addpath('../3rd_party')
addpath('../Tests')