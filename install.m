% Add relevant directories with source code to Matlab path.
%
% Copyright, Henryk Blasinski 2017

close all;
clear all;
clc;

addpath(cmfRootPath);

addpath(fullfile(cmfRootPath,'Utilities'));
addpath(fullfile(cmfRootPath,'Utilities','Natural100'));

addpath(fullfile(cmfRootPath,'Estimation'));
addpath(fullfile(cmfRootPath,'Estimation','Global'));
addpath(fullfile(cmfRootPath,'Estimation','Spatial'));

% If present add a path to the directory containing hardware control
% scripts.

hwPath = fullfile(cmfRootPath,'Hardware');
if exist(hwPath,'dir')
    addpath(genpath(hwPath));
end