function [rawLin, raw, totalGain, mask] = getImage(path)

% function [rawLin, raw, totalGain, mask] = getImage(path)
%
% Read a raw image captured by a Ximea camera and scale by the recorded
% exposure time, gain and remove black level.
%
% It is assumed that the medatada is stored in a text file having the same
% name as the image file.
%
% Inputs:
%  path - path to the image file
%
% Returns:
%  rawLin - raw sensor data scaled by the exposure and gain
%  raw - raw sensor data
%  totalGain - the total scaling parameter applied to raw data to get
%  rawLin data
%  mask - a binary mask indicating pixels that have not been overexposed.
%
% Copyright, Henryk Blasinski 2017

rawInt = imread(path);

mask = rawInt <= 0.95*1024;

[dir, file] = fileparts(path);

fName = fullfile(dir,sprintf('%s.txt',file));
fid = fopen(fName,'r');
fgetl(fid);
res = fscanf(fid,'%f %f %f');
fclose(fid);
shutter = res(1)/1000;
gain = 10^(res(2)/20);
blackLevel = res(3);

totalGain= (shutter*gain);

raw = double(rawInt-blackLevel)/1024;

rawLin = max(double(rawInt-blackLevel)/1024,0)/totalGain;

end

