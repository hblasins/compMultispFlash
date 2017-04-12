function [ measurement, mask,  cp] = readXimeaImageStack( path, ev, nLEDs, varargin )
 
% function [ measurement, mask,  cp] = readXimeaImageStack( path, ev, nLEDs, varargin )
%
% Read image data captured with the Ximea camera
%
% Inputs:
%   path - the top level directory containing image data.
%   ev - exposure value to read (acquisitions with different exposures are
%   placed in sub-directories).
%   nLEDs - number of LEDs used during the capture process.
%
% Inputs (optional):
%   cp - corner points of the Macbeth chart present in the image.
% 
% Outputs
%   measurement - a structure containing image data in different
%   arrangements (vectorized, sampled, ambient, LED only etc.)
%   mask - a binary image indicating pixels that have been overexposed
%   cp - corner points for the macbeth chart.
%
% Copyright, Henryk Blasinski 2017


p = inputParser;
p.addOptional('cp',[]);
p.addOptional('rootPath',cmfRootPath);
p.parse(varargin{:});
inputs = p.Results;


h = 1024;
w = 1280;
wave = 400:10:700;

fName = fullfile(inputs.rootPath,'Parameters','XimeaSpectralQe');
qe = ieReadSpectra(fName,wave);

sensor = sensorCreate('bayer (bggr)');
sensor = sensorSet(sensor,'wave',wave);
sensor = sensorSet(sensor,'size',[h w]);
sensor = sensorSet(sensor,'filter transmissivities',qe);
sensor = sensorSet(sensor,'noise flag',0);
% sensor = sensorSet(sensor,'filterspectra',qe);


exposureDirs = dir(fullfile(path,sprintf('*_EV_%.1f*',ev)));


RAWlin = zeros(h,w,nLEDs+1);
RAWdemosaiced = zeros(h,w,3,nLEDs+1);

patchData = zeros(24,3,nLEDs+1);

globalMask = ones(h,w);
    
subDir = fullfile(path,exposureDirs(1).name);

[RAWlin(:,:,1), raw, totalGain, localMask] = getImage(fullfile(subDir,'Ambient.tiff'));

globalMask = globalMask & localMask;

sensor = sensorSet(sensor,'volts',raw);
sensor = sensorSet(sensor,'name',sprintf('Ambient'));
ieAddObject(sensor);
sensorWindow();

[data, ~, ~, inputs.cp] = macbethSelect(sensor,0,1,inputs.cp);
patchData(:,:,1) = cell2mat(cellfun(@nanmean,data,'UniformOutput',false)')/totalGain;

ip = ipCreate();
ip = ipSet(ip,'sensor conversion method','none');
ip = ipCompute(ip,sensor);
ip = ipSet(ip,'name',sensorGet(sensor,'name'));
ieAddObject(ip);

RAWdemosaiced(:,:,:,1) = ipGet(ip,'sensor channels')/totalGain;


for l=1:nLEDs
    [RAWlin(:,:,l+1), raw, totalGain, localMask] = getImage(fullfile(subDir,sprintf('LED%02i.tiff',l)));
    
    globalMask = globalMask & localMask;
    
    sensor = sensorSet(sensor,'volts',raw);
    sensor = sensorSet(sensor,'name',sprintf('Ambient+LED%i',l));
    
    ieAddObject(sensor);
    sensorWindow();
    
    [data, ~, ~, inputs.cp] = macbethSelect(sensor,1,1,inputs.cp);
    patchData(:,:,l+1) = cell2mat(cellfun(@nanmean,data,'UniformOutput',false)')/totalGain;

    ip = ipCreate();
    ip = ipSet(ip,'sensor conversion method','none');
    ip = ipCompute(ip,sensor);
    ip = ipSet(ip,'name',sensorGet(sensor,'name'));
    ieAddObject(ip);
    
    RAWdemosaiced(:,:,:,l+1) = ipGet(ip,'sensor channels')/totalGain;
    
end


sensorWindow();
ipWindow();
  
RAWlin = RAWlin/max(RAWlin(:));
RAWdemosaiced = RAWdemosaiced/max(RAWdemosaiced(:));

measurement.raw.ambient = RAWlin(:,:,1);
measurement.raw.led = max(RAWlin(:,:,2:end) - repmat(measurement.raw.ambient,[1 1 nLEDs]),0);

measurement.demosaiced.ambient = RAWdemosaiced(:,:,:,1);
measurement.demosaiced.led = max(RAWdemosaiced(:,:,:,2:end) - repmat(measurement.demosaiced.ambient,[1 1 1 nLEDs]),0);


measurement.vectorized.ambient = reshape(measurement.demosaiced.ambient,[h*w,3])';
measurement.vectorized.led = permute(reshape(measurement.demosaiced.led,[h*w, 3, nLEDs]),[2, 3, 1]);


measurement.downsampled.ambient = imresize(measurement.demosaiced.ambient,0.25,'nearest');
measurement.downsampled.led = imresize(measurement.demosaiced.led,0.25,'nearest');

maskDownsampled = logical(imresize(globalMask,0.25,'nearest'));

hd = size(measurement.downsampled.ambient,1);
wd = size(measurement.downsampled.ambient,2);

measurement.downsampled.ambient = reshape(measurement.downsampled.ambient,[hd*wd, 3])';
measurement.downsampled.led = permute(reshape(measurement.downsampled.led,[hd*wd, 3, nLEDs]),[2, 3, 1]);

measurement.downsampled.ambient = measurement.downsampled.ambient(:,maskDownsampled(:));
measurement.downsampled.led = measurement.downsampled.led(:,:,maskDownsampled(:));

measurement.patch.ambient = patchData(:,:,1);
measurement.patch.led = permute(max(patchData(:,:,2:end) - repmat(measurement.patch.ambient,[1 1 nLEDs]),0),[2, 3, 1]);
measurement.patch.ambient = measurement.patch.ambient';

nF = max([measurement.patch.ambient(:); measurement.patch.led(:)]);

measurement.patch.ambient = measurement.patch.ambient/nF;
measurement.patch.led = measurement.patch.led/nF;

mask = globalMask;

cp = inputs.cp;

end

%{
function [rawLin, raw, totalGain, mask] = getImage(path)

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
%}



