% This is a simple validation script to check if the Natural 100
% simulations work as expected and the data reconstruction is correct.
%
% Copyright, Henryk Blasinski 2017

close all;
clear all;
clc;

ieInit;

%%

wave = 400:5:700; 
wave = wave(:);
nWaves = length(wave);

% Define multispectral flash LED spectra
fName = fullfile(cmfRootPath,'Parameters','ximeaLights');
flashSpd = ieReadSpectra(fName,wave);
flashSpd = Energy2Quanta(wave,flashSpd);
flashNorm = flashSpd/max(flashSpd(:));

nChannels = size(flashSpd,2);

ledSets = {[1,4,7],...
           [1,3,4,7],...
           [1,3,4,5,7],...
           [1,2,3,4,5,7],...
           [1,2,3,4,5,6,7]};
nLedSets = length(ledSets);
       
% Define the chromaticity coordinates of the illuminant
xVec = linspace(0,1,10);
yVec = linspace(0,1,10);
nX = length(xVec);
nY = length(yVec);
       
% Define camera models
sensor = sensorCreate('bayer (rggb)');
sensor = sensorSet(sensor,'wave',wave);
sensor = sensorSet(sensor,'noise flag',0);

fName = fullfile(isetRootPath,'data','sensor','colorfilters','NikonD1.mat');
cameraResp = ieReadSpectra(fName,wave);
cameraResp(isnan(cameraResp)) = 0;
sensor = sensorSet(sensor,'filter transmissivities',cameraResp);

%% Render

x = 0.25;
y = 0.25;
[ spd ] = xy2Spectrum( x, y, wave );

measurement = renderData(wave, sensor, spd, flashNorm, 'compact', true);

%% Estimate
%  Plot how well the pixel predictions match the ambient image.

figure;
for s=1:nLedSets
    
    nLEDs = length(ledSets{s});
    currentFlash = flashNorm(:,ledSets{s});
    
    ambient = measurement.patch.ambient/max(measurement.patch.ambient(:));
    data = measurement.patch.led/max(measurement.patch.led(:));
    
    [ ambientEst, ambientWghts, ambientPredictions ] = globalAmbientEst( ambient, data, currentFlash, 'alpha', 0 );
    
    subplot(1,nLedSets,s);
    hold on; grid on; box on;
    plot(squeeze(ambient)',ambientPredictions','.');
    xlabel('Simulated');
    ylabel('Approximation');
    title(sprintf('%i LEDs',nLEDs));
    
end


