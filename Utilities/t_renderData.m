close all;
clear all;
clc;

ieInit;


wave = 400:5:700; wave = wave(:);
nWaves = length(wave);

fName = fullfile(cmfRootPath,'Parameters','ximeaLights');
flashSpd = ieReadSpectra(fName,wave);
flashSpd = Energy2Quanta(wave,flashSpd);
flashNorm = flashSpd/max(flashSpd(:));


sensor = sensorCreate('bayer (rggb)');
sensor = sensorSet(sensor,'wave',wave);
sensor = sensorSet(sensor,'noise flag',0);

cameraResp = sensorGet(sensor,'filter transmissivities');
[ spd ] = xy2Spectrum( 0.25, 0.25, wave );

[measurement, reference] = renderData(wave, sensor, spd, flashNorm, 'target', 'Macbeth');