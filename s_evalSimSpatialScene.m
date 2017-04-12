% This script demonstrates how the complementary spectral flash works
% on simulated data. The data was obtained using RenderToolbox4 and PBRT.
%
% THIS SCRIPT EVALUATES THE SPATIAL METHOD THAT IS NOT DESCRIBED IN THE
% PAPER
%
% Copyright, Henryk Blasinski 2017


close all;
clear all;
clc;

ieInit;
set(groot,'defaultAxesColorOrder',[1 0 0; 0 1 0; 0 0 1]);

wave = 400:10:700;
nWaves = length(wave);


fName = fullfile(cmfRootPath,'Parameters','XimeaSpectralResponsivities');
cameraMat = ieReadColorFilter(wave,fName);

sensor = sensorCreate('bayer (rggb)');
sensor = sensorSet(sensor,'wave',wave);
sensor = sensorSet(sensor,'size',[480 640]);
sensor = sensorSet(sensor,'noise flag',2);
sensor = sensorSet(sensor,'filter transmissivities',cameraMat);


chlorophyllConc = 0.0;
cdomConc = 0.0;
camDist = 1000;
waterDepth = 10000;
smallPart = 0.0;
largePart = 0.0;
ambient = 'D65';

target = 'Macbeth'; %or 'Macbeth' or 'Table' or 'Acropora'


fName = fullfile(slRootPath,'Parameters','ximeaLights');
ledSpectra = ieReadSpectra(fName,wave);

flashNorm = Energy2Quanta(wave,ledSpectra);
flashNorm = flashNorm/max(flashNorm(:));
ledSpectra = [zeros(nWaves,1), 1e7*flashNorm];

ill = illuminantCreate('D65',wave);
surfaceSpectrum = illuminantGet(ill,'photons');
surfaceSpectrum = surfaceSpectrum/max(surfaceSpectrum(:));


nChannels = size(flashNorm,2);
[measurement, reference] = getRenderedData(sensor, wave, surfaceSpectrum, flashNorm,...
                                    'targetDistance',camDist,...
                                    'depth',waterDepth,...
                                    'chlConc',chlorophyllConc,...
                                    'cdomConc',cdomConc,...
                                    'smallPartConc',smallPart,...
                                    'largePartConc',largePart,...
                                    'target',target,...
                                    'rtbResultFolder',fullfile(cmfRootPath,'..','Data','Simulated'));


%% Ambient estimate
%  Estimate the ambient illuminant independently for every pixel, with some
%  spatial smoothing

alpha = 0;
beta = 0;

[ambientApproxWghts, pred, history] = ambientEstADMM( measurement.demosaiced.ambient, ...
    measurement.demosaiced.led, flashNorm, alpha, beta,...
    'maxIter', 100, 'verbose', true , 'rescaleRho',true, 'tol',0);

figure;
hold on; grid on; box on;
plot(measurement.demosaiced.ambient(:),pred(:),'.');
xlabel('Captured');
ylabel('Modeled');


% Plot the (spatially varying) estimated ambient spectrum
ambientEst = flashNorm*reshape(ambientApproxWghts,[640*480, nChannels])';
pointIDs = randi(size(ambientEst,2),1000,1);

figure;
hold on; grid on; box on;
plot(wave,ambientEst(:,pointIDs),'k');
plot(wave,cameraMat,'--');
plot(wave,flashNorm,':');
xlabel('Wavelength, nm');

figure;
plot([history.prRes history.dualRes]);
xlabel('Iteration')
legend('Primal','Dual');
title('ADMM residuals');


figure;
for i=1:nChannels
    subplot(2,4,i);
    imagesc(ambientApproxWghts(:,:,i));
    axis image;
end

%% Depth estimation
% Given the weights on LEDs to approximate the ambient, we should be able
% to derive the depth (larger weights -> more LED light needed to
% illuminate given pixel -> the pixel is further away)
% We compute the 'intensity' as the weight along the principal component of
% the illumination spectrum
wghtsVec = reshape(ambientApproxWghts,[640*480, nChannels])';

[vec, ~, score ] = pca(wghtsVec');

% The 1st principal component describes the global illuminant
meanAmbientEst = flashNorm*vec(:,1);

figure;
hold on; grid on; box on;
plot(wave,meanAmbientEst);
xlabel('Wavelength, nm');
title('Global illuminant');

wref = 1;

comp1Wght = (vec(:,1)'*wghtsVec)/wref;
scaledDistance = 1./abs(sqrt(reshape(comp1Wght,[480, 640])));

figure;
imagesc(scaledDistance,[0 10]); colorbar;
title('Estimated scaled distance');


%% Estimate attenuation on the path between the camera and the target

% We assume we know the distance and the depth.
% distance = 1000;
% depth = 10000;

% exponent = repmat((2*distance/depth*scaledDistance),[1 1 nWaves]);
% attnEst = repmat(shiftdim(meanAmbientEst,-2),[hh ww 1]).^exponent;


%% Complement estimate

alpha = 0;
beta = 0;

desiredIll = 'D65';
ill = illuminantCreate(desiredIll,wave);
desiredSpectrum = illuminantGet(ill,'photons');
desiredSpectrum = desiredSpectrum/max(desiredSpectrum);

[ complWghts, slack, history ] = complementEstADMM( desiredSpectrum, ambientApproxWghts, flashNorm, cameraMat, alpha, beta,...
                                                      'maxIter', 200, 'verbose', true, 'rescaleRho', false, 'flashMode',false );
%{
[ complWghtsAttn, slack, hist ] = complementEstWithAttnADMM( desiredSpectrum, ambientApproxWghts, flashNorm, attnEst, cameraMat, alpha, beta,...
                                                      'maxIter', 200, 'verbose', true, 'rescaleRho', false, 'flashMode',false );
%}
complementEst = flashNorm*reshape(complWghts,[480*640, nChannels])';
figure;
plot(wave,complementEst(:,pointIDs),'k');
xlabel('Wavelength, nm');
ylabel('Intensity');
title('Complementary flash');

figure;
plot(wave,complementEst(:,pointIDs) + ambientEst(:,pointIDs),'r');
xlabel('Wavelength, nm');
ylabel('Intensity');
title('Total illumination');
                                                  

%% Computational synthesis
% Render an image as if captured with the complementary flash.

% No camera-target distance correction
rendering = measurement.demosaiced.ambient;

for i=1:nChannels
    rendering = rendering + repmat(complWghts(:,:,i),[1 1 3]).*measurement.demosaiced.led(:,:,:,i);
end

rendering = rendering/max(rendering(:));

sensor = sensorSet(sensor,'volts',rendering);
sensor = sensorSet(sensor,'name',sprintf('Computational: %s',desiredIll));
ieAddObject(sensor);
sensorWindow;


ip = ipCreate;
ip = ipCompute(ip,sensor);
ip = ipSet(ip,'name','Computational');
ieAddObject(ip);

ipWindow;

%% Camera-target distance correction
%{
rendering = images.ambientDemosaiced;

for i=1:nChannels
    rendering = rendering + repmat(complWghtsAttn(:,:,i),[1 1 nFilters]).*images.ledsDemosaiced(:,:,:,i);
end

rendering = rendering/max(rendering(:));

sensor = sensorSet(sensor,'volts',rendering);
sensor = sensorSet(sensor,'name',sprintf('Computational+Correction: %s',desiredIll));
ieAddObject(sensor);
sensorWindow;


ip = ipCreate;
ip = ipCompute(ip,sensor);
ip = ipSet(ip,'name','Computational+Correction');
ieAddObject(ip);

ipWindow;

%}





