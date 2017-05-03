% This script demonstrates how the complementary spectral flash works
% on real data captured in a swimming pool at about 5m.
%
% Copyright, Henryk Blasinski 2017

close all;
clear all;
clc;

ieInit;

%%
set(groot,'defaultAxesColorOrder',[1 0 0; 0 1 0; 0 0 1]);

wave = 400:5:700;
nWaves = length(wave);

desiredIll = 'D65';

% If the resDir is empty, no fiures will be saved.
% resDir = fullfile(cmfRootPath,'..','Figures');
resDir = [];


fName = fullfile(cmfRootPath,'Parameters','ximeaLights');
flashSpd = ieReadSpectra(fName,wave);
flashSpd = Energy2Quanta(wave,flashSpd);
flashNorm = flashSpd/max(flashSpd(:));
nLEDs = size(flashSpd,2);


fName = fullfile(cmfRootPath,'Parameters','macbethChart');
reflRef = ieReadSpectra(fName,wave);

% Create base sensor model
fName = fullfile(cmfRootPath,'Parameters','XimeaSpectralResponsivities');
cameraResp = ieReadColorFilter(wave,fName);

sensorBase = sensorCreate('bayer (bggr)');
sensorBase = sensorSet(sensorBase,'wave',wave);
sensorBase = sensorSet(sensorBase,'size',[1024 1280]);
sensorBase = sensorSet(sensorBase,'noise flag',0);
sensorBase = sensorSet(sensorBase,'filter transmissivities',cameraResp);

cameraMat = sensorGet(sensorBase,'filter transmissivities');
nFilters = size(cameraMat,2);



dataDir = fullfile(cmfRootPath,'..','Data','Pool','2016-02-05_22-23-22');
cp = [459 644;752 653;754 467;466 458];
[ measurement, mask, cp ] = readXimeaImageStack(dataDir,2,nLEDs,'cp',cp);



%% Ambient estimate

% Global: all image pixels
[ ambientEst, ambientWghts, ambientPredictions ] = globalAmbientEst( measurement.downsampled.ambient,...
    measurement.downsampled.led, ...
    flashNorm,...
    'alpha',1e2);

% Patch: Macbeth patches only
[ ambientEstPatch, ambientWghtsPatch, ambientPredictionsPatch ] = globalAmbientEst( measurement.patch.ambient,...
    measurement.patch.led, ...
    flashNorm,...
    'alpha',1);

% Plot the quality of the approximation in camera RGB space

figure;
hold on; grid on; box on;
plot(squeeze(measurement.patch.ambient)',ambientPredictionsPatch','o');
xlabel('Ambient appearance');
ylabel('Approximated');
title('RGB space');
    
    
%% Render ambient approximation image

[~, ip] = renderFlashImage(zeros(size(measurement.raw.ambient)),measurement.raw.led,ambientWghts,...
    sensorBase,'name','Ambient approximation');


image = ipGet(ip,'data srgb');
figure; imshow(image,'Border','tight');

if isempty(resDir) == false
    fName = fullfile(resDir,sprintf('Pool_ambientApprox.eps'));
    print('-depsc',fName);
end


%% Complement estimate
 
ill = illuminantCreate(desiredIll,wave);
desiredSpectrum = illuminantGet(ill,'photons');

[ flashEst, flashWghts ] = globalComplementEst( desiredSpectrum, ambientEst, flashNorm, cameraResp,...
    'flashMode',true);

% Plot the illuminant spectra
figure;
hold on; grid on; box on;
plot(wave,ambientEst,'LineWidth',2);
plot(wave,flashNorm*flashWghts,'LineWidth',2);
xlabel('Wavelength, nm','Interpreter','Latex','FontSize',6);
set(gca,'TickLabelInterpreter','Latex');
mx = max([ambientEst(:); flashNorm*flashWghts]);
ylim([-0.05*mx 1.05*mx]);
set(gcf,'Units','Centimeters');
set(gca,'FontSize',6);
set(gcf,'PaperPosition',[1 1 4 3.25]);

if isempty(resDir) == false
    fName = fullfile(resDir,sprintf('Pool_spectra.eps'));
    print('-depsc',fName);
end

 
%% Complementary flash
% Render an image as if captured with the complementary flash.

[~, ip] = renderFlashImage(measurement.raw.ambient,measurement.raw.led,flashWghts,...
    sensorBase,'name',sprintf('Computational-complement : %s',desiredIll));


image = ipGet(ip,'data srgb');
figure; imshow(image,'Border','tight');

if isempty(resDir) == false
    fName = fullfile(resDir,sprintf('Pool_flashComplement.eps'));
    print('-depsc',fName);
end


%% Matching flash
 
matchingAmbWghts = ambientWghts./max(ambientWghts);
 
[~, ip] = renderFlashImage(measurement.raw.ambient,measurement.raw.led,matchingAmbWghts,...
     sensorBase,'name',sprintf('Computational-matching : %s',desiredIll));
 
 
 
image = ipGet(ip,'data srgb');
figure; imshow(image,'Border','tight');
 
if isempty(resDir) == false
    fName = fullfile(resDir,sprintf('Pool_flashMatch.eps'));
    print('-depsc',fName);
end



%% Save data
% destDir = fullfile(cmfRootPath,'..','Figures');
destDir = [];
ieImages = vcGetObjects('vcimage');

for i=1:length(ieImages)
   
    image = ipGet(ieImages{i},'data srgb');
    name = ipGet(ieImages{i},'name');
    figure; 
    imshow(image);
    title(name);
    if isempty(destDir) == false
        imwrite(image,fullfile(destDir,sprintf('Pool_%s.png',name)));
    end
end