% This script demonstrates how the complementary spectral flash works
% on simulated data. The data was obtained using RenderToolbox4 and PBRT.
%
% Copyright, Henryk Blasinski 2017

close all;
clear all;
clc;

ieInit;

% Directory where images and figures will be saved. If empty noting is
% saved.
% destDir = fullfile(cmfRootPath,'..','Figures');
destDir = [];


set(groot,'defaultAxesColorOrder',[1 0 0; 0 1 0; 0 0 1]);

wave = 400:5:700;
nWaves = length(wave);

fName = fullfile(cmfRootPath,'Parameters','XimeaSpectralResponsivities');
cameraResp = ieReadColorFilter(wave,fName);

sensor = sensorCreate('bayer (bggr)');
sensor = sensorSet(sensor,'wave',wave);
sensor = sensorSet(sensor,'size',[480 640]);
sensor = sensorSet(sensor,'noise flag',2);
sensor = sensorSet(sensor,'filter transmissivities',cameraResp);

% Define the water properties (Note, data may not have been rendered for
% the properties you've specified).

chlorophyllConc = 0.0;
cdomConc = 0;
camDist = 1000;
waterDepth = 10000;
smallPart = 0.0;
largePart = 0.0;

target = 'Objects'; %or 'Macbeth' or 'Table' or 'Acropora' or 'Objects'

fName = fullfile(cmfRootPath,'Parameters','ximeaLights');
ledSpectra = ieReadSpectra(fName,wave);

flashNorm = Energy2Quanta(wave,ledSpectra);
flashNorm = flashNorm/max(flashNorm(:));

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
%
% Use all pixels in the image to predict the ambient spectrum


[ ambientEst, ambientWghts ] = globalAmbientEst( measurement.vectorized.ambient,...
    measurement.vectorized.led, ...
    flashNorm,...
    'alpha',1);


% Plot the estimated ambient spectrum
figure;
hold on; grid on; box on;
plot(wave,ambientEst);
legend('Estimated');
title('Whole image ambient estimate');


% Render the approximation image

rendering = zeros(size(measurement.raw.ambient));
for i=1:nChannels
    rendering = rendering + ambientWghts(i)*measurement.raw.led(:,:,i);
end

rendering = rendering/max(rendering(:));

sensor = sensorSet(sensor,'volts',rendering);
sensor = sensorSet(sensor,'name','Ambient approximation');
ieAddObject(sensor);
sensorWindow;


ip = ipCreate;
ip = ipCompute(ip,sensor);
ip = ipSet(ip,'name','Ambient approximation');
ieAddObject(ip);

ipWindow;


%% Per-patch ambient estimate


[ ambientEstPatch, ambientWghtsPatch, ambientPredictions ] = globalAmbientEst( measurement.patch.ambient,...
    measurement.patch.led, ...
    flashNorm,...
    'alpha',1);


% Plot the quality of the approximation in the camera RGB space

figure;
hold on; grid on; box on;
plot(squeeze(measurement.patch.ambient)',ambientPredictions','o');
xlabel('Ambient appearance');
ylabel('Approximated');
title('RGB space - individual patches');

% Plot the estimated ambient spectrum
figure;
hold on; grid on; box on;
plot(wave,ambientEstPatch);
legend('Estimated');
title('Patch based ambient estimate');


%% Estimate attenuation on the path between the camera and the target

% We assume we know the distance and the depth.
distance = camDist;
depth = waterDepth;

surfaceGuess = 'D65';
ill = illuminantCreate(surfaceGuess,wave);
spectrum = illuminantGet(ill,'photons');
spectrum = spectrum/max(spectrum(:));

attnEst = ((ambientEst/max(ambientEst(:)))./spectrum).^(2*distance/depth);


%% Complement estimate

desiredIll = 'D65';
ill = illuminantCreate(desiredIll,wave);
desiredSpectrum = illuminantGet(ill,'photons');

[ flashEst, flashWghts ] = globalComplementEst( desiredSpectrum, ambientEst, diag(attnEst)*flashNorm, cameraResp,...
                            'flashMode',true);


% Plot the illuminant spectra
figure;
hold on; grid on; box on;
plot(wave,[ambientEst flashNorm*flashWghts  ambientEst+flashEst])
xlabel('Wavelength');
ylabel('Scaled photons');
legend('Ambient','Comlp flash','Total');


%% Computational synthesis
% Render an image as if captured with the complementary flash.

rendering = measurement.raw.ambient;
for i=1:nChannels
    rendering = rendering + flashWghts(i)*measurement.raw.led(:,:,i);
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


%% Small camera-taregt distance
%  We can assume that the attenuation along this path is

attnEst = ones(nWaves,1);


% Complement estimate

desiredIll = 'D65';
ill = illuminantCreate(desiredIll,wave);
desiredSpectrum = illuminantGet(ill,'photons');

[ flashEst, flashWghts ] = globalComplementEst( desiredSpectrum, ambientEst, diag(attnEst)*flashNorm, cameraResp,...
                                                'flashMode',true);


% Plot the illuminant spectra
figure;
hold on; grid on; box on;
plot(wave,[ambientEst flashNorm*flashWghts  ambientEst+flashEst])
xlabel('Wavelength');
ylabel('Scaled photons');
legend('Ambient','Comlp flash','Total');

%% Plots for print

fs=10;

fg = figure;
hold on; grid off; box off;
plot(wave,ambientEst,'g','LineWidth',2);
ylim([0, max(flashNorm*flashWghts)]);
xlabel('Wavelength, nm','Interpreter','Latex');
ylabel('Scaled photons','Interpreter','Latex');
set(gca,'FontSize',fs+2);
set(gca,'TickLabelInterpreter','Latex');
set(gcf,'Units','Centimeters');
set(gcf,'PaperPosition',[1 1 4.5 2]);
set(gca,'XTick',400:100:800);

if isempty(destDir) == false
    print('-depsc',fullfile(destDir,'Ambient.eps'));
end

figure;
hold on; grid off; box off;
plot(wave,ambientEst,'g','LineWidth',2);
plot(wave,flashNorm*flashWghts,'r--','LineWidth',2);
ylim([0, max(flashNorm*flashWghts)]);
xlabel('Wavelength, nm','Interpreter','Latex');
set(gca,'TickLabelInterpreter','Latex');
ylabel('Scaled photons','Interpreter','Latex');
set(gcf,'Units','Centimeters');
set(gca,'FontSize',fs+2);
set(gcf,'PaperPosition',[1 1 4.5 2]);
set(gca,'XTick',400:100:800);

if isempty(destDir) == false
    print('-depsc',fullfile(destDir,'Ambient+Flash.eps'));
end


%% Computational synthesis
% Render an image as if captured with the complementary flash.

rendering = measurement.raw.ambient;
for i=1:nChannels
    rendering = rendering + flashWghts(i)*measurement.raw.led(:,:,i);
end

rendering = rendering/max(rendering(:));

sensor = sensorSet(sensor,'volts',rendering);
sensor = sensorSet(sensor,'name',sprintf('Computational (noCorrection): %s',desiredIll));
ieAddObject(sensor);
sensorWindow;


ip = ipCreate;
ip = ipCompute(ip,sensor);
ip = ipSet(ip,'name','Computational (no distance corr.)');
ieAddObject(ip);

ipWindow;

%% Save data
% destDir = fullfile(cmfRootPath,'..','Figures');
ieImages = vcGetObjects('vcimage');

for i=1:length(ieImages)
   
    image = ipGet(ieImages{i},'data srgb');
    name = ipGet(ieImages{i},'name');
    figure; 
    imshow(image);
    title(name);
    if isempty(destDir) == false
        imwrite(image,fullfile(destDir,sprintf('%s_%s.png',target,name)));
    end
end
    




