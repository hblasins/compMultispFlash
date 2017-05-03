% This is a highly experimental script using local illuminant estimation
% and global spatial smoothness
%
% USE AT YOUR OWN RISK
%
% Copyright, Henryk Blasinski 2017

close all;
clear all;
clc;

ieInit;
set(groot,'defaultAxesColorOrder',[1 0 0; 0 1 0; 0 0 1]);

wave = 400:5:700;
nWaves = length(wave);

fName = fullfile(cmfRootPath,'Parameters','ximeaLights');
flashSpd = ieReadSpectra(fName,wave);
flashSpd = Energy2Quanta(wave,flashSpd);
flashNorm = flashSpd/max(flashSpd(:));

nChannels = size(flashSpd,2);

fName = fullfile(cmfRootPath,'Parameters','macbethChart');
reflRef = ieReadSpectra(fName,wave);

fName = fullfile(cmfRootPath,'Parameters','XimeaSpectralResponsivities.mat');
camera = ieReadSpectra(fName,wave);

[raw, rawScaled] = readImageSet(fullfile(cmfRootPath,'..','Data','Pool','2016-02-05_22-23-22'),1);
rawScaled = rawScaled/max(rawScaled(:));

hh = size(raw,1);
ww = size(raw,2);

% Create a sensor model
sensor = sensorCreate('bayer (bggr)');
sensor = sensorSet(sensor,'wave',wave);
sensor = sensorSet(sensor,'size',[hh ww]);
sensor = sensorSet(sensor,'noise flag',0);

cameraMat = sensorGet(sensor,'filter transmissivities');
nFilters = size(cameraMat,2);

measVals = zeros(3,nChannels+1,24);
measDiffVals = zeros(3,nChannels,24);

rawLin = zeros(hh,ww,nFilters,nChannels+1);

%% Camera simulation

cp = [454 642;752 656;749 463;466 458];
for i=1:(nChannels + 1)
    
    sensor = sensorSet(sensor,'volts',rawScaled(:,:,i));
   
    switch i
        case 1
            name = 'Ambient';
            sensorAmbient = sensor;
            
        otherwise
            name = sprintf('Ambient + LED %i',i-1);
            
            sensorLed = sensorSet(sensor,'volts',max(sensorGet(sensor,'volts') - sensorGet(sensorAmbient,'volts'),0));
            sensorLed = sensorSet(sensorLed,'name',sprintf('LED %i',i-1)); 
            ieAddObject(sensorLed);

    end
    sensor = sensorSet(sensor,'name',name);
    ieAddObject(sensor);
    
    [mVals, ~, ~, cp] = macbethSelect(sensor,1,1,cp);
    measVals(:,i,:) = cell2mat(cellfun(@(x) nanmean(x)',mVals,'UniformOutput',false));
    
    if i>=2
        
        vcAddObject(sensorLed);
    
        [mVals, ~, ~, cp] = macbethSelect(sensorLed,1,1,cp);
        measDiffVals(:,i-1,:) = cell2mat(cellfun(@(x) nanmean(x)',mVals,'UniformOutput',false));
            
    end
    sensorWindow();
    
    
    ip = ipCreate();
    ip = ipCompute(ip,sensor);
    rawLin(:,:,:,i) = ipGet(ip,'sensor channels');
    ip = ipSet(ip,'name',name);
    ieAddObject(ip);

    
    if i>=2
        ipLed = ipCompute(ip,sensorLed);
        ipLed = ipSet(ipLed,'name',sprintf('LED %i',i-1)); 
    
        ieAddObject(ipLed);
    end
    
    ipWindow();
    
end

% Estimate channel only images from the difference between flash and no-flash
channelMeasVals = measVals;
channelRawLin = rawLin;
for i=2:(nChannels + 1)
     channelMeasVals(:,i,:) = max(channelMeasVals(:,i,:) - channelMeasVals(:,1,:),0);
     channelRawLin(:,:,:,i) = max(channelRawLin(:,:,:,i) - channelRawLin(:,:,:,1),0);
end

channelMeasVals = channelMeasVals/max(channelMeasVals(:));
channelRawLin = channelRawLin/max(channelRawLin(:));
%% Sanity check
%  Make sure that the difference of the estimate is the estimate of the
%  difference

tmp = channelMeasVals(:,2:end,:);

figure;
hold on; grid on; box on;
plot(tmp(:),measDiffVals(:),'.');


%% Ambient estimate
%  Estimate the ambient illuminant independently for every pixel, with some
%  spatial smoothing

% We crop a central portion of the image for computational convenience.

x = [324, 943];
y = [343, 755];
channelRawLin = channelRawLin(y(1):y(2),x(1):x(2),:,:);

hh = size(channelRawLin,1);
ww = size(channelRawLin,2);


alpha = 0.1;
beta = 0.1;

[ambientApproxWghts, hist] = ambientEstADMM( channelRawLin(:,:,:,1), channelRawLin(:,:,:,2:nChannels+1), flashNorm, alpha, beta,...
                                  'maxIter', 100, 'verbose', true );

% Plot the (spatially varying) estimated ambient spectrum

illuminantEst = flashNorm*reshape(ambientApproxWghts,[hh*ww, nChannels])';

figure;
hold on; grid on; box on;
plot(wave,illuminantEst(:,randi(size(illuminantEst,2),100,1)),'k');
plot(wave,cameraMat,'--');
xlabel('Wavelength, nm');

figure;
plot([hist.prRes hist.dualRes]);
xlabel('Iteration')


figure;
imagesc(ambientApproxWghts(:,:,1));


%% Complement estimate

desiredIll = 'D65';
ill = illuminantCreate(desiredIll,wave);
desiredSpectrum = illuminantGet(ill,'photons');
desiredSpectrum = desiredSpectrum/max(desiredSpectrum);

[ complWghts, slack, hist ] = complementEstADMM( desiredSpectrum, ambientApproxWghts, flashNorm, cameraMat, alpha, beta,...
                                                      'maxIter', 100, 'verbose', true, 'rescaleRho', false, 'flashMode', false );


                                                  


%% Computational synthesis
% Render an image as if captured with the complementary flash.

rendering = squeeze(channelRawLin(:,:,:,1));

% Ambient illuminant

sensor = sensorSet(sensor,'volts',rendering);
sensor = sensorSet(sensor,'name','Ambient');
ieAddObject(sensor);

ip = ipCreate;
ip = ipCompute(ip,sensor);
ip = ipSet(ip,'name','Ambient');
ieAddObject(ip);


for i=1:nChannels
    rendering = rendering + repmat(complWghts(:,:,i),[1 1 nFilters]).*squeeze(channelRawLin(:,:,:,i+1));
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


%% Render without spatial smoothing

beta = 0;

ambientApproxWghts = ambientEstADMM( channelRawLin(:,:,:,1), channelRawLin(:,:,:,2:nChannels+1), flashNorm, alpha, beta,...
                                  'maxIter', 100, 'verbose', true );


complWghts = complementEstADMM( desiredSpectrum, ambientApproxWghts, flashNorm, cameraMat, alpha, 10*beta,...
                                                      'maxIter', 100, 'verbose', true, 'rescaleRho', false, 'flashMode', false );


rendering = squeeze(channelRawLin(:,:,:,1));
for i=1:nChannels
    rendering = rendering + repmat(complWghts(:,:,i),[1 1 nFilters]).*squeeze(channelRawLin(:,:,:,i+1));
end

rendering = rendering/max(rendering(:));

sensor = sensorSet(sensor,'volts',rendering);
sensor = sensorSet(sensor,'name',sprintf('Computational: %s (no spatial)',desiredIll));
ieAddObject(sensor);
sensorWindow;


ip = ipCreate;
ip = ipCompute(ip,sensor);
ip = ipSet(ip,'name',sprintf('Computational: %s (no spatial)',desiredIll));
ieAddObject(ip);

ipWindow;

