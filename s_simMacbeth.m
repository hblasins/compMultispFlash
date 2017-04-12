% Simulate the appearance of a Macbeth chart under different illuminants
% and captured using different cameras. For each configuration simulate
% multispectral flash capture and analyze different error metrics.
%
% Copyright, Henryk Blasinski 2017.

close all;
clear all;
clc;

ieInit;

wave = 400:10:700;
nWaves = length(wave);

% Number of test patches.
% Here we use a Macbeth chart so
nSamples = 24; 

% Define camera responsivity functions to be tested
cameras = {'AptinaMT9M031','AptinaMT9M131',...
    'Canon1DMarkIII','Canon5DMarkII','Canon20D','Canon40D','Canon50D','Canon60D','Canon300D','Canon500D','Canon600D'...
    'HasselbladH2',...
    'NikonD1','NikonD3','NikonD3X','NikonD40','NikonD50','NikonD70','NikonD80','NikonD90','NikonD100','NikonD200',...
    'NikonD200IR','NikonD300s','NikonD700','NikonD5100',...
    'NokiaN900',...
    'OlympusE-PL2',...
    'PentaxK-5','PentaxQ',...
    'PhaseOne',...
    'PointGreyGrasshopper50S5C','PointGreyGrasshopper214S5C',...
    'SONYNEX-5N'}; 

nCameras = length(cameras);

% Define the chromaticity coordinates of the illuminants to be tested
xVals = linspace(0,1,20);
nXVals = length(xVals);
yVals = linspace(0,1,20);
nYVals = length(yVals);

% Define LED spectra
fName = fullfile(cmfRootPath,'Parameters','ximeaLights');
flashSpd = ieReadSpectra(fName,wave);
flashSpd = Energy2Quanta(wave,flashSpd);
flashSpdNorm = flashSpd/max(flashSpd(:));

nLEDs = size(flashSpd,2);
nLights = nLEDs + 1;

% Create a sensor
sensor = sensorCreate('bayer (rggb)');
sensor = sensorSet(sensor,'wave',wave);
sensor = sensorSet(sensor,'size',[100 100]);
sensor = sensorSet(sensor,'noise flag',0);
nFilters = 3;



% Initialize placeholders for different metrics


measPixelVals = zeros(nXVals,nYVals,nCameras,3,nLights,nSamples);
measPixelLedVals = zeros(nXVals,nYVals,nCameras,3,nLights,nSamples);


for cc=1:nCameras
    
    for xx=1:nXVals
        for yy=1:nYVals
            
            fName = fullfile(isetRootPath,'data','sensor','colorfilters',cameras{cc});
            cameraResp = ieReadColorFilter(wave,fName);
            cameraResp(isnan(cameraResp)) = 0;
            nFilters = size(cameraResp,2);
            sensor = sensorSet(sensor,'filter spectra',cameraResp);

            
            
            measVals = zeros(nFilters,nLights,nSamples);
            measLedVals = zeros(nFilters,nLights,nSamples);
            measValsLin = zeros(nFilters,nLights,nSamples);
            
            ill = illuminantCreate('D65',wave);
            ambientEnergy = 100*max(illuminantGet(ill,'energy'))*xy2Spectrum(xVals(xx),yVals(yy),wave(:));
            if sum(isnan(ambientEnergy)) > 0,
                continue;
            end
            
            %% Camera capture simulation
            %  Simulate the appearance of objects captured under ambient
            %  and under ambient+individual led.
            
            for i=1:nLights

                % Scene
                scene = sceneCreate('macbeth d65',3,wave);                
                if i==1
                    illEnergy = ambientEnergy;
                else
                    illEnergy = ambientEnergy + Quanta2Energy(wave(:),flashSpd(:,i-1))';
                end
                scene = sceneAdjustIlluminant(scene,illEnergy',0);
                
                % Oi
                oi = oiCompute(oiCreate,scene);
                switch i
                    case 1
                        name = '';
                        oiAmbient = oi;
                        oiAmbient = oiSet(oiAmbient,'name','Ambient');
                    otherwise
                        name = sprintf('LED %i',i-1);
                end
                
                oi = oiSet(oi,'name',['Ambient ' name]);
                
                if i>=2
                    oiLed = oiSet(oi,'photons',max(oiGet(oi,'photons') - oiGet(oiAmbient,'photons'),0));
                    oiLed = oiSet(oiLed,'name',name);
                end
                
                
                % Sensor
                sensor = sensorSetSizeToFOV(sensor,[sceneGet(scene,'fov horizontal') sceneGet(scene,'fov vertical')],scene,oi);
                sensor = sensorSet(sensor,'name',['Ambient ' name]);
                sensor = sensorCompute(sensor,oi);
                
                sensorSize = sensorGet(sensor,'size');
                cp = [0 sensorSize(1);sensorSize(2) sensorSize(1); sensorSize(2) 0; 0 0];
                
                
                [mVals, ~, ~, cp] = macbethSelect(sensor,0,1,cp);
                measVals(:,i,:) = cell2mat(cellfun(@(x) nanmean(x)',mVals,'UniformOutput',false));
                [cameraGain, cameraOffset] = sensorGainAndOffset(0.5,oi,sensor);
                
                measValsLin(:,i,:) = measVals(:,i,:)/cameraGain;
                                
                
                if i>=2
                    sensorLed = sensorSet(sensor,'name',name);
                    sensorLed = sensorCompute(sensorLed,oiLed);
                    [mVals, ~, ~, cp] = macbethSelect(sensorLed,0,1,cp);
                    measLedVals(:,i,:) = cell2mat(cellfun(@(x) nanmean(x)',mVals,'UniformOutput',false));
                    [cameraGainLed, cameraOffsetLed] = sensorGainAndOffset(0.5,oiLed,sensorLed);
                    measLedVals(:,i,:) = measLedVals(:,i,:)/cameraGainLed;
                end
                
                
                
            end
            
            measPixelLedVals(xx,yy,cc,:,:,:) = measLedVals;
            measPixelVals(xx,yy,cc,:,:,:) = measValsLin;
            
        end
    end
end

fName = fullfile(compMultispFlashRootPath,'Results','simGlobalMacbeth.mat');
save(fName);


%{

            
            %% Estimate LED only images
            %  This is done by subtracting the ambient image from ambient+flash images
            %  Note that we normalize the estimated result.
            
            measValsLed = max(measValsLin(:,2:end,:) - repmat(measValsLin(:,1,:),[1, nLEDs, 1]),0);
            measValsLed = measValsLed/max(measValsLed(:));
            
            measValsLin = measValsLin/max(measValsLin(:));
            
            measPixelVals(xx,yy,cc,:,:,:) = measValsLin;
            
            %% Estimate LED weights as seen through the camera
            %  that best approximate the ambient illuminant image. This is
            %  the estimation algorithm idea we propose in the paper.
            
            cvx_begin
                variables ambientApproxWghts(nLEDs,1)
            
                approx = 0;
                for i=1:nLEDs
                    approx = approx + squeeze(measValsLed(:,i,:))*ambientApproxWghts(i);
                end
            
                minimize sum(norms(squeeze(measValsLin(:,1,:)) - approx,2,1))
                subject to
                    flashSpdNorm*ambientApproxWghts >= 0
            cvx_end
            
            
            %% Estimate LED weights measured spectrally
            %  Here we assume that we know the true illuminant spectrum and
            %  we try to optimize the led weights to best approximate this
            %  spectrum.
            
            illPhotons = Energy2Quanta(wave,illEnergy);
            illPhotons = illPhotons/max(illPhotons);
            cvx_begin
                variables ambientSpectralApproxWghts(nLEDs,1)
                minimize norm(illPhotons - flashSpdNorm*ambientSpectralApproxWghts)
                subject to
                    ambientSpectralApproxWghts >= 0
            cvx_end
            
            
            %% Estimate complementary flash to color balance to the desired illuminant
            %  for different desired illuminants an when the weights are
            %  constrained to [0,1] (this is applicable when one wants to
            %  generate an actual flash image), or are constrained to be
            %  non-negative (a computational flash).
            
            for dd=1:nDesiredIlluminants
                
                desiredIll = desiredIlluminants{dd};
                
                cvx_begin
                    variables ambientComplementWghtsConstr(nLEDs,1) scale(1,1)

                    minimize norm(cameraResp'*(flashSpdNorm*(ambientApproxWghts + ambientComplementWghtsConstr)) - scale*cameraResp'*desiredIll)
                    subject to
                        1>= ambientComplementWghtsConstr >= 0
                        scale >= 0
                cvx_end

                cvx_begin
                    variables ambientComplementWghtsUnc(nLEDs,1) scale(1,1)

                    minimize norm(cameraResp'*(flashSpdNorm*(ambientApproxWghts + ambientComplementWghtsUnc)) - scale*cameraResp'*desiredIll)
                    subject to
                        ambientComplementWghtsUnc >= 0
                        scale >= 0
                cvx_end
            
                complementWeightsConstr(xx,yy,cc,dd,:) = ambientComplementWghtsConstr;
                complementWeightsUnc(xx,yy,cc,dd,:) = ambientComplementWghtsUnc;
            
                ambientComplementConstr = flashSpdNorm*(ambientApproxWghts + ambientComplementWghtsConstr);
                ambientComplementConstr = ambientComplementConstr/max(ambientComplementConstr);
            
                ambientComplementUnc = flashSpdNorm*(ambientApproxWghts + ambientComplementWghtsUnc);
                ambientComplementUnc = ambientComplementUnc/max(ambientComplementUnc);
                
                ambientComplementConstrXYZ(xx,yy,cc,dd,:) = ieXYZFromPhotons(ambientComplementConstr,wave);
                ambientComplementUncXYZ(xx,yy,cc,dd,:) = ieXYZFromPhotons(ambientComplementUnc,wave);
                
            end
            
            meas = squeeze(measValsLin(:,1,:));
            approxPixelVals(xx,yy,cc,:,:,:) = approx;
            
            approxWeights(xx,yy,cc,:) = ambientApproxWghts;
            spectralApproxWeights(xx,yy,cc,:) = ambientSpectralApproxWghts;            
            
            ambientEst = flashSpd*ambientApproxWghts;
            ambientEst = ambientEst/max(ambientEst);
            
            ambientSpectrum = Energy2Quanta(wave,ambientEnergy);
            ambientSpectrum = ambientSpectrum/max(ambientSpectrum);

            estXYZ(xx,yy,cc,:) = ieXYZFromPhotons(ambientEst,wave);     
            ambientXYZ(xx,yy,cc,:) = ieXYZFromPhotons(ambientSpectrum,wave);
            
            
            
        end
    end
end

fName = fullfile(slRootPath,'Results','evalApprox.mat');
save(fName);



%% Plot chromaticity results
close all;
clear all;
clc;


fName = fullfile(compMultispFlashRootPath,'Results','evalApprox.mat');
load(fName);

chromaticityPlot;
hold on; grid on; box on;
for xx=1:nXVals
    for yy=1:nYVals
        for cc=1:nCameras
            estXY = estXYZ(xx,yy,cc,1:2)/sum(estXYZ(xx,yy,cc,:));
            plot(estXY(1),estXY(2),'x');
        end
        
        
        cc = 1;
        ambXY = ambientXYZ(xx,yy,cc,1:2)/sum(ambientXYZ(xx,yy,cc,:));
        
        plot(ambXY(1),ambXY(2),'ro','lineWidth',2,'markerSize',5);
    end
    
end

for i=1:nLEDs
    xyz = ieXYZFromPhotons(flashSpd(:,i),wave);
    plot(xyz(1)/sum(xyz),xyz(2)/sum(xyz),'k+','lineWidth',2,'markerSize',10);
end

%% RGB approximation errors

figure;
hold on; grid on; box on;
for cc=1:nCameras
        
    for xx=1:nXVals
        for yy=1:nYVals
            
                measured = squeeze(measPixelVals(xx,yy,cc,:,1,:));
                approximated = squeeze(approxPixelVals(xx,yy,cc,:,1,:));
                plot(measured(:),approximated(:),'k.');
        end
    end
    
end

xlabel('Measured RGB');
ylabel('LED approximated RGB');

%% Complement to desired illuminant error

figure;
for xx=5 %1:nXVals
    for yy=5 %1:nYVals
        
        tmp = ambientXYZ(xx,yy,:,:);
        if sum(tmp(:)) == 0, continue; end;
        
        for dd=1:nDesiredIlluminants
            
            subplot(2,2,dd);
            chromaticityPlot([],'white',500,false);
            hold on; grid on; box on;
            
            % Plot the ambient illuminant
            for cc=1:nCameras
                estXY = estXYZ(xx,yy,cc,1:2)/sum(estXYZ(xx,yy,cc,:));
                plot(estXY(1),estXY(2),'x');
            end
            
            
            cc = 1;
            ambXY = ambientXYZ(xx,yy,cc,1:2)/sum(ambientXYZ(xx,yy,cc,:));
            
            plot(ambXY(1),ambXY(2),'ro','lineWidth',2,'markerSize',5);
            
            
            % Plot the complemented illuminant (Constrained)
            for cc=1:nCameras
                complXY = ambientComplementConstrXYZ(xx,yy,cc,dd,1:2)/sum(ambientComplementConstrXYZ(xx,yy,cc,dd,:));
                plot(complXY(1),complXY(2),'x');
            end
            
            % Plot the complemented illuminant (Unconstrained)
            for cc=1:nCameras
                complXY = ambientComplementUncXYZ(xx,yy,cc,dd,1:2)/sum(ambientComplementUncXYZ(xx,yy,cc,dd,:));
                plot(complXY(1),complXY(2),'d');
            end
            
            % Plot the desired illuminant.
            desXYZ = ieXYZFromPhotons(desiredIlluminants{dd},wave);
            desXY = desXYZ(1:2)/sum(desXYZ);
            
            plot(desXY(1),desXY(2),'m+','lineWidth',3,'markerSize',15);
            
            
            % Plot the individual LED spectra.
            for i=1:nLEDs
                xyz = ieXYZFromPhotons(flashSpd(:,i),wave);
                plot(xyz(1)/sum(xyz),xyz(2)/sum(xyz),'k+','lineWidth',3,'markerSize',15);
            end
            
            
            
            title(sprintf('Desired illuminant %i',dd));
            drawnow;
        end
    end
    
end

%}

