% This script evaluates the computational multispectral flash using the
% Natural100 patches. It produces an image under an illuminant with a
% specified CIE xy coordinates, estimates the image under the complementary
% light and generates the image of the scene image under the desired 
% illuminant to compare.
%
% Copyright, Henryk Blasinski 2017.

close all;
clear all;
clc;

ieInit;
set(groot,'defaultAxesColorOrder',[1 0 0; 0 1 0; 0 0 1]);


wave = 400:5:700; wave = wave(:);
nWaves = length(wave);

fName = fullfile(cmfRootPath,'Parameters','ximeaLights');
flashSpd = ieReadSpectra(fName,wave);
flashSpd = Energy2Quanta(wave,flashSpd);
flashNorm = flashSpd/max(flashSpd(:));

nChannels = size(flashSpd,2);

sensor = sensorCreate('bayer (rggb)');
sensor = sensorSet(sensor,'wave',wave);
sensor = sensorSet(sensor,'noise flag',0);

cameraResp = sensorGet(sensor,'filter transmissivities');


desiredTemp = [10000, 6500, 4000];


%% Camera simulation

for xx=0:0.2:1
    for yy=0:0.2:1
        
        [ spd ] = xy2Spectrum( xx, yy, wave );
        if sum(isnan(spd)) > 0
            continue,
        end
        
        resFig = figure;
        subplot(2,3,1);
        hold on; grid on; box on;
        chromaticityPlot([],'gray',256,false);
        
        spdXYZ = ieXYZFromPhotons(spd,wave);
        spdxy = spdXYZ(1:2)/sum(spdXYZ);
        plot(spdxy(1),spdxy(2),'rd','lineWidth',2,'markerSize',10);
        
        for i=1:nChannels
            xyz = ieXYZFromPhotons(flashSpd(:,i),wave);
            plot(xyz(1)/sum(xyz),xyz(2)/sum(xyz),'k+','lineWidth',2,'markerSize',10);
        end
        
       
        
        measurement = renderData(wave, sensor, spd, flashNorm);
        
        channelRawLin = measurement.raw.data;
        
        %% Ambient estimate
        
        [ ambientEst, ambientWghts, ambientPredictions ] = globalAmbientEst( measurement.patch.ambient, measurement.patch.led, flashNorm, 'alpha',0.1 );
        
        
        err = sqrt(mean((measurement.patch.ambient(:) - ambientPredictions(:)).^2));
        
        figure(resFig);
        subplot(2,3,2);
        hold on; grid on; box on;
        plot(ambientPredictions(:)',squeeze(measurement.patch.ambient(:))','.r');
        xlabel('Weighed sum');
        ylabel('Simulation');
        title(sprintf('Approx. RMSE %.3f',err));
        xlim([0 1]);
        ylim([0 1]);
        
        
        
        % Plot the estimated ambient spectrum
        figure(resFig);
        subplot(2,3,3);
        hold on; grid on; box on;
        plot(wave,ambientEst/max(ambientEst(:)),'r');
        plot(wave,spd/max(spd(:)),'g');
        legend('Estimated','True');
        
        
        %% Render under the desired illuminant
        
        for dd=1:length(desiredTemp)
            
            ill = illuminantCreate('blackbody',wave,desiredTemp(dd));
            desiredSpectrum = illuminantGet(ill,'photons');
            
            desXYZ = ieXYZFromPhotons(desiredSpectrum,wave);
            desxy = desXYZ(1:2)/sum(desXYZ);
            figure(resFig);
            subplot(2,3,1);
            plot(desxy(1),desxy(2),'go','lineWidth',2,'markerSize',10);
            
            camDes = cameraResp'*desiredSpectrum;
            camDes = camDes/max(camDes);
            
            [ flashEst, flashWghts ] = globalComplementEst( desiredSpectrum, ambientEst, flashNorm, cameraResp, 'flashMode', true);
            [ flashUncEst, flashUncWghts ] = globalComplementEst( desiredSpectrum, ambientEst, flashNorm, cameraResp, 'flashMode', false);
            
            
            flash = renderFlashImage(measurement.raw.ambient,measurement.raw.led,flashWghts,sensor);
            flashUnc = renderFlashImage(measurement.raw.ambient,measurement.raw.led,flashUncWghts,sensor);
            
            rendering = sensorGet(flash,'volts');
            rendering = rendering/max(rendering(:));
            renderingUnc = sensorGet(flashUnc,'volts');
            renderingUnc = renderingUnc/max(renderingUnc(:));
            
            figure(resFig);
            flashUncXYZ = ieXYZFromPhotons(flashEst,wave);
            flashUncxy = flashUncXYZ(1:2)/sum(flashUncXYZ);
            plot(flashUncxy(1),flashUncxy(2),'mx','lineWidth',2,'markerSize',10);
            
            flashXYZ = ieXYZFromPhotons(flashUncEst,wave);
            flashxy = flashXYZ(1:2)/sum(flashXYZ);
            plot(flashxy(1),flashxy(2),'c+','lineWidth',2,'markerSize',10);
            
            %% Simulate under the desired illuminant
            
            scene = sceneAdjustIlluminant(measurement.scenes{1},illuminantGet(ill,'energy'),0);
            scene = sceneSet(scene,'name',sprintf('Desired %iK',desiredTemp));
            
            oi = oiCompute(oiCreate,scene);
            oi = oiSet(oi,'name',sprintf('Desired %iK',desiredTemp));
            % ieAddObject(oi);
            % oiWindow();
            
            sensor = sensorCompute(measurement.sensors{1},oi);
            sensor = sensorSet(sensor,'name', sprintf('Desired %iK',desiredTemp));
            [cg, co] = sensorGainAndOffset(0.5,oi,sensor);
            
            desiredSimulation = sensorGet(sensor,'volts')/cg;
            desiredSimulation = desiredSimulation/max(desiredSimulation(:));
            
            % ieAddObject(sensor);
            % sensorWindow();
            
            ip = ipCompute(ipCreate,sensor);
            ip = ipSet(ip,'name',sprintf('Desired %iK',desiredTemp));
            % ieAddObject(ip);
            % ipWindow();
            
            
            figure(resFig);
            subplot(2,3,dd+3);
            hold on; grid on; box on;
            plot(rendering(:),desiredSimulation(:),'c.');
            plot(renderingUnc(:),desiredSimulation(:),'g.');
            axis square;
            
            err = sqrt(mean((rendering(:) - desiredSimulation(:)).^2));
            errUnc = sqrt(mean((renderingUnc(:) - desiredSimulation(:)).^2));
            
            legend(sprintf('Flash %.3f',err),sprintf('Unc. %.3f',errUnc),'location','northWest');
            xlabel('Weighted sum');
            ylabel('Simulation');
            title(sprintf('Desired: %iK',desiredTemp(dd)));
            
            drawnow;
            
        end
    end
end
