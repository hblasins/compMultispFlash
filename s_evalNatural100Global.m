close all;
clear all;
clc;

ieInit;
set(groot,'defaultAxesColorOrder',[1 0 0; 0 1 0; 0 0 1]);


wave = 400:5:700; wave = wave(:);
nWaves = length(wave);

fName = fullfile(slRootPath,'Parameters','ximeaLights');
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
        
                
        for i=1:(nChannels + 1)
            
            if i==1
                name = 'Ambient';
                
                scene = sceneCreateNatural100();
                scene = sceneSet(scene,'wavelength',wave);
                illSpd = spd;
            else
                name = sprintf('Ambient + LED%i',i-1);
                illSpd = spd + flashNorm(:,i-1);
            end
            
            illSpd = illSpd*10^15;
            
            scene = sceneAdjustIlluminant(scene,Quanta2Energy(wave,illSpd),0);
            % ieAddObject(scene);
            % sceneWindow;
            
            oi = oiCompute(oiCreate,scene);
            oi = oiSet(oi,'name',name);
            % ieAddObject(oi);
            % oiWindow();
            
            if i==1
                sensor = sensorSet(sensor,'pixel size', oiGet(oi,'spatial resolution'));
                sensor = sensorSetSizeToFOV(sensor,[sceneGet(scene,'hfov') sceneGet(scene,'vfov')],scene,oi);
                sensor = pixelCenterFillPD(sensor,1);
            end
            sensor = sensorSet(sensor,'name',name);
            sensor = sensorCompute(sensor,oi);
            
            [cg, co] = sensorGainAndOffsetV2(0.5,oi,sensor);
            
            if i==1
                tmp = sensorGet(sensor,'volts');
                rawLin = zeros([size(tmp),nChannels + 1]);
            end
            
            rawLin(:,:,i) = sensorGet(sensor,'volts')/cg;
            
            % vcAddObject(sensor);
            % sensorWindow();
            
            
            ip = ipCreate();
            ip = ipCompute(ip,sensor);
            ip = ipSet(ip,'name',name);
            % ieAddObject(ip);
            % ipWindow();
            
        end
        
        
        % Estimate channel only images from the difference between flash and no-flash
        channelRawLin = rawLin;
        for i=2:(nChannels + 1)
            channelRawLin(:,:,i) = max(channelRawLin(:,:,i) - channelRawLin(:,:,1),0);
        end
        
        channelRawLin = channelRawLin/max(channelRawLin(:));
        
        %% Ambient estimate
        
        lambda = 0;
        R = [eye(nWaves-1) zeros(nWaves-1,1)] - [zeros(nWaves-1,1) eye(nWaves-1)];
        
        ambientAppearance = channelRawLin(:,:,1);
        cvx_begin
        variables ambientApproxWghts(nChannels,1)
        
        approx = 0;
        for i=1:nChannels
            approx = approx + squeeze(channelRawLin(:,:,i+1))*ambientApproxWghts(i);
        end
        
        minimize sum(norms(ambientAppearance - approx,2,1)) + lambda * norm(R*flashNorm*ambientApproxWghts,2)
        subject to
        flashNorm*ambientApproxWghts >= 0
        cvx_end
        
        err = sqrt(mean((ambientAppearance(:) - approx(:)).^2));
        
        figure(resFig);
        subplot(2,3,2);
        hold on; grid on; box on;
        plot(approx(:),ambientAppearance(:),'.r');
        xlabel('Weighed sum');
        ylabel('Simulation');
        title(sprintf('Approx. RMSE %3.f',err));
        
        % Plot the quality of the approximation in camera RGB space
        %{
        figure;
        hold on; grid on; box on;
        plot(ambientAppearance',approx','o');
        xlabel('Ambient appearance');
        ylabel('Approximated');
        title('RGB space');
        
        % Plot the estimated ambient spectrum
        
        illuminantEst = flashNorm*ambientApproxWghts;
        
        figure;
        hold on; grid on; box on;
        plot(wave,illuminantEst);
        legend('Estimated');
        plot(wave,cameraResp,'--');
        %}
        
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
            
            cvx_begin
            variables flashCompWeightsUnc(nChannels,1) en(1,1)
            minimize norm( cameraResp'*(flashNorm*(flashCompWeightsUnc+ambientApproxWghts))- en*camDes)
            subject to
            flashCompWeightsUnc >= 0
            en >= 0
            cvx_end
            
            cvx_begin
            variables flashCompWeights(nChannels,1) en(1,1)
            minimize norm( cameraResp'*(flashNorm*(flashCompWeights+ambientApproxWghts))- en*camDes)
            subject to
            1 >= flashCompWeights >= 0
            en >= 0
            cvx_end
            
            
            rendering = squeeze(channelRawLin(:,:,1));
            for i=1:nChannels
                rendering = rendering + flashCompWeights(i)*squeeze(channelRawLin(:,:,i+1));
            end
            rendering = rendering/max(rendering(:));
            
            renderingUnc = squeeze(channelRawLin(:,:,1));
            for i=1:nChannels
                renderingUnc = renderingUnc + flashCompWeightsUnc(i)*squeeze(channelRawLin(:,:,i+1));
            end
            renderingUnc = renderingUnc/max(renderingUnc(:));
            
            figure(resFig);
            flashUnc = flashNorm*(flashCompWeightsUnc+ambientApproxWghts);
            flashUncXYZ = ieXYZFromPhotons(flashUnc,wave);
            flashUncxy = flashUncXYZ(1:2)/sum(flashUncXYZ);
            plot(flashUncxy(1),flashUncxy(2),'mx','lineWidth',2,'markerSize',10);
            
            flash = flashNorm*(flashCompWeights+ambientApproxWghts);
            flashXYZ = ieXYZFromPhotons(flash,wave);
            flashxy = flashXYZ(1:2)/sum(flashXYZ);
            plot(flashxy(1),flashxy(2),'c+','lineWidth',2,'markerSize',10);
            
            %% Simulate under teh desired illuminant
            
            scene = sceneAdjustIlluminant(scene,illuminantGet(ill,'energy'),0);
            scene = sceneSet(scene,'name',sprintf('Desired %iK',desiredTemp));
            
            oi = oiCompute(oi,scene);
            oi = oiSet(oi,'name',sprintf('Desired %iK',desiredTemp));
            % ieAddObject(oi);
            % oiWindow();
            
            sensor = sensorCompute(sensor,oi);
            sensor = sensorSet(sensor,'name', sprintf('Desired %iK',desiredTemp));
            [cg, co] = sensorGainAndOffsetV2(0.5,oi,sensor);
            
            desiredSimulation = sensorGet(sensor,'volts')/cg;
            desiredSimulation = desiredSimulation/max(desiredSimulation(:));
            
            % ieAddObject(sensor);
            % sensorWindow();
            
            ip = ipCompute(ip,sensor);
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
            
            legend(sprintf('Flash RMSE: %.3f',err),sprintf('Unc. RMSE: %.3f',errUnc),'location','northWest');
            xlabel('Weighted sum');
            ylabel('Simulation');
            title(sprintf('Desired: %iK',desiredTemp(dd)));
            
            drawnow;
            
        end
    end
end
