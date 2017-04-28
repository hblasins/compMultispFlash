% Simulate the appearance of a Natural100 chart under ambient illuminant
% and a number of narrowband LED lights. Use ambient illuminant with different
% CIE xy chromaticity coordinates.
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
       
% Define the chromaticity coordinates of the illuminant
xVec = linspace(0,1,10);
yVec = linspace(0,1,10);
nX = length(xVec);
nY = length(yVec);
       
% Define camera models

cameras = {'AptinaMT9M031','AptinaMT9M131',...
           'Canon1DMarkIII','Canon5DMarkII','Canon20D','Canon40D','Canon50D','Canon60D','Canon300D','Canon500D','Canon600D',...
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

%% Camera simulation

measurement = cell(nX,nY,nCameras);
for xx=1:nX
    for yy=1:nY
        
        [ spd ] = xy2Spectrum( xVec(xx), yVec(yy), wave );
        if sum(isnan(spd)) > 0
            continue,
        end
        
        for cam=1:nCameras
            
            sensor = sensorCreate('bayer (rggb)');
            sensor = sensorSet(sensor,'wave',wave);
            sensor = sensorSet(sensor,'noise flag',0);
            
            fName = fullfile(isetRootPath,'data','sensor','colorfilters',cameras{cam});
            cameraResp = ieReadSpectra(fName,wave);
            cameraResp(isnan(cameraResp)) = 0;
            sensor = sensorSet(sensor,'filter transmissivities',cameraResp);
            
            measurement{xx,yy,cam} = renderData(wave, sensor, spd, flashNorm, 'compact',true);
            
            %{
            
            rawLin = zeros(10,11,3,nChannels+1);
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
                
                [cg, co] = sensorGainAndOffset(0.5,oi,sensor);
                
                sz = sensorGet(sensor,'size');
                cp = [1 sz(1); sz(2) sz(1); sz(2) 1; 1 1];
                
                vals = chartSelect(sensor,0,1,10,11,cp);
                tmp = cellfun(@nanmean,vals,'UniformOutput', false)';
                tmp = cell2mat(tmp);
                rawLin(:,:,:,i) = reshape(tmp/cg,[10 11 3]);
                
                % vcAddObject(sensor);
                % sensorWindow();
                
                
                ip = ipCreate();
                ip = ipCompute(ip,sensor);
                ip = ipSet(ip,'name',name);
                % ieAddObject(ip);
                % ipWindow();
                
            end
            
            % Estimate channel only images from the difference between flash and no-flash
            
            for i=2:(nChannels + 1)
                rawLin(:,:,:,i) = max(rawLin(:,:,:,i) - rawLin(:,:,:,1),0);
            end
            rawLin = rawLin/max(rawLin(:));
            channelRawLin(xx,yy,:,:,:,:,cam) = rawLin;
            
            %}
            
            % Clear ISET data
            vcDeleteSomeObjects('scene',1:length(vcGetObjects('scene'))); 
            vcDeleteSomeObjects('oi',1:length(vcGetObjects('oi'))); 
            vcDeleteSomeObjects('sensor',1:length(vcGetObjects('sensor'))); 
            vcDeleteSomeObjects('ip',1:length(vcGetObjects('ip'))); 
            
        end

    end
end

%%
fName = fullfile(cmfRootPath,'Results','simNatural100V2.mat');
save(fName);
