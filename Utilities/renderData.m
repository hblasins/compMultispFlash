function [measurement, reference] = renderData(wave, sensor, illuminantSpectrum, ledSpectra, varargin)

p = inputParser;
p.addOptional('target','Natural100');
p.addOptional('compact',false);
p.parse(varargin{:});
inputs = p.Results;

nLeds = size(ledSpectra,2);
nChannels = nLeds + 1;

measurement = [];
reference = [];
 
ips = cell(nChannels,1);
sensors = cell(nChannels,1);
scenes = cell(nChannels,1);
ois = cell(nChannels,1);
 
 for i=1:(nChannels)
     
     %% Scene
     switch inputs.target
         case 'Natural100'
             scenes{i} = sceneCreateNatural100();
             patchH = 10;
             patchW = 11;
             
         case 'Macbeth'
             scenes{i} = sceneCreate('macbeth',32);
             patchH = 4;
             patchW = 6;
     end
     scenes{i} = sceneSet(scenes{i},'wavelength',wave);
     
     
     if i==1
         name = 'Ambient';
         illSpd = illuminantSpectrum;
     else
         name = sprintf('Ambient + LED%i',i-1);
         illSpd = illuminantSpectrum + ledSpectra(:,i-1);
     end
     scenes{i} = sceneSet(scenes{i},'name',name);
     illSpd = illSpd*10^15;
     
     scenes{i} = sceneAdjustIlluminant(scenes{i},Quanta2Energy(wave,illSpd),0);
     
     ieAddObject(scenes{i});
     sceneWindow;
     
     %% OI
     
     ois{i} = oiCompute(oiCreate,scenes{i});
     ois{i} = oiSet(ois{i},'name',name);
     
     ieAddObject(ois{i});
     oiWindow();
     
     %% Sensor
     sensors{i} = sensorSet(sensor,'pixel size', oiGet(ois{i},'spatial resolution'));
     sensors{i} = sensorSetSizeToFOV(sensors{i},[sceneGet(scenes{i},'hfov') sceneGet(scenes{i},'vfov')],scenes{i},ois{i});
     sensors{i} = pixelCenterFillPD(sensors{i},1);
     sensors{i} = sensorSet(sensors{i},'name',name);
     sensors{i} = sensorCompute(sensors{i},ois{i});
     
     % cg = sensorGainAndOffset(0.5,ois{i},sensors{i});
     cg = sensorGet(sensors{i},'exposure time');
     
     if i==1
         tmp = sensorGet(sensors{i},'volts');
         rawData = zeros([size(tmp),nChannels]);
         demosaicedData = zeros([size(tmp),3,nChannels]);
         sampleData = zeros(3,nChannels,patchH*patchW);
     end
     
     rawData(:,:,i) = sensorGet(sensors{i},'volts')/cg;
     
     rows = sensorGet(sensors{i},'rows');
     cols = sensorGet(sensors{i},'cols');
     
     cp = [1 rows;
           cols rows;
           cols 1;
           1 1];
     
     tmp = chartSelect(sensors{i},0,1,patchH,patchW,cp);
     tmp = cell2mat(cellfun(@nanmean,tmp,'UniformOutput',false)')';
     sampleData(:,i,:) = tmp/cg;
     
     vcAddObject(sensors{i});
     sensorWindow();
     
     %% ISP
     
     ips{i} = ipCreate();
     ips{i} = ipCompute(ips{i},sensors{i});
     ips{i} = ipSet(ips{i},'name',name);
     
     demosaicedData(:,:,:,i) = ipGet(ips{i},'sensor channels')/cg;
     
     ieAddObject(ips{i});
     ipWindow();
     
 end

if inputs.compact == false

    measurement.sensors = sensors;
    measurement.ips = ips;
    measurement.ois = ois;
    measurement.scenes = scenes;

    measurement.demosaiced.data = demosaicedData/max(demosaicedData(:));
    measurement.raw.data = rawData/max(rawData(:));

    measurement.demosaiced.ambient = measurement.demosaiced.data(:,:,:,1);
    measurement.demosaiced.led = max(measurement.demosaiced.data(:,:,:,2:end) - repmat(measurement.demosaiced.ambient,[1 1 1, nLeds]),0);

    measurement.raw.ambient = measurement.raw.data(:,:,1);
    measurement.raw.led = max(measurement.raw.data(:,:,2:end) - repmat(measurement.raw.ambient,[1 1 nLeds]),0);
    
    measurement.vectorized.ambient = reshape(measurement.demosaiced.ambient,[sensorGet(sensors{1},'rows')*sensorGet(sensors{1},'cols'), 3])';
    measurement.vectorized.led = permute(reshape(measurement.demosaiced.led,[sensorGet(sensors{1},'rows')*sensorGet(sensors{1},'cols'), 3, nLeds]),[2, 3, 1]);
    
end

measurement.patch.data = sampleData/max(sampleData(:));
measurement.patch.ambient = measurement.patch.data(:,1,:);
measurement.patch.led = max(measurement.patch.data(:,2:end,:) - repmat(measurement.patch.ambient,[1 nLeds,1]),0);



 
 