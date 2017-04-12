function [ sensors, ips, rawData, demosaicedData, patchData ] = getSensorData( sensor, oi)

cp = [4 81;97 81;97 20;4 19]/100;
cp = [cp(:,1)*sensorGet(sensor,'col'), cp(:,2)*sensorGet(sensor,'row')];
nSamples = 24;

nOis = length(oi);
sensors = cell(nOis,1);
ips = cell(nOis,1);

patchData = zeros(3,nOis,nSamples);
rawData = zeros(sensorGet(sensor,'row'),sensorGet(sensor,'col'),nOis);
demosaicedData = zeros(sensorGet(sensor,'row'),sensorGet(sensor,'col'),3,nOis);


%% Camera simulation

for i=1:nOis
   
    
    
    ieAddObject(oi{i});
    oiWindow();
    
    sensors{i} = sensorSet(sensor,'pixel size', oiGet(oi{i},'spatial resolution'));
    sensors{i} = sensorSet(sensors{i},'size',oiGet(oi{i},'size'));
    sensors{i} = pixelCenterFillPD(sensors{i},1);
    sensors{i} = sensorSet(sensors{i},'name',oiGet(oi{i},'name'));
    sensors{i} = sensorCompute(sensors{i},oi{i});
    
    vcAddObject(sensors{i});
    
    [mVals, ~, ~, cp] = macbethSelect(sensors{i},1,1,cp);
    patchData(:,i,:) = cell2mat(cellfun(@(x) nanmean(x)',mVals,'UniformOutput',false));
        
    [cg, co] = sensorGainAndOffset(0.5,oi{i},sensors{i});
    rawData(:,:,i) = sensorGet(sensors{i},'volts')/cg;
    patchData(:,i,:) = patchData(:,i,:)/cg;
        
     
    sensorWindow();
    
    
    ips{i} = ipCreate();
    ips{i} = ipCompute(ips{i},sensors{i});
    ips{i} = ipSet(ips{i},'name',oiGet(oi{i},'name'));
    ieAddObject(ips{i});
    
    demosaicedData(:,:,:,i) = ipGet(ips{i},'sensor channels');
    
    ipWindow();
    
end

%{

% Estimate channel only images from the difference between flash and no-flash
channelMeasVals = measVals;
channelRawLin = rawLin;
for i=2:(nChannels + 1)
     channelMeasVals(:,i,:) = max(channelMeasVals(:,i,:) - channelMeasVals(:,1,:),0);
     channelRawLin(:,:,i) = max(channelRawLin(:,:,i) - channelRawLin(:,:,1),0);
end
channelMeasVals = channelMeasVals/max(channelMeasVals(:));
channelRawLin = channelRawLin/max(channelRawLin(:));


patches.ambient = squeeze(channelMeasVals(:,1,:));
patches.leds = channelMeasVals(:,2:end,:);

images.ambient = squeeze(channelRawLin(:,:,1));
images.ambientDemosaiced = squeeze(rawDemosaiced(:,:,:,1));
images.ambientDemosaicedVectorized = reshape(images.ambientDemosaiced,[size(images.ambient,1)*size(images.ambient,2),3])';

images.leds = squeeze(channelRawLin(:,:,2:end));
images.ledsDemosaiced = squeeze(rawDemosaiced(:,:,:,2:end));
images.ledsDemosaicedVectorized = reshape(images.ledsDemosaiced,[size(images.ambient,1)*size(images.ambient,2),3, 7]);
images.ledsDemosaicedVectorized = permute(images.ledsDemosaicedVectorized, [2, 3, 1]);

% data.ambientDemosaiced = 
% data.ledsDemosaiced = 
% data.ledsDemosaiced = permute(data.ledsDemosaiced,[2, 3, 1]);
%}
end

