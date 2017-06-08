%%% RenderToolbox3 Copyright (c) 2012-2013 The RenderToolbox3 Team.
%%% About Us://github.com/DavidBrainard/RenderToolbox3/wiki/About-Us
%%% RenderToolbox3 is released under the MIT License.  See LICENSE.txt.

%% Scene description

% Henryk Blasinski
close all;
clear all;
clc;

ieInit;

%% Simulation parameters

model = 'Acropora';

filmDiag = 20;
targetDistance = 1000;
depth = 10000;


chlConc = 10.000;
cdomConc = 0.000;
smallPartConc = 0.000;
largePartConc = 0.000;

nSamples = 128;


%%
switch (model)
    
    case 'Pocillopora'
        parentSceneFile = fullfile('Models',model,'model.obj');
        objectSize = 45;
        lookAt = [0; 0; 1];
        viewAxis = [1 0 0];
        hints.recipeName = 'Coral-Pocillopora';
        
    case 'Acropora'
        parentSceneFile = fullfile('Models',model,'model.obj');
        objectSize = 40;
        lookAt = [0; 0; 10];
        viewAxis = [0.5 -1 0.5];
        hints.recipeName = 'Coral-Acropora';
        
        
    case 'Table'
        parentSceneFile = fullfile('Models',model,'model.obj');
        objectSize = 40;
        lookAt = [0; 0; 4];
        viewAxis = [0.5; -1; 0.5];
        hints.recipeName = 'Coral-Table';
        
end



% Renderer options.
hints.imageWidth = 640;
hints.imageHeight = 480;
hints.renderer = 'PBRT'; % We're only using PBRT right now
hints.copyResources = 1;
hints.batchRenderStrategy = RtbAssimpStrategy(hints);
hints.batchRenderStrategy.remodelPerConditionAfterFunction = @MakeCoralMexximpRemodeller;
hints.batchRenderStrategy.converter.remodelAfterMappingsFunction = @MakeCoralPBRTRemodeller;

% Change the docker container
hints.batchRenderStrategy.renderer.pbrt.dockerImage = 'vistalab/pbrt-v2-spectral';

resourceFolder = rtbWorkingFolder('folderName','resources',...
    'rendererSpecific',false,...
    'hints',hints);

conditionsFile = fullfile(resourceFolder,'CoralConditions.txt');

root = rtbWorkingFolder('hints',hints);



[wavesDaylight, daylight] = rtbReadSpectrum('D65.spd');
daylight = daylight/max(daylight(:));



fName = fullfile('..','Parameters','ximeaLights.mat');
[ledSpectra, wavesLed] = ieReadSpectra(fName,[]);
ledSpectra = Energy2Quanta(wavesLed,ledSpectra);
ledSpectra = ledSpectra/max(ledSpectra(:));

daylight = interp1(wavesDaylight,daylight,wavesLed);

ledSpectra = 10000000*[zeros(length(wavesLed),1), ledSpectra, daylight];

nLeds = size(ledSpectra,2);


%% Save data to resources folder

for cC = chlConc
    for cdom = cdomConc
        
        [sig_a, waves] = createAbsorptionCurve(cC,cdom);
        absorptionFile = fullfile(resourceFolder, sprintf('abs_%.3f_%.3f.spd',cC,cdom));
        rtbWriteSpectrumFile(waves, sig_a, absorptionFile);
    end
end

for sp = smallPartConc
    for lp = largePartConc
        
        [phase, sig_s, waves] = calculateScattering(sp,lp);
        rtbWriteSpectrumFile(waves,sig_s,fullfile(resourceFolder,sprintf('scat_%.3f_%.3f.spd',sp,lp)));
        WritePhaseFile(waves,phase,fullfile(resourceFolder,sprintf('phase_%.3f_%.3f.spd',sp,lp)));
        
    end
end

rtbWriteSpectrumFile(wavesLed,daylight,fullfile(resourceFolder,'SunLight.spd'));

for ledID=1:nLeds
    rtbWriteSpectrumFile(wavesLed,ledSpectra(:,ledID),fullfile(resourceFolder,sprintf('LED%i.spd',ledID-1)));
end



%% Choose files to render


[scene, elements] = mexximpCleanImport(parentSceneFile,...
    'ignoreRootTransform',false,...
    'toReplace',{'png','jpg'},...
    'targetFormat','exr',...
    'flipUVs',true,...
    'exrToolsImage','hblasins/imagemagic-docker',...
    'workingFolder',resourceFolder);


%% Start rendering
condition = {'surface','water'};
illumination = {'ambient','noAmbient'};
nConditions = length(targetDistance)*length(depth)*length(chlConc)*length(cdomConc)*...
    length(smallPartConc)*length(largePartConc)*nLeds*2;

names = {'imageName','mode','illumination','ledID','pixelSamples','volumeStep','filmDist','filmDiag','cameraDistance','depth',...
    'chlConc','cdomConc','smallPartConc','largePartConc','lookAt','viewAxis'};

values = cell(nConditions,numel(names));
cntr = 1;
for td = targetDistance
    for de = depth
        for cC = chlConc
            for cdom = cdomConc
                for sp = smallPartConc
                    for lp = largePartConc
                        for ledID=0:(nLeds-1)
                            for c=1:length(condition)
                                for a=1:length(illumination)
                                    
                                    c1 = strcmp(condition{c},'water') && strcmp(illumination{a},'ambient') && ledID>=0 && ledID<=7;
                                    c2 = strcmp(condition{c},'water') && strcmp(illumination{a},'noAmbient') && ledID>=1 && ledID<=7;
                                    c3 = strcmp(condition{c},'surface') && strcmp(illumination{a},'ambient') && ledID==0;
                                    c4 = strcmp(condition{c},'surface') && strcmp(illumination{a},'noAmbient') && ledID==8;
                                    
                                    if c1==0 && c2==0 && c3==0 && c4==0
                                        continue;
                                    end
                                    
                                    % Generate condition entries
                                    values(cntr,1) = {sprintf('%s_%s_%s_LED%i_%i_%i_%.3f_%.3f_%.3f_%.3f',model,condition{c},illumination{a},ledID,td,de,cC,cdom,sp,lp)};
                                    values(cntr,2) = condition(c);
                                    values(cntr,3) = illumination(a);
                                    values(cntr,4) = num2cell(ledID,1);
                                    values(cntr,5) = num2cell(nSamples,1);
                                    values(cntr,6) = num2cell(50,1);
                                    values(cntr,7) = num2cell(((filmDiag*td)/objectSize),1);
                                    values(cntr,8) = num2cell(filmDiag,1);
                                    values(cntr,9) = num2cell(td,1);
                                    values(cntr,10) = num2cell(de,1);
                                    values(cntr,11) = num2cell(cC,1);
                                    values(cntr,12) = num2cell(cdom,1);
                                    values(cntr,13) = num2cell(sp,1);
                                    values(cntr,14) = num2cell(lp,1);
                                    values(cntr,15) = {mat2str(lookAt)};
                                    values(cntr,16) = {mat2str(viewAxis)};
                                    
                                    cntr = cntr+1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

rtbWriteConditionsFile(conditionsFile,names,values);

nativeSceneFiles = rtbMakeSceneFiles(scene, 'hints', hints, ...
    'conditionsFile',conditionsFile);
%%
radianceDataFiles = rtbBatchRender(nativeSceneFiles, 'hints', hints);

for i=1:length(radianceDataFiles)
    radianceData = load(radianceDataFiles{i});
    
    % Create an oi
    oiParams.lensType = 'pinhole';
    oiParams.filmDistance = values{i,4};
    oiParams.filmDiag = 20;
    
    
    oi = BuildOI(radianceData.multispectralImage, [], oiParams);
    oi = oiSet(oi,'name',values{i,1});
    
    ieAddObject(oi);
    oiWindow;
end

