% Evaluate the illuminant estimation algorithm using pixel data from
% different number of distinct patches of uniform reflectance. To simplify
% the analysis we pick the illuminant spectrum in terms of CIE xy
% chromaticity, a camera model, and regularization.
%
% Copyright, Henryk Blaisnski 2017


close all;
clear all;
clc;

%%
fName = fullfile(cmfRootPath,'Results','simNatural100V2.mat');
load(fName);

% Pick number of LEDs
ledSets = {[1,2,3,4,5,6,7]};       
nLedSets = length(ledSets);

% Illuminant chromaticity
xChromaticities = [3];
yChromaticities = [3];

% Camera models
cameraIDs = [3];

% Number of bootstrap samples
nDraws = 10;

% Number of patches in an evaluation experiment
numPatches = [1:2:10, 20, 30];
nPatchSets = length(numPatches);

alphaVec = [0.06]; 
nAlpha = length(alphaVec);



%% Estimate the illuminant
%  From measurements

channelEstLin =  cell(nX,nY,nCameras,nLedSets,nAlpha,nDraws,nPatchSets);
wghtEst = cell(nX,nY,nCameras,nLedSets,nAlpha,nDraws,nPatchSets);
selectedPatches = cell(nDraws,nPatchSets);

% Random draw
for d=1:nDraws 
    % Patch set
    for p=1:nPatchSets 
        
        nPatches = numPatches(p);
        patchIndices = randi(110,[nPatches, 1]);
        
        selectedPatches{d,p} = patchIndices;
        
        %LED set
        for s=1:nLedSets 
            
            ledIndices = ledSets{s};
            nLEDs = length(ledIndices);
            
            currentFlash = flashNorm(:,ledIndices);
            
            % Chromaticity (we only analyze a subset)
            for xx=xChromaticities 
                for yy=yChromaticities
                    
                    [ spd ] = xy2Spectrum( xVec(xx), yVec(yy), wave );
                    if sum(isnan(spd)) > 0
                        continue;
                    end
                    
                    % Regularization
                    for a=1:nAlpha
                        % Camera model
                        for cam=cameraIDs
                            
                            
                             measurementSet = measurement{xx,yy,cam};
                             measurementSet = [measurementSet.patch];
                             ambient = measurementSet.ambient(:,:,patchIndices);
                             
                             measurementSet = [measurementSet.led];
                             measurementSet = measurementSet(:,ledIndices,patchIndices);

                            [ ambientEst, ambientWghts, ambientPredictions ] = globalAmbientEst( ambient, measurementSet, currentFlash, 'alpha', alphaVec(a) );

                            
                            channelEstLin{xx,yy,cam,s,a,d,p} = ambientPredictions;
                            wghtEst{xx,yy,cam,s,a,d,p} = ambientWghts;
                            
                        end
                    end
                    
                end
            end
        end
        
    end
end

% Save data
fName = fullfile(cmfRootPath,'Results','numPatches.mat');
save(fName);


%% Compute estimation errors


estimationError = zeros(nX,nY,cam,s,a,d,p);
referenceError = zeros(nX,nY,nLedSets,1);

% Random draw
for d=1:nDraws 
    % Patch set
    for p=1:nPatchSets 
        
        patchIndices = selectedPatches{d,p};
        
        %LED set
        for s=1:nLedSets 
            
            ledIndices = ledSets{s};
            nLEDs = length(ledIndices);
            
            % Chromaticity (we only analyze a subset)
            for xx=xChromaticities 
                for yy=yChromaticities
                    
                    [ spd ] = xy2Spectrum( xVec(xx), yVec(yy), wave );
                    if sum(isnan(spd)) > 0
                        continue;
                    end
                    
                    cvx_begin
                        variables vv(nLEDs,1)
                        minimize norm(spd - flashNorm(:,ledIndices)*vv)
                        subject to
                            flashNorm(:,ledIndices)*vv >= 0
                    cvx_end
                    
                    fNorm = flashNorm(:,ledIndices)*vv;
                    fNorm = fNorm/max(fNorm);
                    
                    cvx_begin
                        variable g(1,1)
                        minimize norm(fNorm - g*spd)
                    cvx_end
                    
                    referenceError(xx,yy,s) = rms(g*spd - fNorm);
                    
                    % Regularization
                    for a=1:nAlpha
                        % Camera model
                        for cam=cameraIDs
                            
                            estimatedSpectrum = flashNorm(:,ledIndices)*wghtEst{xx,yy,cam,s,a,d,p};
                            estimatedSpectrum = estimatedSpectrum/max(estimatedSpectrum);
                            
                            cvx_begin
                                variable g(1,1)
                                minimize norm(estimatedSpectrum - g*spd)
                            cvx_end
                            
                            err = rms(estimatedSpectrum - g*spd);
                            
                            estimationError(xx,yy,cam,s,a,d,p) = err;
                            
                        end
                    end
                    
                end
            end
        end
        
    end
end

% Save data
fName = fullfile(cmfRootPath,'Results','numPatches.mat');
save(fName);
