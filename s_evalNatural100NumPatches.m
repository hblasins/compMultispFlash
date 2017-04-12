close all;
clear all;
clc;

ledSets = {[1,2,3,4,5,6,7]};
       
nLedSets = length(ledSets);

xChromaticities = [3];
yChromaticities = [3];

cameraIDs = [3];

nDraws = 10;

numPatches = [1:2:10, 20, 30];
nPatchSets = length(numPatches);

alphaVec = [0.06]; % logspace(-3,1,10);
nAlpha = length(alphaVec);

fName = fullfile(compMultispFlashRootPath,'Results','simNatural100.mat');
load(fName);

%% Estimate the illuminant
%  From measurements

lambda = 0;
R = [eye(nWaves-1) zeros(nWaves-1,1)] - [zeros(nWaves-1,1) eye(nWaves-1)];

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
            
            % Chromaticity (we only analyze a subset)
            for xx=xChromaticities 
                for yy=yChromaticities
                    
                    [ spd ] = xy2Spectrum( xVec(xx), yVec(yy), wave );
                    if sum(isnan(spd)) > 0
                        continue;
                    end
                    
                    measurementSet = squeeze(channelRawLin(xx,yy,:,:,:,:,:));
                    sz = size(measurementSet);
                    
                    tmp = reshape(measurementSet,[10*11, sz(3:end)]);
                    channelRawLinSubset = tmp(patchIndices,:,:,:);
                    
                    
                    % Regularization
                    for a=1:nAlpha
                        % Camera model
                        for cam=cameraIDs
                            
                            
                            ambientAppearance = channelRawLinSubset(:,:,1,cam);
                            ledAppearance = channelRawLinSubset(:,:,ledIndices+1,cam);
                            
                            cvx_begin
                                variables ambientApproxWghts(nLEDs,1)
                            
                                approx = 0;
                                for i=1:nLEDs
                                    approx = approx + ledAppearance(:,:,i)*ambientApproxWghts(i);
                                end
                            
                                minimize norm(ambientAppearance(:) - approx(:)) + alphaVec(a) * norm(R*flashNorm(:,ledIndices)*ambientApproxWghts,2)
                                subject to
                                    flashNorm(:,ledIndices)*ambientApproxWghts >= 0
                            cvx_end
                            
                            channelEstLin{xx,yy,cam,s,a,d,p} = approx;
                            wghtEst{xx,yy,cam,s,a,d,p} = ambientApproxWghts;
                            
                        end
                    end
                    
                end
            end
        end
        
    end
end


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

