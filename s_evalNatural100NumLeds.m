% This script reads in the simulated patch data from the Natural 100 test
% chart and performs the spectral estimation using different number of
% LEDs.
%
% Copytight, Henryk Blasinski 2017.

close all;
clear all;
clc;

%%

ledSets = {[1,4,7],...
           [1,3,4,7],...
           [1,3,4,5,7],...
           [1,2,3,4,5,7],...
           [1,2,3,4,5,6,7]};
nLedSets = length(ledSets);

alphaVec = logspace(-3,1,10);
nAlpha = length(alphaVec);

fName = fullfile(cmfRootPath,'Results','simNatural100V2.mat');
load(fName);

%% Estimate the illuminant
%  from measurements

channelEstLin =  zeros(nX,nY,10*11,3,nCameras,nLedSets,nAlpha);
wghtEst = zeros(nX,nY,nChannels,nCameras,nLedSets,nAlpha);

for a=1:nAlpha
    for xx=1:nX
        for yy=1:nY
            
            [ spd ] = xy2Spectrum( xVec(xx), yVec(yy), wave );
            if sum(isnan(spd)) > 0
                continue,
            end
            
            % figure;
            for cam=1:nCameras
                
                currentMeas = measurement{xx,yy,cam};
                
                for s=1:nLedSets
                    
                    nLEDs = length(ledSets{s});
                    currentFlash = flashNorm(:,ledSets{s});
                    
                    ambient = currentMeas.patch.ambient;
                    data = currentMeas.patch.led(:,ledSets{s},:);
                    
                    ambient = ambient/max(ambient(:));
                    data = data/max(data(:));
                    
                    
                    [ ambientEst, ambientWghts, ambientPredictions ] = globalAmbientEst( ambient, data, currentFlash, 'alpha', alphaVec(a) );
                    
                    %{
                    subplot(2,nLedSets,s);
                    hold on; grid on; box on;
                    plot(wave,ambientEst);
                    xlabel('Wavelength, nm');
                    title(sprintf('%i LEDs',nLEDs));
                                        
                    
                    subplot(2,nLedSets,s+nLedSets);
                    hold on; grid on; box on;
                    plot(squeeze(ambient)',ambientPredictions','.');
                    xlabel('Simulated');
                    ylabel('Approximation');
                    %}
                    
                    channelEstLin(xx,yy,:,:,cam,s,a) = ambientPredictions';
                    wghtEst(xx,yy,1:nLEDs,cam,s,a) = ambientWghts;
                    
                end
            end
            
        end
    end
    
    % Save data
    fName = fullfile(cmfRootPath,'Results','uniformApproxV3.mat');
    save(fName);
    
end



