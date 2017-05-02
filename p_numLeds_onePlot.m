% This script is very similar to p_numLeds.m wit the exception that all
% plots are generated in one figure (for presentation purposes).
%
% Copyright, Henryk Blasinski 2017.

close all;
clear all;
clc;

% Provide non-empty destPath to save figure files.
destPath = [];
% destPath = fullfile(slRootPath,'Figures');

fName = fullfile(cmfRootPath,'Results','uniformApproxV3.mat');
load(fName);

set(groot,'defaultAxesColorOrder',[1 0 0; 0 1 0; 0 0 1]);

%% Version 1: (Figure. 5)
%  Fix the illuminant chromaticity, vary the number and positions of leds,
%  camera spectral responsivity. Show pixel approximation scatter plots,
%  illuminant spectra estimates and illuminant chromaticity estimates.

% Pick regularization strength
figure;
set(gcf,'Position',[55 183 1373 912]);

for a=1
    
    % Pick illuminant chromaticity index
    for xx=5
        for yy=5
            
            spd = xy2Spectrum( xVec(xx), yVec(yy), wave );
            if sum(isnan(spd)) > 0
                continue,
            end
            spdXYZ = ieXYZFromPhotons(spd,wave);
            spdxy = spdXYZ(1:2)/sum(spdXYZ);
            
            
            for s=1:nLedSets
                
                nLEDs = length(ledSets{s});
                
                spdEst = zeros(nWaves,nCameras);
                spdEstXYZ = zeros(3,nCameras);
                spdEstxy = zeros(2,nCameras);
                for c=1:nCameras
                    spdEst(:,c) = flashSpd(:,ledSets{s})*squeeze(wghtEst(xx,yy,1:nLEDs,c,s,a));
                    spdEst(:,c) = spdEst(:,c)/max(spdEst(:,c));
                    spdEstXYZ(:,c) = ieXYZFromPhotons(spdEst(:,c)*10^10,wave);
                    spdEstxy(:,c) = spdEstXYZ(1:2,c)/sum(spdEstXYZ(:,c));
                    
                end
                
                subplot(3,nLedSets,s);
                hold on; grid on; box on;
                chromaticityPlot([],'gray',256,false);
                plot(spdxy(1),spdxy(2),'rd','lineWidth',2,'markerSize',10);
                title(sprintf('%i LEDs',nLEDs));
                
                
                for i=1:nLEDs
                    xyz = ieXYZFromPhotons(flashSpd(:,ledSets{s}(i)),wave);
                    plot(xyz(1)/sum(xyz),xyz(2)/sum(xyz),'+','lineWidth',2,'markerSize',10,'color',[0.5 0.5 0.5]);
                end
                
                plot(spdEstxy(1,:),spdEstxy(2,:),'kx');
                xlabel('CIE x');
                ylabel('CIE y');
                
                appear = [measurement{xx,yy,:}];
                appear = [appear(:).patch];
                appear = [appear(:).ambient];
                appear = appear./repmat(max(max(appear,[],1),[],3),[3 1 10*11]);
                appear = permute(appear,[3,2,1]);
                appear = reshape(appear,[10*11*nCameras,3]);
                
                
                approx = squeeze(channelEstLin(xx,yy,:,:,:,s,a));
                approx = permute(approx,[1,3,2]);
                approx = reshape(approx,[11*10*nCameras, 3]);
                
                
                subplot(3,nLedSets,s+nLedSets);
                hold on; grid on; box on;
                plot(approx,appear,'.');
                
                xlim([0 1]);
                ylim([0 1]);
                xlabel('Approximated');
                ylabel('Measured');
                %}
                
                
                avgSpdEst = mean(spdEst,2);
                minSpdEst = min(spdEst,[],2);
                maxSpdEst = max(spdEst,[],2);
                
                spdTrue = spd/max(spd(:));
                
                gf = spdTrue(:)\avgSpdEst(:);
                
                err = sqrt(mean((spdTrue(:)*gf - avgSpdEst(:)).^2));
                
                
                
                subplot(3,nLedSets,s+2*nLedSets);
                hold on; grid on; box on;
                fl = fill([wave; flipud(wave)],[maxSpdEst; flipud(minSpdEst)],[1 0.8 0.8]);
                pl = plot(wave,[spdTrue*gf, avgSpdEst],'LineWidth',2);
                set(pl(1),'Color','g');
                set(pl(2),'Color','r');
                set(fl,'EdgeColor','none');
                xlabel('Wavelength, nm');
                xlim([min(wave), max(wave)]);
                
                
                % Direct least squares fit, we compute the approximation
                % error up to a multiplicative scale on the estimate (i.e.
                % first we find a gain parameter that minimizes the error,
                % and then compute the error).
                cvx_begin quiet
                    variable lsWeights(nLEDs,1)
                    minimize norm(spdTrue*gf - flashNorm(:,ledSets{s})*lsWeights)
                    subject to
                        flashNorm(:,ledSets{s})*lsWeights >= 0
                cvx_end
                
                errLs = sqrt(mean((spdTrue(:)*gf - flashNorm(:,ledSets{s})*lsWeights).^2));
                
                plot(wave,flashNorm(:,ledSets{s})*lsWeights,'--b','LineWidth',2);
                
                xlabel('Wavelength, nm');
                ylabel('Intensity, au');  
            end
        end
    end
end





%% Version 2 (Figure 6)
%  Fix the number of LEDs, vary the illuminant chromaticity and camera
%  spectral responsivity curves.

close all;

% Select the chromaticity indices of interest.
xxSel = [2 2 3 3 6];
yySel = [2 6 3 6 3];

a = 1; % Smoothness penalty
s = 5; % The fifth led set containing 7 LEDs

figure;
set(gcf,'Position',[55 183 1373 912]);

for i=1:length(xxSel)
    xx=xxSel(i);
    yy=yySel(i);
    
    spd = xy2Spectrum( xVec(xx), yVec(yy), wave );
    if sum(isnan(spd)) > 0
        continue,
    end
    spdXYZ = ieXYZFromPhotons(spd,wave);
    spdxy = spdXYZ(1:2)/sum(spdXYZ);
    
    
    nLEDs = length(ledSets{s});
    
    spdEst = zeros(nWaves,nCameras);
    spdEstXYZ = zeros(3,nCameras);
    spdEstxy = zeros(2,nCameras);
    for c=1:nCameras
        spdEst(:,c) = flashSpd(:,ledSets{s})*squeeze(wghtEst(xx,yy,1:nLEDs,c,s,a));
        spdEst(:,c) = spdEst(:,c)/max(spdEst(:,c));
        spdEstXYZ(:,c) = ieXYZFromPhotons(spdEst(:,c)*10^10,wave);
        spdEstxy(:,c) = spdEstXYZ(1:2,c)/sum(spdEstXYZ(:,c));
        
    end
    
    subplot(3,length(xxSel),i);
    hold on; grid on; box on;
    chromaticityPlot([],'gray',256,false);
    plot(spdxy(1),spdxy(2),'rd','lineWidth',2,'markerSize',10);
    
    
    for j=1:nLEDs
        xyz = ieXYZFromPhotons(flashSpd(:,ledSets{s}(j)),wave);
        plot(xyz(1)/sum(xyz),xyz(2)/sum(xyz),'+','lineWidth',2,'markerSize',10,'color',[0.5 0.5 0.5]);
    end
    
    plot(spdEstxy(1,:),spdEstxy(2,:),'kx');
    xlabel('CIE x');
    ylabel('CIE y');
    
    appear = [measurement{xx,yy,:}];
    appear = [appear(:).patch];
    appear = [appear(:).ambient];
    appear = appear./repmat(max(max(appear,[],1),[],3),[3 1 10*11]);
    appear = permute(appear,[3,2,1]);
    appear = reshape(appear,[10*11*nCameras,3]);
    
    
    approx = squeeze(channelEstLin(xx,yy,:,:,:,s,a));
    approx = permute(approx,[1,3,2]);
    approx = reshape(approx,[11*10*nCameras, 3]);
    
    
    
    subplot(3,length(xxSel),i + length(xxSel));
    hold on; grid on; box on;
    plot(approx,appear,'.');
    
    xlim([0 1]);
    ylim([0 1]);
    xlabel('Approximated');
    ylabel('Measured');
    set(gca,'YTick',0:0.2:1);
    
    avgSpdEst = mean(spdEst,2);
    minSpdEst = min(spdEst,[],2);
    maxSpdEst = max(spdEst,[],2);
    
    spdTrue = spd/max(spd(:));
    
    gf = spdTrue(:)\avgSpdEst(:);
    
    err = sqrt(mean((spdTrue(:)*gf - avgSpdEst(:)).^2));
    
    
    
    subplot(3,length(xxSel),i+2*length(xxSel));
    hold on; grid on; box on;
    fl = fill([wave; flipud(wave)],[maxSpdEst; flipud(minSpdEst)],[1 0.8 0.8]);
    pl = plot(wave,[spdTrue*gf, avgSpdEst],'LineWidth',2);
    set(pl(1),'Color','g');
    set(pl(2),'Color','r');
    set(fl,'EdgeColor','none');
    xlabel('Wavelength, nm');
    xlim([min(wave), max(wave)]);
    
    
    % Direct least squares fit, we compute the approximation
    % error up to a multiplicative scale on the estimate (i.e.
    % first we find a gain parameter that minimizes the error,
    % and then compute the error).
    cvx_begin
    variable lsWeights(nLEDs,1)
    minimize norm(spdTrue*gf - flashNorm(:,ledSets{s})*lsWeights)
    subject to
    flashNorm(:,ledSets{s})*lsWeights >= 0
    cvx_end
    
    errLs = sqrt(mean((spdTrue(:)*gf - flashNorm(:,ledSets{s})*lsWeights).^2));
    
    plot(wave,flashNorm(:,ledSets{s})*lsWeights,'--b','LineWidth',2);
    
    xlabel('Wavelength, nm');
    ylabel('Intensity, au');
    
end



