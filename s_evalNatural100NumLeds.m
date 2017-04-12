close all;
clear all;
clc;

ledSets = {[1,4,7],...
           [1,3,4,7],...
           [1,3,4,5,7],...
           [1,2,3,4,5,7],...
           [1,2,3,4,5,6,7]};
nLedSets = length(ledSets);

alphaVec = logspace(-3,1,10);
nAlpha = length(alphaVec);

fName = fullfile(compMultispFlashRootPath,'Results','simNatural100.mat');
load(fName);

%% Estimate the illuminant
%  From measurements

lambda = 0;
R = [eye(nWaves-1) zeros(nWaves-1,1)] - [zeros(nWaves-1,1) eye(nWaves-1)];

channelEstLin =  zeros(nX,nY,10,11,3,nCameras,nLedSets,nAlpha);
wghtEst = zeros(nX,nY,nChannels,nCameras,nLedSets,nAlpha);

for a=1:nAlpha
    for xx=1:nX
        for yy=1:nY
            
            [ spd ] = xy2Spectrum( xVec(xx), yVec(yy), wave );
            if sum(isnan(spd)) > 0
                continue,
            end
            
            for cam=1:nCameras
                for s=1:nLedSets
                    
                    nLEDs = length(ledSets{s});
                    
                    
                    ambientAppearance = channelRawLin(xx,yy,:,:,:,1,cam);
                    
                    cvx_begin
                        variables ambientApproxWghts(nLEDs,1)
                    
                        approx = 0;
                        for i=1:nLEDs
                            approx = approx + squeeze(channelRawLin(xx,yy,:,:,:,ledSets{s}(i)+1,cam))*ambientApproxWghts(i);
                        end
                    
                        minimize norm(ambientAppearance(:) - approx(:)) + alphaVec(a) * norm(R*flashNorm(:,ledSets{s})*ambientApproxWghts,2)
                        subject to
                            flashNorm(:,ledSets{s})*ambientApproxWghts >= 0
                    cvx_end
                    
                    channelEstLin(xx,yy,:,:,:,cam,s,a) = approx;
                    wghtEst(xx,yy,1:nLEDs,cam,s,a) = ambientApproxWghts;
                    
                end
            end
            
        end
    end
end

% Save data
fName = fullfile(compMultispFlashRootPath,'Results','uniformApprox.mat');
save(fName);

%% Display Plots

close all;

fName = fullfile(compMultispFlashRootPath,'Results','uniformApprox.mat');
load(fName);

set(groot,'defaultAxesColorOrder',[1 0 0; 0 1 0; 0 0 1]);



pixelRMSE = zeros(nX,nY,nCameras,nLedSets,nAlpha);
chromaticityRMSE = zeros(nX,nY,nCameras,nLedSets,nAlpha);
estimateRMSE = zeros(nX,nY,nCameras,nLedSets,nAlpha);
spectrumRMSE = zeros(nX,nY,nLedSets,nAlpha);

validChromaticity = zeros(nX,nY);

generateFigures = true;

% We pick a specific chromaticity xx,yy and a tuning parameter a.
for a=5 %:nAlpha
    for xx=5 %1:nX
        for yy=5 %1:nY
            
            spd = xy2Spectrum( xVec(xx), yVec(yy), wave );
            if sum(isnan(spd)) > 0
                continue,
            end
            validChromaticity(xx,yy) = 1;
            spdXYZ = ieXYZFromPhotons(spd,wave);
            spdxy = spdXYZ(1:2)/sum(spdXYZ);
            
            if generateFigures;
                resFig = figure;
                set(gcf,'Position',[539 515 1184 600]);
            end
            
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
                
                chromaticityRMSE(xx,yy,:,s,a) = mean((spdEstxy - repmat(spdxy',[1,nCameras])).^2,1);
                
                if generateFigures
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
                    
                end
                
                
                
                
                
                appear = squeeze(channelRawLin(xx,yy,:,:,:,1,:));
                appear = reshape(appear,[10*11, 3, nCameras]);
                approx = squeeze(channelEstLin(xx,yy,:,:,:,:,s,a));
                approx = reshape(approx,[11*10, 3, nCameras]);
                
                pixelRMSE(xx,yy,:,s,a) = sqrt(sum(sum((appear - approx).^2,1),2)/(10*11*3));
                
                appear = permute(appear,[1 3 2]);
                appear = reshape(appear,[10*11*nCameras, 3]);
                
                approx = permute(approx,[1 3 2]);
                approx = reshape(approx,[10*11*nCameras, 3]);
                
                
                err = sqrt(mean((appear(:) - approx(:)).^2));
                
                if generateFigures
                    set(gca,'ColorOrder',[1 0 0; 0 1 0; 0 0 1]);
                    subplot(3,length(ledSets),length(ledSets) + s);
                    hold on; grid on; box on;
                    plot(approx,appear,'.');
                    xlabel('Weighed sum');
                    ylabel('Simulation');
                    xlim([0 1]);
                    ylim([0 1]);
                    title(sprintf('RMSE %.4f',err));
                end
                
                
                
                avgSpdEst = mean(spdEst,2);
                minSpdEst = min(spdEst,[],2);
                maxSpdEst = max(spdEst,[],2);
                
                spdTrue = spd/max(spd(:));
                
                gf = spdTrue(:)\avgSpdEst(:);
                
                err = sqrt(mean((spdTrue(:)*gf - avgSpdEst(:)).^2));
                
                
                estimateRMSE(xx,yy,:,s,a) = mean(sqrt(mean((repmat(spdTrue(:)*gf,[1, nCameras]) - spdEst).^2)));
                
                
                if generateFigures
                    subplot(3,length(ledSets),2*length(ledSets) + s);
                    hold on; grid on; box on;
                    fl = fill([wave; flipud(wave)],[maxSpdEst; flipud(minSpdEst)],[1 0.8 0.8]);
                    pl = plot(wave,[spdTrue*gf, avgSpdEst],'LineWidth',2);
                    set(pl(1),'Color','g');
                    set(pl(2),'Color','r');
                    set(fl,'EdgeColor','none');
                    xlabel('Wavelength, nm');
                    xlim([min(wave), max(wave)]);
                end
                
                % Direct least squares fit
                cvx_begin
                variable lsWeights(nLEDs,1)
                minimize norm(spdTrue*gf - flashNorm(:,ledSets{s})*lsWeights)
                subject to
                flashNorm(:,ledSets{s})*lsWeights >= 0
                cvx_end
                
                errLs = sqrt(mean((spdTrue(:)*gf - flashNorm(:,ledSets{s})*lsWeights).^2));
                
                if generateFigures
                    plot(wave,flashNorm(:,ledSets{s})*lsWeights,'--b','LineWidth',2);
                    title(sprintf('RMSE %.2f (%.2f)',err,errLs));
                end
                
                spectrumRMSE(xx,yy,s) = sqrt(mean((spdTrue*gf - flashNorm(:,ledSets{s})*lsWeights).^2));
                
                
                
                % drawnow;
                % print('-dpng',fullfile('./','Images',sprintf('spectra-%i-%i-%i.png',xx,yy,s)));
                
            end
            
            if generateFigures
                drawnow;
            end
        end
    end
end


%% Analyze errors

% pixelRMSE = zeros(nX,nY,nCameras,nLedSets,nAlpha);
% chromaticityRMSE = zeros(nX,nY,nCameras,nLedSets,nAlpha);
% estimateRMSE = zeros(nX,nY,nCameras,nLedSets,nAlpha);
% spectrumRMSE = zeros(nX,nY,nLedSets,nAlpha);


for i=1:nLedSets
   
    % Pixel
    subPixelRMSE = pixelRMSE(:,:,:,i,:);
    subValid = repmat(validChromaticity,[1 1 nCameras 1 1]);
    subValid = logical(subValid(:));
    
    subPixelRMSE = reshape(subPixelRMSE,[nX*nY*nCameras, nAlpha]);
    subPixelRMSE = subPixelRMSE(subValid,:);
    
    [avgSubPixelRMSE, aValPixel] = min(mean(subPixelRMSE));
    
    % Chromaticity
    
    subChrRMSE = chromaticityRMSE(:,:,:,i,:);
    subValid = repmat(validChromaticity,[1 1 nCameras 1 1]);
    subValid = logical(subValid(:));
    
    subChrRMSE = reshape(subChrRMSE,[nX*nY*nCameras, nAlpha]);
    subChrRMSE = subChrRMSE(subValid,:);
    
    [avgSubChrRMSE, aValChr] = min(mean(subChrRMSE));
    
    % Chromaticity
    
    subEstRMSE = estimateRMSE(:,:,:,i,:);
    subValid = repmat(validChromaticity,[1 1 nCameras 1 1]);
    subValid = logical(subValid(:));
    
    subEstRMSE = reshape(subEstRMSE,[nX*nY*nCameras, nAlpha]);
    subEstRMSE = subEstRMSE(subValid,:);
    
    [avgSubEstRMSE, aValEst] = min(mean(subEstRMSE));
   
    
    % Best LS
    
    subSpectrRMSE = spectrumRMSE(:,:,i,1);
    subValid = repmat(validChromaticity,[1 1 1 1]);
    subValid = logical(subValid(:));
    
    subSpectrRMSE = reshape(subSpectrRMSE,[nX*nY, 1]);
    subSpectrRMSE = subSpectrRMSE(subValid,:);
    
    avgSubSpectrRMSE = mean(subSpectrRMSE);
    
    fprintf('%10i | %10f | %10f | %10f | %10f \n',i+2,avgSubPixelRMSE,avgSubChrRMSE,avgSubEstRMSE,avgSubSpectrRMSE);
    fprintf('%10i | %10f | %10f | %10f  \n',i+2,aValPixel,aValChr,aValEst);
    
end


