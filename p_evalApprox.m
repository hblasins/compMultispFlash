% Plot the data generated with the s_evalApprox.m script

close all;
clear all;
clc;


fName = fullfile(cmfRootPath,'Results','evalApprox.mat');
load(fName);

%% Ambient illuminant approximation
%  Compare the ambient illuminant chromaticity and the chromaticity of the
%  estimated ambient for different ambient chromaticities (xy) as well as
%  across cameras with different responsivity functions.

chromaticityPlot;
hold on; grid on; box on;
for xx=1:nXVals
    for yy=1:nYVals
        for cc=1:nCameras
            estXY = estXYZ(xx,yy,cc,1:2)/sum(estXYZ(xx,yy,cc,:));
            plot(estXY(1),estXY(2),'x');
        end
        
        
        cc = 1;
        ambXY = ambientXYZ(xx,yy,cc,1:2)/sum(ambientXYZ(xx,yy,cc,:));
        
        plot(ambXY(1),ambXY(2),'ro','lineWidth',2,'markerSize',5);
    end
    
end

% Add LED chromaticities
for i=1:nLEDs
    xyz = ieXYZFromPhotons(flashSpd(:,i),wave);
    plot(xyz(1)/sum(xyz),xyz(2)/sum(xyz),'k+','lineWidth',2,'markerSize',10);
end

%% RGB approximation errors
%  Give an idea of how well can the LED only RGB images approximate the a
%  ambient illuminant RGB image.

figure;
hold on; grid on; box on;
for cc=1:nCameras
        
    for xx=1:nXVals
        for yy=1:nYVals
            
                measured = squeeze(measPixelVals(xx,yy,cc,:,1,:));
                approximated = squeeze(approxPixelVals(xx,yy,cc,:,1,:));
                plot(measured(:),approximated(:),'k.');
        end
    end
    
end

xlabel('Measured RGB');
ylabel('LED approximated RGB');

%% Complement to desired illuminant error
%  Given a particular ambient illuminant chromaticity, and a desited
%  illuminant estimate how well can the effective illuminant (ambient+flash)
%  chromaticity approximate the desired illuminant chromaticity

figure;
for xx=5 %1:nXVals
    for yy=5 %1:nYVals
        
        tmp = ambientXYZ(xx,yy,:,:);
        if sum(tmp(:)) == 0, continue; end;
        
        for dd=1:nDesiredIlluminants
            
            subplot(2,2,dd);
            chromaticityPlot([],'white',500,false);
            hold on; grid on; box on;
            
            % Plot the ambient illuminant
            for cc=1:nCameras
                estXY = estXYZ(xx,yy,cc,1:2)/sum(estXYZ(xx,yy,cc,:));
                plot(estXY(1),estXY(2),'*');
            end
            
            
            cc = 1;
            ambXY = ambientXYZ(xx,yy,cc,1:2)/sum(ambientXYZ(xx,yy,cc,:));
            
            plot(ambXY(1),ambXY(2),'ro','lineWidth',2,'markerSize',10);
            
            
            % Plot the complemented illuminant (Estimated weights are constrained to [0, 1])
            for cc=1:nCameras
                complXY = ambientComplementConstrXYZ(xx,yy,cc,dd,1:2)/sum(ambientComplementConstrXYZ(xx,yy,cc,dd,:));
                plot(complXY(1),complXY(2),'x');
            end
            
            % Plot the complemented illuminant (Estimated weights are constrained to [0, Inf])
            % This is the 'full' computational mode since the flash would
            % be delivering more energy than it can produce.
            for cc=1:nCameras
                complXY = ambientComplementUncXYZ(xx,yy,cc,dd,1:2)/sum(ambientComplementUncXYZ(xx,yy,cc,dd,:));
                plot(complXY(1),complXY(2),'d');
            end
            
            % Plot the desired illuminant.
            desXYZ = ieXYZFromPhotons(desiredIlluminants{dd},wave);
            desXY = desXYZ(1:2)/sum(desXYZ);
            
            plot(desXY(1),desXY(2),'m+','lineWidth',3,'markerSize',15);
            
            
            % Plot the individual LED spectra.
            for i=1:nLEDs
                xyz = ieXYZFromPhotons(flashSpd(:,i),wave);
                plot(xyz(1)/sum(xyz),xyz(2)/sum(xyz),'k+','lineWidth',3,'markerSize',15);
            end
            
            title(sprintf('Desired illuminant %i',dd));
            drawnow;
        end
    end
    
end



