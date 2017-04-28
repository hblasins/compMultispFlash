% Display plots showing the variability in illuminant estimate as a
% function of the number of different reflectance spectra (patches) used
% for estimation. This script reproduces Figure 5 from the manuscript.
%
% Copyright, Henryk Blasinski 2017.

close all;
clear all;
clc;

destPath = [];
% destPath = fullfile(cmfRootPath,'..','Figures');

fName = fullfile(cmfRootPath,'Results','numPatches.mat');
load(fName);


%% Display Plots

lw = 1.5;

figure;
hold on; grid on; box on;

set(gca,'TickLabelInterpreter','Latex');
set(gcf,'Units','Centimeters');
set(gca,'FontSize',6);
set(gcf,'PaperPosition',[1 1 5 2]);

selXchrom = [3];
selYchrom = [3];
selCameras = 3;

for s=1:nLedSets 

    data = estimationError(selXchrom,selYchrom,selCameras,s,:,:,:);
    data = min(data,5);
    sz = size(data);

    reshapedData = reshape(data,[sz(1)*sz(2)*sz(3)*sz(4)*sz(5)*sz(6), sz(7)]);

    avgErr = mean(reshapedData);
    stdErr = std(reshapedData);


    errorbar(numPatches,avgErr,stdErr,'LineWidth',lw);
    plot(numPatches,ones(nPatchSets,1)*referenceError(selXchrom,selYchrom,s),'--','LineWidth',lw);
    
end



ylim([0.1 0.5]);
xlim([0 max(numPatches)+1]);
ylabel('Illuminant estimate RMS error','Interpreter','latex');
xlabel('Number of surface reflectance spectra','Interpreter','latex');
legend({'Proposed','Error bound'},'Location','NorthEast','Interpreter','latex')


if isempty(destPath) == false
    fName = fullfile(destPath,'numPatches.eps');
    print('-depsc',fName);
end






