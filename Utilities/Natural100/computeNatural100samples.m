function [reflectances,wavelength] = computeNatural100samples(plotFlag)

if ieNotDefined('plotFlag'), plotFlag = false; end

rng(1);
Marks = {'+','o','x','s','d','v'};
Marks = Marks(randperm(length(Marks)));

reflectances = [];
rng(1);

% Compute principal components
load MunsellSamples_Vhrel
munsellData = data;
load DupontPaintChip_Vhrel
dupontData = data;
load Objects_Vhrel.mat
objectsData = data;
data = [munsellData,dupontData,objectsData];
dataMean = mean(data,2);
dataM = bsxfun(@minus,data,dataMean);
[U,S,~] = svd(data,'econ');
U = U(:,1:7);
S = S(1:7,1:7);
% Projection matrix
P = U';
% P = diag(1./diag(S))*U'; % Normalize variances

if plotFlag
    close all
    figure, set(gcf,'Position',[0 0 800 1500])
    
    plot(wavelength,bsxfun(@plus,U,-(0:6)*.8),'LineWidth',2)    
    set(gca,'YTick',[])
    axis tight, set(gca,'FontSize',13,'FontWeight','b')
    legend('Component 1','Component 2','Component 3','Component 4','Component 5','Component 6','Component 7','Location','northeastoutside')
        
    figure, set(gcf,'Position',[0 0 800 400])
    plot(wavelength,dataMean,'k--','LineWidth',2), axis tight, set(gca,'FontSize',13,'FontWeight','b'), ylim([0 Inf])
    legend('PCA Mean','Location','northeastoutside')
end

   
% Clothes reflectances
load Clothes_Vhrel
data = removeOutOfGamut(data,wavelength,{'D65','Tungsten','Fluorescent'});
dataM = bsxfun(@minus,data,dataMean);
N = 20;
r = zeros(N,length(wavelength));
c = P * dataM;
cPCA = zeros(N,7);
Z = linkage(c');
T = cluster(Z,'maxclust',N);
for i = 1:N
    k = find(T == i);
    k = k(randperm(length(k),1));
    r(i,:) = data(:,k);
    cPCA(i,:) = c(:,k);
end

if plotFlag
    plot(wavelength,r(1:7,:)','LineWidth',2), axis tight, set(gca,'FontSize',13,'FontWeight','b'), ylim([0 Inf])
    legend('Clothes 1','Clothes 2','Clothes 3','Clothes 4','Clothes 5','Clothes 6','Clothes 7','Location','northeastoutside')
    set(gcf,'Position',[0 0 700 400])
    scatter(c(1,:),c(2,:),50,'LineWidth',2), grid on, axis equal, set(gca,'FontSize',13,'FontWeight','b')
    xlabel('PCA Component 1'), ylabel('PCA Component 2')
    for i = 1:N
        k = find(T == i);
        scatter(c(1,k),c(2,k),50,Marks{mod(i-1,6)+1},'LineWidth',2), grid on, axis equal, set(gca,'FontSize',13,'FontWeight','b')
        xlabel('PCA Component 1'), ylabel('PCA Component 2'), hold all
    end
    hold off
    set(gcf,'Position',[0 0 800 400])
end
    

XYZ = ieXYZFromEnergy([r;ones(1,length(wavelength))],wavelength);
LAB = xyz2lab(XYZ(1:end-1,:),XYZ(end,:));
[~,I] = sortrows(LAB,[2 3 1]);
r = r(I,:);
reflectances = [reflectances;r];% size(reflectances)

% Food reflectances
load Food_Vhrel
data = removeOutOfGamut(data,wavelength,{'D65','Tungsten','Fluorescent'});
dataM = bsxfun(@minus,data,dataMean);
N = 20;
r = zeros(N,length(wavelength));
c = P * dataM;
cPCA = zeros(N,7);
Z = linkage(c');
T = cluster(Z,'maxclust',N);
for i = 1:N
    k = find(T == i);
    k = k(randperm(length(k),1));
    r(i,:) = data(:,k);
    cPCA(i,:) = c(:,k);
end

if plotFlag
    plot(wavelength,r(1:7,:)','LineWidth',2), axis tight, set(gca,'FontSize',13,'FontWeight','b'), ylim([0 Inf])
    legend('Food 1','Food 2','Food 3','Food 4','Food 5','Food 6','Food 7','Location','northeastoutside')
    set(gcf,'Position',[0 0 700 400])
    scatter(c(1,:),c(2,:),50,'LineWidth',2), grid on, axis equal, set(gca,'FontSize',13,'FontWeight','b')
    xlabel('PCA Component 1'), ylabel('PCA Component 2')
    for i = 1:N
        k = find(T == i);
        scatter(c(1,k),c(2,k),50,Marks{mod(i-1,6)+1},'LineWidth',2), grid on, axis equal, set(gca,'FontSize',13,'FontWeight','b')
        xlabel('PCA Component 1'), ylabel('PCA Component 2'), hold all
    end
    hold off
    set(gcf,'Position',[0 0 800 400])
end

XYZ = ieXYZFromPhotons([r;ones(1,length(wavelength))],wavelength);
LAB = xyz2lab(XYZ(1:end-1,:),XYZ(end,:));
[~,I] = sortrows(LAB,[2 3 1]);
r = r(I,:);
reflectances = [reflectances;r];% size(reflectances)

% Nature data
load Food_Vhrel
foodData = data;
load reflectances/Nature_Vhrel
data = setdiff(data',foodData','rows')';
data = removeOutOfGamut(data,wavelength,{'D65','Tungsten','Fluorescent'});
dataM = bsxfun(@minus,data,dataMean);
N = 20;
r = zeros(N,length(wavelength));
c = P * dataM;
cPCA = zeros(N,7);
Z = linkage(c');
T = cluster(Z,'maxclust',N);
for i = 1:N
    k = find(T == i);
    k = k(randperm(length(k),1));
    r(i,:) = data(:,k);
    cPCA(i,:) = c(:,k);
end

if plotFlag
    plot(wavelength,r(1:7,:)','LineWidth',2), axis tight, set(gca,'FontSize',13,'FontWeight','b'), ylim([0 Inf])
    legend('Nature 1','Nature 2','Nature 3','Nature 4','Nature 5','Nature 6','Nature 7','Location','northeastoutside')
    set(gcf,'Position',[0 0 700 400])
    scatter(c(1,:),c(2,:),50,'LineWidth',2), grid on, axis equal, set(gca,'FontSize',13,'FontWeight','b')
    xlabel('PCA Component 1'), ylabel('PCA Component 2')
    for i = 1:N
        k = find(T == i);
        scatter(c(1,k),c(2,k),50,Marks{mod(i-1,6)+1},'LineWidth',2), grid on, axis equal, set(gca,'FontSize',13,'FontWeight','b')
        xlabel('PCA Component 1'), ylabel('PCA Component 2'), hold all
    end
    hold off
    set(gcf,'Position',[0 0 800 400])
end

XYZ = ieXYZFromPhotons([r;ones(1,length(wavelength))],wavelength);
LAB = xyz2lab(XYZ(1:end-1,:),XYZ(end,:));
[~,I] = sortrows(LAB,[2 3 1]);
r = r(I,:);
reflectances = [reflectances;r];% size(reflectances)

% Hair samples
load Hair_Vhrel
data = removeOutOfGamut(data,wavelength,{'D65','Tungsten','Fluorescent'});
dataM = bsxfun(@minus,data,dataMean);
N = 5;
r = zeros(N,length(wavelength));
c = P * dataM;
cPCA = zeros(N,7);
Z = linkage(c');
T = cluster(Z,'maxclust',N);
for i = 1:N
    k = find(T == i);
    k = k(randperm(length(k),1));
    r(i,:) = data(:,k);
    cPCA(i,:) = c(:,k);
end

if plotFlag
    plot(wavelength,r(1:5,:)','LineWidth',2), axis tight, set(gca,'FontSize',13,'FontWeight','b'), ylim([0 Inf])
    legend('Hair 1','Hair 2','Hair 3','Hair 4','Hair 5','Location','northeastoutside')
    set(gcf,'Position',[0 0 700 400])
    scatter(c(1,:),c(2,:),50,'LineWidth',2), grid on, axis equal, set(gca,'FontSize',13,'FontWeight','b')
    xlabel('PCA Component 1'), ylabel('PCA Component 2')
    for i = 1:N
        k = find(T == i);
        scatter(c(1,k),c(2,k),50,Marks{mod(i-1,6)+1},'LineWidth',2), grid on, axis equal, set(gca,'FontSize',13,'FontWeight','b')
        xlabel('PCA Component 1'), ylabel('PCA Component 2'), hold all
    end
    hold off
    set(gcf,'Position',[0 0 800 400])
end

XYZ = ieXYZFromPhotons([r;ones(1,length(wavelength))],wavelength);
LAB = xyz2lab(XYZ(1:end-1,:),XYZ(end,:));
[~,I] = sortrows(LAB,[2 3 1]);
r = r(I,:);
reflectances = [reflectances;r];% size(reflectances)

% Skin samples
load reflectances/Skin_Vhrel
data = removeOutOfGamut(data,wavelength,{'D65','Tungsten','Fluorescent'});
dataM = bsxfun(@minus,data,dataMean);
N = 15;
r = zeros(N,length(wavelength));
c = P * dataM;
cPCA = zeros(N,7);
Z = linkage(c');
T = cluster(Z,'maxclust',N);
for i = 1:N
    k = find(T == i);
    k = k(randperm(length(k),1));
    r(i,:) = data(:,k);
    cPCA(i,:) = c(:,k);
end

if plotFlag
    plot(wavelength,r(1:7,:)','LineWidth',2), axis tight, set(gca,'FontSize',13,'FontWeight','b'), ylim([0 Inf])
    legend('Skin 1','Skin 2','Skin 3','Skin 4','Skin 5','Skin 6','Skin 7','Location','northeastoutside')
    set(gcf,'Position',[0 0 700 400])
    scatter(c(1,:),c(2,:),50,'LineWidth',2), grid on, axis equal, set(gca,'FontSize',13,'FontWeight','b')
    xlabel('PCA Component 1'), ylabel('PCA Component 2')
    for i = 1:N
        k = find(T == i);
        scatter(c(1,k),c(2,k),50,Marks{mod(i-1,6)+1},'LineWidth',2), grid on, axis equal, set(gca,'FontSize',13,'FontWeight','b')
        xlabel('PCA Component 1'), ylabel('PCA Component 2'), hold all
    end
    hold off
    set(gcf,'Position',[0 0 800 400])
end

XYZ = ieXYZFromPhotons([r;ones(1,length(wavelength))],wavelength);
LAB = xyz2lab(XYZ(1:end-1,:),XYZ(end,:));
[~,I] = sortrows(LAB,[2 3 1]);
r = r(I,:);
reflectances = [reflectances;r];% size(reflectances)

% Dupont Paint samples
load DupontPaintChip_Vhrel
data = removeOutOfGamut(data,wavelength,{'D65','Tungsten','Fluorescent'});
dataM = bsxfun(@minus,data,dataMean);
N = 20;
r = zeros(N,length(wavelength));
c = P * dataM;
cPCA = zeros(N,7);
Z = linkage(c');
T = cluster(Z,'maxclust',N);
for i = 1:N
    k = find(T == i);
    k = k(randperm(length(k),1));
    r(i,:) = data(:,k);
    cPCA(i,:) = c(:,k);
end

if plotFlag
    plot(wavelength,r(1:7,:)','LineWidth',2), axis tight, set(gca,'FontSize',13,'FontWeight','b'), ylim([0 Inf])
    legend('Dupont 1','Dupont 2','Dupont 3','Dupont 4','Dupont 5','Dupont 6','Dupont 7','Location','northeastoutside')
    set(gcf,'Position',[0 0 700 400])
    scatter(c(1,:),c(2,:),50,'LineWidth',2), grid on, axis equal, set(gca,'FontSize',13,'FontWeight','b')
    xlabel('PCA Component 1'), ylabel('PCA Component 2')
    for i = 1:N
        k = find(T == i);
        scatter(c(1,k),c(2,k),50,Marks{mod(i-1,6)+1},'LineWidth',2), grid on, axis equal, set(gca,'FontSize',13,'FontWeight','b')
        xlabel('PCA Component 1'), ylabel('PCA Component 2'), hold all
    end
    hold off
    set(gcf,'Position',[0 0 800 400])
end

XYZ = ieXYZFromPhotons([r;ones(1,length(wavelength))],wavelength);
LAB = xyz2lab(XYZ(1:end-1,:),XYZ(end,:));
[~,I] = sortrows(LAB,[2 3 1]);
r = r(I,:);
reflectances = [reflectances;r];% size(reflectances)

