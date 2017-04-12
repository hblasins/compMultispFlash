close all;
clear all;
clc;

% Prepare synthetic data
nChannels = 5;
nFilters = 3;
imHeight = 5;
imWidth = 7;
nWaves = 30;

Led = rand(nWaves,nChannels);

w = rand(imHeight,imWidth,nChannels);
L = rand(imHeight,imWidth,nFilters,nChannels);
m = zeros(imHeight, imWidth, nFilters);
for i=1:nChannels
    m = m + L(:,:,:,i).*repmat(w(:,:,i),[1, 1, nFilters]);
end

%% Solver tuning parameters
alpha = 1;
beta = 1;

%% CVX

estCVX = ambientEstCVX( m, L, Led, alpha, beta );

%% Implement ADMM solver

[ estADMM, pred, hist ] = ambientEstADMM( m, L, Led, alpha, beta, 'verbose', true, 'maxIter', 1000, 'tol', 0);


%% Compare results

figure;
subplot(1,2,1);
imagesc(estCVX(:,:,1)); axis image;
title('CVX');
subplot(1,2,2);
imagesc(estADMM(:,:,1)); axis image;
title('ADMM');


figure;
hold on; grid on; box on;
plot(estCVX(:),estADMM(:),'x','lineWidth',2,'markerSize',3);
xlabel('CVX');
ylabel('ADMM');
title('ADMM solver accuracy: weights');

figure;
hold on; grid on; box on;
plot([hist.prRes hist.dualRes]);
set(gca,'yscale','log');
xlabel('Iteration');
ylabel('Residal');
title('ADMM solver convergence');
legend('Primal','Dual');




