close all;
clear all;
clc;

%% Generate synthetic data

nWaves = 30;
nLeds = 1;
nChannels = 10;
imHeight = 13;
imWidth = 7;

camera = rand(nWaves,nChannels);
leds = rand(nWaves,nLeds);
weights = rand(imHeight,imWidth,nLeds);
ill0 = rand(nWaves,1);

%% Solver tuning parameters
alpha = 0;
beta = 0;
fm = true;

%% CVX solver

[estCVX, slackCVX, pred] = complementEstCVXuniformSlack( ill0, weights, leds, camera, alpha, beta, 'flashMode',fm);


%% ADMM

[ estADMMunif, slackADMMunif, histUnif ] = complementEstADMMuniformSlack( ill0, weights, leds, camera, alpha, beta,...
                                            'verbose', true, 'maxIter', 1000, 'tol', 0, 'flashMode',fm, 'rescaleRho',false);
% [ estADMM, slackADMM, hist ] = complementEstADMM( ill0, weights, leds, camera, alpha, beta, 'verbose', true, 'maxIter', 200, 'tol', 0 );

%% Compare results

figure;
hold on; grid on; box on;
plot(estCVX(:),estADMMunif(:),'x','lineWidth',2,'markerSize',3);
xlabel('CVX');
ylabel('ADMM');
title('ADMM solver accuracy: weights');

figure;
hold on; grid on; box on;
plot(slackCVX(:), slackADMMunif(:),'x','lineWidth',2,'markerSize',3);
xlabel('CVX');
ylabel('ADMM');
title('ADMM solver accuracy: slack var.');

figure;
hold on; grid on; box on;
plot([histUnif.prRes histUnif.dualRes]);
set(gca,'yscale','log');
xlabel('Iteration');
ylabel('Residual');
title('ADMM solver convergence');
legend('Primal','Dual');


