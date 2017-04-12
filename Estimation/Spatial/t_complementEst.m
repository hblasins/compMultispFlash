close all;
clear all;
clc;

%% Generate synthetic data

nWaves = 30;
nLeds = 5;
nChannels = 3;
imHeight = 13;
imWidth = 17;

camera = rand(nWaves,nChannels);
leds = rand(nWaves,nLeds);
weights = rand(imHeight,imWidth,nLeds);
ill0 = rand(nWaves,1);

%% Solver tuning parameters
alpha = 0.1;
beta = 0.1;

%% CVX solver

[estCVX, slackCVX] = complementEstCVX( ill0, weights, leds, camera, alpha, beta);

%% ADMM

[ estADMM, slackADMM, hist ] = complementEstADMM( ill0, weights, leds, camera, alpha, beta, 'verbose', true, 'maxIter', 1000, 'tol', 0 );

%% Compare results

figure;
hold on; grid on; box on;
plot(estCVX(:),estADMM(:),'x','lineWidth',2,'markerSize',3);
xlabel('CVX');
ylabel('ADMM');
title('ADMM solver accuracy: weights');

figure;
hold on; grid on; box on;
plot(slackCVX(:), slackADMM(:),'x','lineWidth',2,'markerSize',3);
xlabel('CVX');
ylabel('ADMM');
title('ADMM solver accuracy: slack var.');

figure;
hold on; grid on; box on;
plot([hist.prRes hist.dualRes]);
set(gca,'yscale','log');
xlabel('Iteration');
ylabel('Residal');
title('ADMM solver convergence');
legend('Primal','Dual');


