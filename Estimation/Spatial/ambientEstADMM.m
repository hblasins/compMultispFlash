function [ X, pred, hist ] = ambientEstADMM( ambient, ledMeasurements, leds, alpha, beta, varargin )

% [ X, hist ] = spectralEstADMM( meas, cameraMat, basisFcns, alpha, beta, gamma, ... )
%
% This function computes the reflectance basis function weights that best
% explain the data in meas using the ADMM approach. See the following paper
% for details
%
% H. Blasinski, J. Farrell and B. Wandell; 'Iterative algorithm for
% spectral estimation with spatial smoothing,' ICIP2015, Quebec City.
%
%
% Parameters:
%  meas - a h x w x c 3D matrix, where the first two dimensions represent
%      image height and width. The last dimension is the number of imaging
%      channels
%
%  cameraMat - a c x nWaves matrix where each row represents the spectral
%      responsivity of a given camera channel. It is assumed that the
%      spectrum is discretized to nWaves bins.
%
%  basisFcns - a nWaves x nBasis matrix where each column is a reflectance
%      basis function discretized to nWaves spectral bins.
%
%  alpha - a scalar controling the spectral smoothess of the solution.
%  
%  beta - a scalar controlling the spatial smoothness of the solution.
%
%  gamma - a scalar controlling whether the solution is bounded to the
%     [0,1] interval. If gamma is equal to zero this constraint is NOT
%     enforced. If gamma is different from zero then the constraint IS
%     enfoced and the value of the parameter does not matter.
%
%
% Optional Parameters (parameter - value pairs):
%
%  'maxIter' - a scalar defining the maximal number of ADMM iterations
%     (default = 100).
%
%  'tol' - a scalar defining the desired tolerance, which once reached
%     terminated the main ADMM loop (default = 1e-3).
%
%  'rhoInit' - the initial value of the proximal term weight rho (see paper
%     for details, default = 0.1);
%
%  'rescaleRho' - a boolean value enableing or disabling dynamic rho update
%     heuristic (default = true).
%
%  'mu','tauIncr','tauDecr' - scalars controling the dynamic rho update
%     heuristic, see Boyd 2011 for details (default, mu = 10, tauIncr =
%     tauDect = 5);
%
%  'verbose' - a boolean value turning on/off ADMM convergence summary at
%     every iteration (default = true).
%
%  'reference' - a reference solution computed using other means. If this
%     argument is present, then the hist output will contain a convergence
%     to reference plot. This input is useful for algorithmic analysis
%     purposes (default = []).
%
%  'tolConv' - the tolerance for the convergence of the solution to the
%     provided reference (default = 0).
%
%  'gradient' - selects between isotropic and anisotropic graidents,
%     allowable values: 'isotropic','anisotropic' (default = 'anisotropic').
%
%
% Returns:
%  X - a h x w x nBasis matrix of estimated spectral reflectance weights.
%  
%  hist - a structure containing the state of the ADMM solver at every
%     iteration. History contains for example primal and dual residuals,
%     number of conjugate gradient iterations, values of parameter rho etc.
%
%
% Copyright, Henryk Blasinski 2015.

p = inputParser;
p.addParameter('maxIter',100,@isscalar);
p.addParameter('tol',1e-3,@isscalar);
p.addParameter('mu',10,@isscalar);
p.addParameter('tauIncr',5,@isscalar);
p.addParameter('tauDecr',5,@isscalar);
p.addParameter('rhoInit',0.1,@isscalar);
p.addParameter('rescaleRho',false,@islogical);
p.addParameter('reference',[]);
p.addParameter('tolConv',0);
p.addParameter('verbose',false);
p.addParameter('gradient','anisotropic');
p.parse(varargin{:});
inputs = p.Results;

nChannels = size(leds,2);
nWaves = size(leds,1);
h = size(ambient,1);
w = size(ambient,2);
nFilters = size(ambient,3);

ambient = reshape(ambient,h*w,nFilters)';

ledMeasurements = reshape(ledMeasurements,[h*w, nFilters, nChannels]);
ledMeasurements = permute(ledMeasurements,[2, 3, 1]);

ledMeasTledMeas = zeros(nChannels,nChannels,h*w);
for i=1:h*w
    ledMeasTledMeas(:,:,i) = ledMeasurements(:,:,i)'*ledMeasurements(:,:,i);
end

Z1 = zeros(nChannels,(h-1)*w + (w-1)*h);
Z1minus = zeros(size(Z1));
U1 = zeros(nChannels,(h-1)*w + (w-1)*h);

Z2 = zeros(nWaves,h*w);
Z2minus = zeros(size(Z2));
U2 = zeros(nWaves,h*w);

hist.prRes = zeros(inputs.maxIter,1);
hist.dualRes = zeros(inputs.maxIter,1);
hist.rho = inputs.rhoInit*ones(inputs.maxIter+1,1);
hist.pcg.iter = zeros(inputs.maxIter,1);
hist.pcg.acc = zeros(inputs.maxIter,1);


if isempty(inputs.reference) == false
    hist.conv = zeros(inputs.maxIter,1);
end

X = [];

for i=1:inputs.maxIter
    t1 = tic;
      
    
    % X update using cg
    t2 = tic;
    Atb = applyAt(ambient,Z1-U1,Z2-U2,h,w,ledMeasurements,leds, beta, hist.rho(i)/2);
    AtAhndl = @(x) applyAtA(x, h, w, ledMeasTledMeas, leds, alpha, beta, hist.rho(i)/2);
    [X, ~, hist.pcg.acc(i), hist.pcg.iter(i)] = pcg(AtAhndl,Atb,1e-6,10000,[],[],X(:));
    t3 = toc(t2);
    
    X = reshape(X,nChannels,h*w);
    wghtImg = reshape(X',[h,w,nChannels]);

    
    
    
    % Anisotropic TV: Soft thresholding on Z1
    if beta > 0
        
        [grX, grY] = computeGradient(wghtImg);
        dX = [reshape(grX,h*(w-1),nChannels)' reshape(grY,(h-1)*w,nChannels)'];
        tmp = U1 + dX;
        
        switch inputs.gradient
            case 'isotropic'
                
                Vx = [zeros(h,1,nChannels) reshape(tmp(:,1:h*(w-1))',[h,w-1,nChannels])];
                Vy = [zeros(1,w,nChannels); reshape(tmp(:,h*(w-1)+1:end)',[h-1,w,nChannels])];
    
                normV = sqrt(Vx.^2 + Vy.^2);
                nu = beta./(hist.rho(i)*normV - beta);
                cond = normV > beta/hist.rho(i);
    
                Z1x = zeros(h,w,nChannels); 
                Z1y = zeros(h,w,nChannels);
                Z1x(cond) = 1./(1+nu(cond)) .* Vx(cond);
                Z1y(cond) = 1./(1+nu(cond)) .* Vy(cond);
                
                Z1x = Z1x(:,2:end,:);
                Z1y = Z1y(2:end,:,:);
                
                Z1 = [reshape(Z1x,h*(w-1),nChannels)' reshape(Z1y,(h-1)*w,nChannels)'];
                
            otherwise % 'anisotropic'

                Z1 = sign(tmp).*max(abs(tmp) - beta/hist.rho(i),0);
        end
            
    else
        Z1 = zeros(size(Z1));
    end
    
    % Box constraint: Projection on Z2
    bFX = leds*X;
    tmp = bFX + U2;
    Z2 = max(tmp,0);
   
    
    
    % Scaled dual variable update
    if beta > 0
        res1 = dX - Z1;
    else
        res1 = zeros(size(Z1));
    end
    res2 = bFX - Z2;
    U1 = U1 + res1;
    U2 = U2 + res2;
    
    
    % Residual computation
    hist.prRes(i) = sqrt(norm(res1,'fro')^2 + norm(res2,'fro')^2);
    hist.dualRes(i) = hist.rho(i)*sqrt(norm(Z1-Z1minus,'fro')^2 + norm(Z2-Z2minus,'fro')^2);
    
    
    % Re-scale the parameter rho
    if hist.prRes(i) > inputs.mu*hist.dualRes(i) && inputs.rescaleRho == true
        hist.rho(i+1) = hist.rho(i)*inputs.tauIncr;
    end;
    if hist.dualRes(i) > inputs.mu*hist.prRes(i) && inputs.rescaleRho == true
        hist.rho(i+1) = hist.rho(i)/inputs.tauDecr;
    end;
    
    if inputs.verbose == true
        fprintf('ADMM iter %i (%f), primal res %e, dual res %e \n',i,toc(t1),hist.prRes(i),hist.dualRes(i));
        fprintf('     -> PCG (%f), err %f, nIter %i\n',t3,hist.pcg.acc(i),hist.pcg.iter(i));
    end
    
    if max(hist.prRes(i),hist.dualRes(i)) < inputs.tol
        break;
    end;
    
    if isempty(inputs.reference) == false
        hist.conv(i) = norm(wghtImg(:) - inputs.reference(:));
        if inputs.verbose == true
            fprintf('     -> True error %f\n',hist.conv(i));
        end
        
        if hist.conv(i) < inputs.tolConv, break; end;
    end
    
    % Now we need to re-scale the scaled dual variable U as well as the
    % dual residuals
    U1 = U1*hist.rho(i)/hist.rho(i+1);
    U2 = U2*hist.rho(i)/hist.rho(i+1);
    
    Z1minus = Z1;
    Z2minus = Z2;
    
end

% If we terminate earlier, remove the unused portions of the vectors.
hist.prRes = hist.prRes(1:i);
hist.dualRes = hist.dualRes(1:i);
hist.rho = hist.rho(1:i);
hist.pcg.Iter = hist.pcg.iter(1:i);
hist.pcg.acc = hist.pcg.acc(1:i);

if isempty(inputs.reference) == false
    hist.conv = hist.conv(1:i);
end

pred = zeros(size(ambient));
for i=1:nFilters
    pred(i,:) = sum(squeeze(ledMeasurements(i,:,:)).*X);
end

pred = reshape(pred',[h, w ,nFilters]);

X = reshape(X',[h w nChannels]);



end

%% Helper functions called by the least-squares conjugate gradient solver.

function res = applyAt( meas, spatialSmooth, nonnegativity, h, w, ledMeasurements, leds, beta, rho)

nChannels = size(leds,2);

tmp = repmat(meas,[1 1 nChannels]);
tmp = permute(tmp,[1 3 2]);
s1 = squeeze(sum(tmp.*ledMeasurements));

%{
s1 = zeros(size(ledMeasurements,2),size(meas,2));
for i=1:size(meas,2)
    s1(:,i) = ledMeasurements(:,:,i)'*meas(:,i);
end
%}
if beta > 0
    
    % Derivatives in the x direction
    spatialXSm = spatialSmooth(:,1:h*(w-1));
    spatialXSm = reshape(spatialXSm',[h, w-1, nChannels]);
    
    dx2Img = imfilter(spatialXSm,[-1 1 0]);
    dx2Img = dx2Img(:,2:w-1,:);
    dx2Img = cat(2,spatialXSm(:,1,:),dx2Img,-spatialXSm(:,w-1,:,:));
    
    dx2ImgMat = reshape(dx2Img,h*w,nChannels)';
    
    % derivatives in the y direction
    spatialYSm = spatialSmooth(:,h*(w-1)+1:h*(w-1) + w*(h-1));
    spatialYSm = reshape(spatialYSm', [h-1,w,nChannels]);
    
    dy2Img = imfilter(spatialYSm,[-1 1 0]');
    dy2Img = dy2Img(2:h-1,:,:);
    dy2Img = cat(1,spatialYSm(1,:,:),dy2Img,-spatialYSm(h-1,:,:));
    
    dy2ImgMat = reshape(dy2Img,h*w,nChannels)';
    
    s34 = rho*(dy2ImgMat + dx2ImgMat);
else
    s34 = 0;
end

s5 = rho*leds'*nonnegativity;

res = s1 + s34 + s5;
res = res(:);

end

function res = applyAtA( meas, h,w, ledMeasTledMeas, leds, alpha, beta, rho)

% meas will be passed as a vector
nChannels = size(leds,2);
nWaves = size(leds,1);

meas = reshape(meas,nChannels,h*w);
Rlambda = [eye(nWaves-1) zeros(nWaves-1,1)] - [zeros(nWaves-1,1) eye(nWaves-1)];

% Measurement approximation
%{
s1 = zeros(size(meas));
for i=1:size(meas,2);
    s1(:,i) = ledMeasurements(:,:,i)'*ledMeasurements(:,:,i)*meas(:,i);
end
%}

tmp = repmat(meas,[1 1 nChannels]);
tmp = permute(tmp,[1 3 2]);
s1 = squeeze(sum(tmp.*ledMeasTledMeas));

% Smoothness
if alpha > 0
    M2 = alpha*leds'*(Rlambda'*Rlambda)*leds;
else
    M2 = 0;
end
s2 = M2*meas;

% Now the spatial filtering
if beta > 0
    img = reshape(meas',h,w,nChannels);
    
    % Derivatives in the x direction
    firstLastColumn = imfilter(img(:,[1 2 w-1 w],:),[1 -1 0],'same');
    interior = imfilter(img,[-1 2 -1],'same');
    interior(:,1,:) = firstLastColumn(:,2,:);
    interior(:,w,:) = -firstLastColumn(:,4,:);
    
    dx2ImgMat = reshape(interior,h*w,nChannels)';
    
    % Derivatives in the y direction
    firstLastRow = imfilter(img([1 2 h-1 h],:,:), [1 -1 0]','same');
    interior = imfilter(img,[-1 2 -1]','same');
    interior(1,:,:) = firstLastRow(2,:,:);
    interior(h,:,:) = -firstLastRow(4,:,:);
    
    dy2ImgMat = reshape(interior,h*w,nChannels)';
    
    
    s34 = rho*(dx2ImgMat + dy2ImgMat);
else
    s34 = 0;
end

% Box constraint
M3 = rho*(leds'*leds);
s5 = M3*meas;


res = s1 + s2 + s34 + s5 ;
res = res(:);

end

function [grX, grY] = computeGradient( in )

    % Derivatives in the x direction
    grX = imfilter(in,[1 -1 0]);
    grX = grX(:,2:end,:);

    % derivatives in the y direction
    grY = imfilter(in,[1 -1 0]');
    grY = grY(2:end,:,:);

end


