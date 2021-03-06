function [ X, slack, hist ] = complementEstADMMuniformSlack( illRef, ambientWghts, leds, cameraMat, alpha, beta, varargin )


p = inputParser;
p.addParameter('maxIter',100,@isscalar);
p.addParameter('tol',1e-3,@isscalar);
p.addParameter('mu',10,@isscalar);
p.addParameter('tauIncr',5,@isscalar);
p.addParameter('tauDecr',5,@isscalar);
p.addParameter('rhoInit',1,@isscalar);
p.addParameter('rescaleRho',false,@islogical);
p.addParameter('reference',[]);
p.addParameter('tolConv',0);
p.addParameter('verbose',false);
p.addParameter('gradient','anisotropic');
p.addParameter('flashMode',true);
p.parse(varargin{:});
inputs = p.Results;

nChannels = size(leds,2);
nWaves = size(leds,1);
h = size(ambientWghts,1);
w = size(ambientWghts,2);

illRefCam = cameraMat'*illRef;

ambientWghts = reshape(ambientWghts,h*w,nChannels)';
ambientCam = cameraMat'*leds*ambientWghts;

Z1 = zeros(nChannels,(h-1)*w + (w-1)*h);
Z1minus = zeros(size(Z1));
U1 = zeros(nChannels,(h-1)*w + (w-1)*h);

Z2 = zeros(nChannels,h*w);
Z2minus = zeros(size(Z2));
U2 = zeros(nChannels,h*w);

Z2slack = 0;
Z2slackminus = 0;
U2slack = 0;

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
    Atb = applyAt(ambientCam, (Z1-U1), (Z2-U2), (Z2slack-U2slack),h,w,illRefCam,cameraMat'*leds, alpha, beta, hist.rho(i)/2);
    AtAhndl = @(x) applyAtA(x, h, w, illRefCam, cameraMat, leds, alpha, beta, hist.rho(i)/2);
    [X, ~, hist.pcg.acc(i), hist.pcg.iter(i)] = pcg(AtAhndl,Atb,1e-8,10000,[],[],X(:));
    t3 = toc(t2);
    
   
    
    %% Test 
    %{
    AA = [];
    for j=1:h*w
        AA = blkdiag(AA,cameraMat'*leds);
    end
    
    AA = [AA, -repmat(illRefCam,h*w,1)];
    bb = ambientCam(:);
    
    cvx_begin
        variables xt(h*w*nChannels+1,1)
        minimize sum_square(AA*xt + bb) + hist.rho(i)/2*sum_square(xt(1:h*w*nChannels) - (Z2(:) - U2(:))) + hist.rho(i)/2*sum_square(xt(end) - (Z2slack - U2slack));
    cvx_end
      
    [X, xt]
    
    norm(X-xt)
    
    X = xt;
    
    %}
   
    
    
    
    
    slack = X(end);
    
    wghts = reshape(X(1:end-1),nChannels,h*w);
    wghtImg = reshape(wghts',[h,w,nChannels]);

    
    
    
    % Anisotropic TV: Soft thresholding on Z1
    % We apply smoothes to weight estimates only and we ignore slack
    % variables
    if beta ~= 0
        
        [grX, grY] = computeGradient(wghtImg(:,:,1:nChannels));
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
        % Z1 = zeros(size(Z1));
        Z1 = 0;
    end
    
    % Box constraint: Projection on Z2. If in flash mode then weights are
    % constrained from above by 1 and still have to be non-negative.
    % We allow slack variable to be greater than 1.
    tmp = wghts + U2;
    if inputs.flashMode
        tmp = min(tmp,1);
    end
    Z2 = max(tmp,0);
    res2 = wghts - Z2;
    
    
    tmp = slack + U2slack;
    Z2slack = max(tmp,0);
    res2slack = slack - Z2slack;
    
    
    % Scaled dual variable update
    if beta ~= 0
        res1 = dX - Z1;
    else
        res1 = 0;
    end
    
    U1 = U1 + res1;
    U2 = U2 + res2;
    U2slack = U2slack + res2slack;
    
    
    % Residual computation
    hist.prRes(i) = sqrt(norm(res1,'fro')^2 + norm(res2,'fro')^2 + norm(res2slack,'fro')^2);
    hist.dualRes(i) = hist.rho(i)*sqrt(norm(Z1-Z1minus,'fro')^2 + norm(Z2-Z2minus,'fro')^2 + norm(Z2slack-Z2slackminus,'fro')^2);
    
    
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
    U2slack = U2slack*hist.rho(i)/hist.rho(i+1);
    
    Z1minus = Z1;
    Z2minus = Z2;
    Z2slackminus = Z2slack;
    
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

slack = X(end);

X = reshape(X(1:end-1),nChannels,h*w);
if inputs.flashMode
     X = min(max(reshape(X',[h,w,nChannels]),0),1);
else
     X = max(reshape(X',[h,w,nChannels]),0);
end



end

%% Helper functions called by the least-squares conjugate gradient solver.

function res = applyAt( ambientCam, spatialSmooth, nonnegativity, slack, h, w, illRefCam, cameraMatTleds, alpha, beta, rho)

nChannels = size(cameraMatTleds,2);

% individual
%{
A = [-cameraMatTleds, illRefCam];
s1 = A'*ambientCam;
%}
% shared

%% This contains slack
s1 = -cameraMatTleds'*ambientCam;
sd = illRefCam'*sum(ambientCam,2) + rho*slack;

% Removing slac variable
% s1 = -cameraMatTleds'*ambientCam;
% sd = [];

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

% Box constraint
s5 = rho*nonnegativity;


res = s1 + s34 + s5;
res = [res(:); sd ];

end

function res = applyAtA( meas, h, w, illRefCam, cameraMat, leds, alpha, beta, rho)

% meas will be passed as a vector
nChannels = size(leds,2);
nWaves = size(leds,1);

% With slack
slack = meas(end);
meas = reshape(meas(1:end-1),nChannels ,h*w);



cameraMatTledsTillRefCam = -(cameraMat'*leds)'*illRefCam;
sc = slack*repmat(cameraMatTledsTillRefCam,[1 h*w]);
% sc = 0;

A = -cameraMat'*leds;
M1 = A'*A;

newSlack = -illRefCam'*(cameraMat'*leds)*sum(meas,2) + h*w*(illRefCam'*illRefCam)*slack + rho*slack;
% newSlack = [];

% Smoothness
if alpha > 0
    Rlambda = [eye(nWaves-1) zeros(nWaves-1,1)] - [zeros(nWaves-1,1) eye(nWaves-1)];
    
    M2 = alpha*leds'*(Rlambda'*Rlambda)*leds;
else
    M2 = 0;
end
s12 = (M1 + M2)*meas + sc;

% Now the spatial filtering
if beta > 0
    tmpMeas = meas(1:nChannels,:);
    img = reshape(tmpMeas',h,w,nChannels);
    
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
% We add a zero row and column to account for the slack variable.
% We care abot nonnegativity of weight coefficients, not the spectrum.
s5 = rho*meas;


res = s12 + s34 + s5 ;
res = [res(:); newSlack];

end

function [grX, grY] = computeGradient( in )

    % Derivatives in the x direction
    grX = imfilter(in,[1 -1 0]);
    grX = grX(:,2:end,:);

    % derivatives in the y direction
    grY = imfilter(in,[1 -1 0]');
    grY = grY(2:end,:,:);

end

