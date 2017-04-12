function [ X, hist ] = ambientEstCVX( ambient, ledMeasurements, leds, alpha, beta )

imWidth = size(ledMeasurements,2);
imHeight = size(ledMeasurements,1);

nWaves = size(leds,1);
nChannels = size(leds,2);

% Prepare matrices
mVec = permute(ambient,[3 1 2]);
mVec = mVec(:);

LL = [];
for ww=1:imWidth
    for hh=1:imHeight
        subL = squeeze(ledMeasurements(hh,ww,:,:));
        LL = blkdiag(LL,subL);
    end
end

Rlambda = [eye(nWaves-1) zeros(nWaves-1,1)] - [zeros(nWaves-1,1) eye(nWaves-1)];

Rrows = [eye(imHeight-1) zeros(imHeight-1,1)] - [zeros(imHeight-1,1) eye(imHeight-1)];
Rcols = [eye(imWidth-1) zeros(imWidth-1,1)] - [zeros(imWidth-1,1) eye(imWidth-1)];

cvx_begin
    variables wEstCvx(nChannels*imHeight*imWidth,1)
    
    anisotropicGrad = 0;
    if beta > 0
        w = reshape(wEstCvx,[nChannels, imHeight, imWidth]);
        w1 = permute(w,[2, 3, 1]);
        w2 = permute(w,[3, 2, 1]);
    
        for i=1:nChannels
            anisotropicGrad = anisotropicGrad + beta*( sum(sum(abs(Rrows*w1(:,:,i)))) + sum(sum(abs(Rcols*w2(:,:,i)))) );
        end
    end
    
    minimize sum(square(mVec - LL*wEstCvx)) + alpha*sum(sum(square(Rlambda*leds*reshape(wEstCvx,[nChannels, imHeight*imWidth])))) + anisotropicGrad
    subject to
        leds*reshape(wEstCvx,[nChannels, imHeight*imWidth]) >= 0
cvx_end

X = reshape(wEstCvx,[nChannels, imHeight, imWidth]);
X = permute(X,[2 3 1]);
