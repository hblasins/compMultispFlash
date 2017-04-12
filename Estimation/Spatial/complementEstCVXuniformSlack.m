function [ X, slack, pred ] = complementEstCVXuniformSLack( illRef, ambientWghts, leds, cameraMat, alpha, beta, varargin)

p = inputParser;
p.addParameter('flashMode',true);
p.parse(varargin{:});
inputs = p.Results;


nChannels = size(leds,2);
nFilters = size(cameraMat,2);
nWaves = size(leds,1);
imHeight = size(ambientWghts,1);
imWidth = size(ambientWghts,2);

ambient = cameraMat'*leds*(reshape(ambientWghts,[imHeight*imWidth, nChannels])');
ambient = ambient(:);

appearanceRef = cameraMat'*illRef;
appearanceRef = repmat(shiftdim(appearanceRef,-2),[imHeight imWidth, 1]);
appearanceRef = appearanceRef(:);

A = cameraMat'*leds;

AA = [];
for i=1:imHeight*imWidth
    AA = blkdiag(AA,A);
end

Rlambda = [eye(nWaves-1) zeros(nWaves-1,1)] - [zeros(nWaves-1,1) eye(nWaves-1)];
Rrows = [eye(imHeight-1) zeros(imHeight-1,1)] - [zeros(imHeight-1,1) eye(imHeight-1)];
Rcols = [eye(imWidth-1) zeros(imWidth-1,1)] - [zeros(imWidth-1,1) eye(imWidth-1)];

cvx_begin
    variables estCvx(imHeight*imWidth*nChannels,1) slack(1,1)
    
    anisotropicGrad = 0;
    w = reshape(estCvx,[nChannels, imHeight, imWidth]);

    if beta > 0
        w1 = permute(w,[2, 3, 1]);
        w2 = permute(w,[3, 2, 1]);
    
    
        for i=1:nChannels
            anisotropicGrad = anisotropicGrad + beta*( sum(sum(abs(Rrows*w1(:,:,i)))) + sum(sum(abs(Rcols*w2(:,:,i)))) );
        end
    
    end
    
    
    minimize sum(sum_square(ambient(:) + AA*estCvx - slack*appearanceRef(:) )) + alpha*sum(sum_square(Rlambda*leds*reshape(estCvx,[nChannels,imHeight*imWidth]))) + anisotropicGrad
    
    subject to
        estCvx >= 0
        if inputs.flashMode
           estCvx <= 1 
        end
        
        slack >= 0
cvx_end

pred = [];
% pred = (ambient(:) + AA*estCvx - slack*appearanceRef(:))/slack;
% pred = reshape(pred, [nFilters, imHeight, imWidth]);
% pred = permute(pred,[2, 3, 1]);

estCvx = reshape(estCvx,[nChannels, imHeight, imWidth]);
X = permute(estCvx,[2, 3, 1]);
