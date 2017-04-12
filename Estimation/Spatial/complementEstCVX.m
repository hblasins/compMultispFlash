function [ X, slack ] = complementEstCVX( illRef, ambientWghts, leds, cameraMat, alpha, beta, varargin)

p = inputParser;
p.addParameter('flashMode',true);
p.parse(varargin{:});
inputs = p.Results;


nChannels = size(leds,2);
nWaves = size(leds,1);
imHeight = size(ambientWghts,1);
imWidth = size(ambientWghts,2);

ambient = cameraMat'*leds*reshape(ambientWghts,[imHeight*imWidth, nChannels])';
ambient = ambient(:);

A = [-cameraMat'*leds cameraMat'*illRef];

AA = [];
for i=1:imHeight*imWidth
    AA = blkdiag(AA,A);
end

Rlambda = [eye(nWaves-1) zeros(nWaves-1,1)] - [zeros(nWaves-1,1) eye(nWaves-1)];
Rrows = [eye(imHeight-1) zeros(imHeight-1,1)] - [zeros(imHeight-1,1) eye(imHeight-1)];
Rcols = [eye(imWidth-1) zeros(imWidth-1,1)] - [zeros(imWidth-1,1) eye(imWidth-1)];

cvx_begin
    variables estCvx(imHeight*imWidth*(nChannels+1),1)
    
    anisotropicGrad = 0;
    w = reshape(estCvx,[nChannels + 1, imHeight, imWidth]);

    if beta > 0
        w1 = permute(w,[2, 3, 1]);
        w2 = permute(w,[3, 2, 1]);
    
    
        for i=1:nChannels
            anisotropicGrad = anisotropicGrad + beta*( sum(sum(abs(Rrows*w1(:,:,i)))) + sum(sum(abs(Rcols*w2(:,:,i)))) );
        end
    
    end
    
    
    minimize sum(sum_square(ambient(:) - AA*estCvx)) + alpha*sum(sum_square(Rlambda*[leds, zeros(nWaves,1)]*reshape(estCvx,[nChannels+1,imHeight*imWidth]))) + anisotropicGrad
    
    
    subject to
        estCvx >= 0
        if inputs.flashMode
           w(1:nChannels,:,:) <= 1 
        end
cvx_end

estCvx = reshape(estCvx,[nChannels+1, imHeight, imWidth]);
X = permute(estCvx,[2, 3, 1]);

slack = X(:,:,end);
X = X(:,:,1:nChannels);
