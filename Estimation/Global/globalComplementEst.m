function [ flashEst, flashWghts ] = globalComplementEst( desiredIll, ambientEstimate, leds, cameraResp, varargin)

p = inputParser;
p.addRequired('desiredIll');
p.addRequired('ambientWghts');
p.addRequired('leds');
p.addRequired('cameraResp');
p.addOptional('flashMode',true);

p.parse(desiredIll, ambientEstimate, leds, cameraResp, varargin{:});
inputs = p.Results;

nChannels = size(inputs.leds,2);

camDes = inputs.cameraResp'*desiredIll;
camDes = camDes/max(camDes);

cvx_begin
    variables flashCompWeights(nChannels,1) en(1,1)
    minimize norm( inputs.cameraResp'*(inputs.leds*(flashCompWeights) + ambientEstimate)- en*camDes)
    subject to
       if inputs.flashMode
           1>= flashCompWeights >= 0
       else
           flashCompWeights >= 0
       end
       en >= 0
cvx_end

flashEst = inputs.leds*flashCompWeights;
flashWghts = flashCompWeights;


end

