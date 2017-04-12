function [ ambientEst, ambientWghts, ambientPredictions ] = globalAmbientEst( ambient, ledMeasurements, leds, varargin )

p = inputParser;
p.addRequired('ambient');
p.addRequired('ledMeasurements');
p.addRequired('leds');
p.addOptional('alpha',1);

p.parse(ambient,ledMeasurements,leds, varargin{:});
inputs = p.Results;

nWaves = size(inputs.leds,1);
nChannels = size(inputs.leds,2);

R = [eye(nWaves-1) zeros(nWaves-1,1)] - [zeros(nWaves-1,1) eye(nWaves-1)];

cvx_begin
    variables ambientApproxWghts(nChannels,1)

    approx = 0;
    for i=1:nChannels
        approx = approx + squeeze(ledMeasurements(:,i,:))*ambientApproxWghts(i);
    end

    minimize sum(norms(squeeze(ambient) - approx,2,1)) + inputs.alpha * norm(R*inputs.leds*ambientApproxWghts,2)
    subject to
        inputs.leds*ambientApproxWghts >= 0
        % ambientApproxWghts >= 0
cvx_end

ambientWghts = ambientApproxWghts;
ambientEst = inputs.leds*ambientApproxWghts;
ambientPredictions = approx;

end

