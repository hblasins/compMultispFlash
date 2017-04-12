function [ sensor, ip ] = renderFlashImage( baseImage, flashImages, weights, sensor, varargin )

% [ sensor, ip ] = renderFlashImage( baseImage, flashImages, weights, sensor, ... )
%
% Render the computational flash image given a set of LED only images,
% computed weights and the base image. Returs ISET sensor and ip structures
% containing the computational image data.
%
% Optional arguments:
%    ev - a double indicating the fraction of pixels in the image that will
%    be saturated during histogram streching (Default = 0.01)
%    name - a string used in the sensor and ip name.
%
% Copyright, Henryk Blasinski 2017

p = inputParser;
p.addOptional('ev',0.01);
p.addOptional('name','');
p.parse(varargin{:});
inputs = p.Results;

for i=1:length(weights)
    baseImage = baseImage + weights(i)*flashImages(:,:,i);
end

baseImage = imageExpose(baseImage,inputs.ev);

sensor = sensorSet(sensor,'volts',baseImage);
sensor = sensorSet(sensor,'name',inputs.name);
ieAddObject(sensor);
sensorWindow;


ip = ipCreate;
% ip = ipSet(ip,'sensor conversion method','none');
ip = ipCompute(ip,sensor);
ip = ipSet(ip,'name',inputs.name);
ieAddObject(ip);

ipWindow;

end

