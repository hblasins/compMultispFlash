function [ measurement, reference ] = getRenderedData( sensor, wave, surfaceSpectrum, ledSpectra, varargin)

p = inputParser;
p.addRequired('sensor')
p.addRequired('wave');
p.addRequired('surfaceSpectrum');
p.addRequired('ledSpectra');
p.addOptional('targetDistance',1000);
p.addOptional('depth',10000);
p.addOptional('chlConc',0);
p.addOptional('cdomConc',0);
p.addOptional('smallPartConc',0);
p.addOptional('largePartConc',0);
p.addOptional('target','Macbeth');
p.addOptional('rtbResultFolder',fullfile('/','home','hblasins','Documents','MATLAB','render_toolbox'));

p.parse(sensor, wave, surfaceSpectrum, ledSpectra, varargin{:});

nLeds = size(ledSpectra,2);


% Read in the oi data
[oi, referenceOi] = getOiData(varargin{:});
  

% Simulate sensor
[ measurement.sensors, measurement.ips, measurement.raw.data, measurement.demosaiced.data, measurement.patch.data ] = getSensorData( sensor, oi);
[ reference.sensors, reference.ips, reference.raw.data, reference.demosaiced.data, reference.patch.data ] = getSensorData( sensor, referenceOi);


measurement.demosaiced.data = measurement.demosaiced.data/max(measurement.demosaiced.data(:));
measurement.patch.data = measurement.patch.data/max(measurement.patch.data(:));
measurement.raw.data = measurement.raw.data/max(measurement.raw.data(:));

measurement.demosaiced.ambient = measurement.demosaiced.data(:,:,:,1);
measurement.demosaiced.led = max(measurement.demosaiced.data(:,:,:,2:end) - repmat(measurement.demosaiced.ambient,[1 1 1, nLeds]),0);

measurement.raw.ambient = measurement.raw.data(:,:,1);
measurement.raw.led = max(measurement.raw.data(:,:,2:end) - repmat(measurement.raw.ambient,[1 1 nLeds]),0);

measurement.patch.ambient = measurement.patch.data(:,1,:);
measurement.patch.led = max(measurement.patch.data(:,2:end,:) - repmat(measurement.patch.ambient,[1 nLeds,1]),0);


measurement.vectorized.ambient = reshape(measurement.demosaiced.ambient,[sensorGet(sensor,'rows')*sensorGet(sensor,'cols'), 3])';
measurement.vectorized.led = permute(reshape(measurement.demosaiced.led,[sensorGet(sensor,'rows')*sensorGet(sensor,'cols'), 3, nLeds]),[2, 3, 1]);

end

