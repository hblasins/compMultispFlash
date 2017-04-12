function [ result ] = imageExpose( img, percent )

% [ result ] = imageExpose( img, percent )
%
% Perform image histogram streching so that percent of pixels are
% saturated.
%
% Copyright, Henryk Blaisnski 2017


[h, bins] = hist(img(:),256);
chist= cumsum(h)/sum(h);

loc = find(chist >= (1-percent),1,'first');

nF = bins(loc);

result = img/nF;

end

