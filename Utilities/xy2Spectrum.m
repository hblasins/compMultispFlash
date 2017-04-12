function [ spd ] = xy2Spectrum( x, y, wave )

% [ spd ] = xy2Spectrum( x, y, wave )
%
% Generate the smoothest spectral power distribution with particular CIE xy
% chromaticity coordinates and wavelength sampling. Returns a vector of
% NaNs if given chromaticity coordinates are not feasible.
%
% Copyright, Henryk Blasinski 2017

nWave = length(wave);

R = [eye(nWave-1) zeros(nWave-1,1)] - [zeros(nWave-1,1) eye(nWave-1)];
cie = ieReadSpectra(fullfile(isetRootPath,'data','human','XYZ'),wave);

% Solve a cvx problem to find the spectrum

cvx_begin quiet
    variable spd(nWave,1)
    minimize norm(R*spd,2)
    subject to
        cie(:,1)'*spd == x;
        cie(:,2)'*spd == y;
        cie(:,3)'*spd == 1-x-y;
        spd >= 0
cvx_end

% If the problem is infeasible return NaNs
if strcmp(cvx_status,'Infeasible')
    spd = NaN;
else
    spd = Quanta2Energy(wave,spd*10^10);
    spd = spd/max(spd);
end

spd = spd(:);

end

