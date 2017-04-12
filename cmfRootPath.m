function rootPath=cmfRootPath()

% function rootPath=compMultispFlashRootPath()
%
% Returns the absolute path for the directory containing the computational
% multispectral flash code.
%
% Copyright, Henryk Blasinski 2017

rootPath=which('cmfRootPath');

[rootPath,fName,ext]=fileparts(rootPath);

return
